#include <sycl/sycl.hpp>
#include "numeric"
#include <vector>
#include <tuple>
#include <iterator>
#include <type_traits>
#include <algorithm>
#include <fullerenes/sycl-headers/execution-compat.hh>
#include <numeric>
#include <fullerenes/kernel-headers/spherical-projection-functor.hh>
#include <fullerenes/kernel-headers/sycl-parallel-primitives.hh>
#include "kernel.cc"
#include "forcefield-includes.cc"

using namespace sycl;
template <typename T>
struct CuDeque
{
private:

    int front, back, q_size, capacity;
    local_accessor<T, 1> array;

public:
    //This is a custom implementation of a deque class that is meant to be used on the device. 
//The deque is implemented using a circular buffer in a shared memory array. 
//The class is used to store the work queue of the thread blocks, and the class provides the
//necessary functions to allow the work queue to be used as a work stealing queue. 
//The class is implemented as a template class to allow the user to specify the type of data that the
//deque will hold.
CuDeque(const local_accessor<T,1> memory, const int capacity): array(memory), front(-1), back(0), q_size(0), capacity(capacity) {}
    
    /**
     * @brief  Returns the size of the queue. 
     * @param  None
     * @retval Size of the queue.
     */
    int size(){return q_size;}

    /**
     * @brief  Returns true if the queue is empty and false otherwise. 
     * @param  None
     * @retval True if the queue is empty and false otherwise.
     */
    bool empty(){
        return (front == -1);
    }

    /**
     * @brief  Returns true if the queue is full and false otherwise. 
     * @param  None
     * @retval True if the queue is full and false otherwise.
     */
    bool full(){
        return (front == 0 && back == capacity-1) || (front == back+1);
    }
    
    /**
     * @brief  This function is used to pop the first element of the queue.
     * @param  None
     * @retval First element of the queue
     */
    T pop_front(){
        if (empty()){ assert(false); return T();} 
        T return_val = array[front];
        if(front == back) {
            front = -1;
            back = -1;
        } 
        else if (front == capacity-1) front = 0;
        else front = front+1;
        q_size--;
        return return_val;
    }

    /**
     * @brief Returns the last element of the queue and removes it from the queue.
     * @param None
     * @return The last element of the queue
     */
    T pop_back(){
        if (empty()){ assert(false); return T();}
        T return_val = array[back];
        if(front == back) {
            front = -1;
            back = -1;
        } 
        else if (back == 0) back = capacity-1;
        else back = back-1;
        q_size--;
        return return_val;
    }

    /** @brief Insert a value into the back of the queue 
     *  @param val the value to insert
     */
    void push_back(T val){
        assert(!full());
        if (front == -1) {
            front = 0;
            back = 0;
        }
        else if (back == capacity-1) back = 0;
        else back = back+1;
        array[back] = val;
        q_size++;
    }

    /** @brief Insert a value into the front of the queue
     *  @param val the value to insert
     */
    void push_front(T val){
        assert(!full());
        if (front == -1) {
            front = 0;
            back = 0;
        }
        else if (front == 0) front = capacity-1;
        else front = front-1;
        array[front] = val;
        q_size++;
    }
};

template <typename K>
K multiple_source_shortest_paths(const sycl::group<1>& cta, const Span<std::array<K,3>> cubic_neighbours,const local_accessor<int,1>& distances, const local_accessor<K,1>& smem){
    INT_TYPEDEFS(K);
    constexpr int unreachable_distance = std::numeric_limits<int>::max();
    auto N = cta.get_local_linear_range();
    auto tid = cta.get_local_linear_id();
    DeviceCubicGraph FG(cubic_neighbours);
    std::array<K,6> outer_face; memset(outer_face.data(), 0, 6*sizeof(node_t));
    uint8_t Nface = FG.get_face_oriented(0,FG[0][0], outer_face);
    distances[tid] = unreachable_distance;
    sycl::group_barrier(cta);
    if(tid < Nface) distances[outer_face[tid]] = 0;
    sycl::group_barrier(cta);
    if(tid == 0){
        CuDeque<K> work_queue(smem, N);
        for (size_t i = 0; i < Nface; i++) work_queue.push_back(outer_face[i]);
        while(!work_queue.empty()){
            auto v = work_queue.pop_front();
            for (size_t i = 0; i < 3; i++){
                auto w = FG[v][i];
                if(distances[w] == unreachable_distance){
                    distances[w] = distances[v] + 1;
                    work_queue.push_back(w);
                }
            }
        }
    }
    sycl::group_barrier(cta);
    node_t distance = distances[tid];
    sycl::group_barrier(cta);
    return distance;
}

template<typename T, typename K>
SyclEvent spherical_projection(SyclQueue& Q, FullereneBatchView<T,K>& batch){
    TEMPLATE_TYPEDEFS(T,K);
    constexpr real_t scalerad = 4.0;
    SyclEventImpl projection_done = Q->submit([&](handler& h) {
        auto N = batch.N_;
        auto capacity = batch.capacity();

        local_accessor<node_t, 1>   work_queue_memory(N*2, h);
        local_accessor<int, 1>      smem(N, h); //Has to be int for atomic operations
        local_accessor<coord2d, 1>  atomic_coordinate_memory(N, h);
        local_accessor<coord3d, 1>  xyz_smem(N, h);


        h.parallel_for(sycl::nd_range(sycl::range(N*capacity), sycl::range(N)), [=](nd_item<1> nditem) {
            auto cta = nditem.get_group();
            auto tid = nditem.get_local_linear_id();
            auto isomer_idx = nditem.get_group_linear_id();
            Fullerene<T,K> full = batch[isomer_idx];
            auto cubic_neighbours = full.d_.A_cubic_;
            auto xys_acc = full.d_.X_cubic_.template as_span<std::array<T,2>>();
            auto xyz_acc = full.d_.X_cubic_.template as_span<std::array<T,3>>();

            if(isomer_idx >= capacity) assert(false);
            if ( full.m_.flags_.get().is_set(StatusEnum::FULLERENEGRAPH_PREPARED)){
            atomic_coordinate_memory[tid] = {0.0, 0.0};
            NodeNeighbours node_graph(cubic_neighbours, (K)tid);
            node3 neighbours = node_graph.cubic_neighbours;
            node_t distance = multiple_source_shortest_paths(cta, cubic_neighbours, smem, work_queue_memory);
            node_t d_max = reduce_over_group(cta, distance, maximum<node_t>{});
            smem[tid] = 0;
            sycl::group_barrier(cta);
            sycl::atomic_ref<int, sycl::memory_order::seq_cst, sycl::memory_scope::work_group> atomic_same_dist(smem[distance]);
            atomic_same_dist.fetch_add(1);
            sycl::group_barrier(cta);
            node_t num_same_dist = smem[distance];
            sycl::group_barrier(cta);
            coord2d xys = xys_acc[tid];
            sycl::atomic_ref<real_t, sycl::memory_order::seq_cst, sycl::memory_scope::work_group> atomic_coord_x(atomic_coordinate_memory[distance][0]);
            sycl::atomic_ref<real_t, sycl::memory_order::seq_cst, sycl::memory_scope::work_group> atomic_coord_y(atomic_coordinate_memory[distance][1]);
            atomic_coord_x.fetch_add(xys[0]);
            atomic_coord_y.fetch_add(xys[1]);
            sycl::group_barrier(cta);

            coord2d centroid = atomic_coordinate_memory[distance]/real_t(num_same_dist);
            coord2d xy = xys - centroid;
            real_t dtheta = real_t(M_PI)/real_t(d_max+1);
            real_t phi = dtheta*(distance+0.5);
            real_t theta = sycl::atan2(xy[0],xy[1]);
            coord3d xyz = {sycl::cos(theta)*sycl::sin(phi), sycl::sin(theta)*sycl::sin(phi), sycl::cos(phi)};
            real_t xsum = sycl::reduce_over_group(cta, xyz[0], sycl::plus<real_t>{});
            real_t ysum = sycl::reduce_over_group(cta, xyz[1], sycl::plus<real_t>{});
            real_t zsum = sycl::reduce_over_group(cta, xyz[2], sycl::plus<real_t>{});
            coord3d cm = {xsum/real_t(N), ysum/real_t(N), zsum/real_t(N)};

            xyz -= cm;
            real_t Ravg = 0.0;
            xyz_smem[tid] = xyz;
            sycl::group_barrier(cta);
            real_t local_Ravg = 0.0;
            for (size_t i = 0; i < 3; i++){ local_Ravg += norm(xyz_smem[tid] - xyz_smem[neighbours[i]]); }
            Ravg = sycl::reduce_over_group(cta, local_Ravg, sycl::plus<real_t>{})/real_t(3*N);
            xyz *= scalerad*real_t(1.5)/Ravg;
            xyz_acc[tid] = xyz;
            full.m_.flags_.get().set(StatusEnum::NOT_CONVERGED);
            } 
        });
    });
    return SyclEvent(std::move(projection_done));
}

template <typename K>
void multiple_source_shortest_paths(const Span<std::array<K,3>> neighbours, const std::vector<K>& sources, Span<K> distances, const unsigned int max_depth = INT_MAX)
{
    std::vector<K> queue_buf(neighbours.size());
    Deque<K> queue(queue_buf);
        
    for(K s: sources){
        distances[s] = 0;
        queue.push_back(s);
    }

    while(!queue.empty()){
        K v = queue.pop_front();
        for(K w: neighbours[v]){
            const edge_t edge(v,w);
            if(distances[w] == std::numeric_limits<K>::max()){ // Node not previously visited
                distances[w] = distances[v] + 1;
                if(distances[w] < max_depth) queue.push_back(w);
            }
        }
    }
}

template<typename T, typename K>
SyclEvent spherical_projection_impl( SyclQueue& Q,
                                Span<std::array<T,2>> xys,
                                Span<std::array<T,3>> X,
                                Span<std::array<K,3>> cubic_neighbours,
                                Span<K> distances,
                                Span<K> reduce_in,
                                Span<K> reduce_out,
                                Span<K> output_keys,
                                Span<std::array<T,2>> sorted_xys)
{
    //MSSPs
    auto N = X.size();
    Q.wait();
    std::fill(FULLERENE_PAR_UNSEQ distances.begin(), distances.end(), std::numeric_limits<K>::max());
    DeviceCubicGraph FG(cubic_neighbours);
    std::array<K,6> outer_face;
    auto face_size = FG.get_face_oriented(0, FG[0][0], outer_face);

    multiple_source_shortest_paths(cubic_neighbours, std::vector<K>(outer_face.data(), outer_face.data() + face_size), distances);

    //Compute maximum topological distance
    Q.wait();
    K d_max = std::reduce(FULLERENE_PAR_UNSEQ distances.begin(), distances.end(), K{0}, [](K a, K b){ return std::max(a, b); });
    //Count number of nodes at each distance
    Q.wait();
    std::copy(FULLERENE_PAR_UNSEQ xys.subspan(0,N).begin(), xys.subspan(0,N).end(), sorted_xys.begin());
    std::iota(reduce_in.subspan(0, N).begin(), reduce_in.subspan(0, N).end(), K{0});
    std::sort(FULLERENE_PAR_UNSEQ reduce_in.subspan(0, N).begin(), reduce_in.subspan(0, N).end(), [distances](K a, K b){return distances[a] < distances[b];});
    std::transform(FULLERENE_PAR_UNSEQ reduce_in.subspan(0, N).begin(), reduce_in.subspan(0, N).end(), sorted_xys.begin(), [xys](K idx){return xys[idx];});
    std::transform(FULLERENE_PAR_UNSEQ reduce_in.subspan(0, N).begin(), reduce_in.subspan(0, N).end(), reduce_out.begin(), [distances](K idx){return distances[idx];});
    auto summed_coordinates = reduce_in.template as_span<std::array<T,2>>().subspan(0, d_max + 1);
    primitives::reduce_by_segment(Q, reduce_out.subspan(0,N), sorted_xys, output_keys, summed_coordinates, std::equal_to<K>{}, [](std::array<T,2> a, std::array<T,2> b){return a + b;});
    Q.wait();
    std::fill(FULLERENE_PAR_UNSEQ sorted_xys.template as_span<K>().subspan(0, N).begin(), sorted_xys.template as_span<K>().subspan(0, N).end(), K{1});
    //Compute number of nodes at each distance and store in sorted_xys at indices N to  N + d_max + 1
    auto num_nodes_at_distance = sorted_xys.subspan(N, d_max + 1).template as_span<K>().subspan(0, d_max + 1);
    
    primitives::reduce_by_segment(Q, reduce_out.subspan(0,N), sorted_xys.template as_span<K>(), output_keys, num_nodes_at_distance);
    //Compute the centroid of the nodes at each distance
    auto centroids = reduce_out.template as_span<std::array<T,2>>().subspan(0, d_max + 1);
    Q.wait();
    std::transform(FULLERENE_PAR_UNSEQ summed_coordinates.begin(), summed_coordinates.end(), num_nodes_at_distance.begin(), centroids.begin(), [](std::array<T,2> a, K b){return a/T(b);});
    //Shift the coordinates of the nodes at each distance by the centroid
    Q.wait();
    std::transform(FULLERENE_PAR_UNSEQ xys.subspan(0,N).begin(), xys.subspan(0,N).end(), distances.begin(), sorted_xys.begin(), [centroids](std::array<T,2> a, K b){return a - centroids[b];});
    //Compute the spherical coordinates of the nodes at each distance
    Q.wait();
    std::transform(FULLERENE_PAR_UNSEQ sorted_xys.subspan(0,N).begin(), sorted_xys.subspan(0,N).end(), distances.begin(), X.begin(), [d_max](std::array<T,2> xy, K dist){
        T dtheta = T(M_PI)/T(d_max+1);
        T phi = dtheta*(dist+0.5);
        T theta = sycl::atan2(xy[0],xy[1]);
        return std::array<T,3>{sycl::cos(theta)*sycl::sin(phi), sycl::sin(theta)*sycl::sin(phi), sycl::cos(phi)};
    });
    //Compute center of mass
    Q.wait();
    std::array<T,3> cm = std::reduce(FULLERENE_PAR_UNSEQ X.begin(), X.end(), std::array<T,3>{0.0, 0.0, 0.0}, [](std::array<T,3> a, std::array<T,3> b){return a + b;});
    cm /= T(N);
    //Shift the coordinates by the center of mass
    std::transform(FULLERENE_PAR_UNSEQ X.begin(), X.end(), X.begin(), [cm](std::array<T,3> a){return a - cm;});
    //Compute the average distance between nodes
    T Ravg = std::transform_reduce(FULLERENE_PAR_UNSEQ X.begin(), X.end(), cubic_neighbours.template as_span<std::array<K,3>>().begin(), T{0.0}, std::plus<T>{}, [X](std::array<T,3> a, std::array<K,3> neighbours){
        T local_Ravg = 0.0;
        for (size_t i = 0; i < 3; i++){ local_Ravg += norm(a - X[neighbours[i]]); }
        return local_Ravg;
    }) / T(3*N);
    //Scale the coordinates
    T scalerad = 4.0;
    std::transform(FULLERENE_PAR_UNSEQ X.begin(), X.end(), X.begin(), [scalerad, Ravg](std::array<T,3> a){return a*scalerad*T(1.5)/Ravg;});
    auto ret_event = Q.get_event();
    return ret_event;
}

template <typename T, typename K>
SyclEvent SphericalProjectionFunctor<T,K>::compute(SyclQueue& Q, Fullerene<T,K> fullerene
                                                    , Span<K> topological_distances,
                                                    Span<K> reduce_in,
                                                    Span<K> reduce_out,
                                                    Span<K> output_keys,
                                                    Span<std::array<T,2>> sorted_xys){
    if (fullerene.m_.flags_.get() & (int)StatusEnum::FULLERENEGRAPH_PREPARED){
        auto ret_event = spherical_projection_impl<T,K>(Q,
                                    fullerene.d_.X_cubic_.template as_span<std::array<T,2>>(),
                                    fullerene.d_.X_cubic_,
                                    fullerene.d_.A_cubic_,
                                    topological_distances,
                                    reduce_in,
                                    reduce_out,
                                    output_keys,
                                    sorted_xys);
        fullerene.m_.flags_.get().set(StatusEnum::NOT_CONVERGED);
        return ret_event;
    } else return Q.get_event();
}

template <typename T, typename K>
SyclEvent SphericalProjectionFunctor<T,K>::compute(SyclQueue& Q, FullereneBatchView<T,K> batch){
    return spherical_projection<T,K>(Q, batch);
}

template struct SphericalProjectionFunctor<float,uint16_t>;
template struct SphericalProjectionFunctor<float,uint32_t>;
template struct SphericalProjectionFunctor<double,uint16_t>;
template struct SphericalProjectionFunctor<double,uint32_t>;