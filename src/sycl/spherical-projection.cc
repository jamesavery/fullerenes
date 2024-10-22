#include <CL/sycl.hpp>
#include "numeric"
#include <vector>
#include <tuple>
#include <iterator>
#include <type_traits>
#include <fullerenes/sycl-headers/spherical-projection-kernel.hh>
#include "kernel.cc"
#include "primitives.cc"
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
    auto N = cta.get_local_linear_range();
    auto tid = cta.get_local_linear_id();
    auto isomer_idx = cta.get_group_linear_id();
    DeviceCubicGraph FG(cubic_neighbours);
    std::array<K,6> outer_face; memset(outer_face.data(), 0, 6*sizeof(node_t));
    uint8_t Nface = FG.get_face_oriented(0,FG[0][0], outer_face);
    distances[tid] = std::numeric_limits<K>::max();
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
                if(distances[w] == std::numeric_limits<K>::max()){
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
        auto Nf = batch.Nf_;
        auto capacity = batch.capacity();
        auto max_iter = N * 10;

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
            sycl::atomic_ref<int, sycl::memory_order::relaxed, sycl::memory_scope::work_group> atomic_same_dist(smem[distance]);
            atomic_same_dist.fetch_add(1);
            sycl::group_barrier(cta);
            node_t num_same_dist = smem[distance];
            sycl::group_barrier(cta);
            coord2d xys = xys_acc[tid];
            sycl::atomic_ref<real_t, sycl::memory_order::relaxed, sycl::memory_scope::work_group> atomic_coord_x(atomic_coordinate_memory[distance][0]);
            sycl::atomic_ref<real_t, sycl::memory_order::relaxed, sycl::memory_scope::work_group> atomic_coord_y(atomic_coordinate_memory[distance][1]);
            atomic_coord_x.fetch_add(xys[0]);
            atomic_coord_y.fetch_add(xys[1]);
            sycl::group_barrier(cta);

            coord2d centroid = atomic_coordinate_memory[distance]/real_t(num_same_dist);
            coord2d xy = xys - centroid;
            real_t dtheta = real_t(M_PI)/real_t(d_max+1);
            real_t phi = dtheta*(distance+0.5);
            real_t theta = sycl::atan2(xy[0],xy[1]);
            coord2d spherical_coords = {theta, phi};
            coord3d xyz = {cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)};
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
    Deque<K> queue(neighbours.size());
        
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
    primitives::fill(Q, distances, std::numeric_limits<K>::max());
    DeviceCubicGraph FG(cubic_neighbours);
    std::array<K,6> outer_face;
    auto face_size = FG.get_face_oriented(0, FG[0][0], outer_face);

    multiple_source_shortest_paths(cubic_neighbours, std::vector<K>(outer_face.data(), outer_face.data() + face_size), distances);

    //Compute maximum topological distance
    K d_max = primitives::reduce(Q, distances, K{0}, Max{});
    //Count number of nodes at each distance
    primitives::copy(Q, xys.subspan(0,N), sorted_xys);
    primitives::iota(Q, reduce_in.subspan(0, N), K{0});
    primitives::sort(Q, reduce_in.subspan(0, N), [distances](K a, K b){return distances[a] < distances[b];}); 
    primitives::transform(Q, reduce_in.subspan(0, N), sorted_xys, [xys](K idx){return xys[idx];});
    primitives::transform(Q, reduce_in.subspan(0, N), reduce_out, [distances](K idx){return distances[idx];});
    auto summed_coordinates = reduce_in.template as_span<std::array<T,2>>().subspan(0, d_max + 1);
    primitives::reduce_by_segment(Q, reduce_out.subspan(0,N), sorted_xys, output_keys, summed_coordinates, Equal{}, [](std::array<T,2> a, std::array<T,2> b){return a + b;});
    primitives::fill(Q, sorted_xys.template as_span<K>().subspan(0, N), K{1});
    //Compute number of nodes at each distance and store in sorted_xys at indices N to  N + d_max + 1
    auto num_nodes_at_distance = sorted_xys.subspan(N, d_max + 1).template as_span<K>().subspan(0, d_max + 1);
    primitives::reduce_by_segment(Q, reduce_out.subspan(0,N), sorted_xys.template as_span<K>(), output_keys, num_nodes_at_distance);
    //Compute the centroid of the nodes at each distance
    auto centroids = reduce_out.template as_span<std::array<T,2>>().subspan(0, d_max + 1);
    primitives::transform(Q, summed_coordinates, num_nodes_at_distance, centroids, [](std::array<T,2> a, K b){return a/T(b);});
    //Shift the coordinates of the nodes at each distance by the centroid
    primitives::transform(Q, xys.subspan(0,N), distances, sorted_xys, [centroids](std::array<T,2> a, K b){return a - centroids[b];});
    //Compute the spherical coordinates of the nodes at each distance
    primitives::transform(Q, sorted_xys.subspan(0,N), distances, X, [d_max](std::array<T,2> xy, K dist){
        T dtheta = T(M_PI)/T(d_max+1);
        T phi = dtheta*(dist+0.5);
        T theta = sycl::atan2(xy[0],xy[1]);
        return std::array<T,3>{cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)};
    });
    //Compute center of mass
    std::array<T,3> cm = primitives::reduce(Q, X, std::array<T,3>{0.0, 0.0, 0.0}, [](std::array<T,3> a, std::array<T,3> b){return a + b;});
    cm /= T(N);
    //Shift the coordinates by the center of mass
    primitives::transform(Q, X, X, [cm](std::array<T,3> a){return a - cm;});
    //Compute the average distance between nodes
    T Ravg = primitives::transform_reduce(Q, X, cubic_neighbours.template as_span<std::array<K,3>>(), T{0.0}, Plus{}, [X](std::array<T,3> a, std::array<K,3> neighbours){
        T local_Ravg = 0.0;
        for (size_t i = 0; i < 3; i++){ local_Ravg += norm(a - X[neighbours[i]]); }
        return local_Ravg;
    }) / T(3*N);
    //Scale the coordinates
    T scalerad = 4.0;
    primitives::transform(Q, X, X, [scalerad, Ravg](std::array<T,3> a){return a*scalerad*T(1.5)/Ravg;});
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
//template struct SphericalProjectionFunctor<double,uint16_t>;
//template struct SphericalProjectionFunctor<double,uint32_t>;

//template void spherical_projection<float, uint16_t>(sycl::queue& Q, IsomerBatch<float,uint16_t>& batch, const LaunchPolicy policy);
//template void spherical_projection<double, uint16_t>(sycl::queue& Q, IsomerBatch<double,uint16_t>& batch, const LaunchPolicy policy);