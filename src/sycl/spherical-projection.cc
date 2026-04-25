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
K multiple_source_shortest_paths(const sycl::group<1>& cta, const std::span<std::array<K,3>> cubic_neighbours,const local_accessor<int,1>& distances, const local_accessor<K,1>& smem){
    INT_TYPEDEFS(K);
    auto N = cta.get_local_linear_range();
    auto tid = cta.get_local_linear_id();
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


// ---------------------------------------------------------------------------
// View-based batch implementation (Phase 7)
// ---------------------------------------------------------------------------
template <typename T, typename K>
static SyclEvent spherical_projection_view_batch_impl(
    SyclQueue& Q,
    batch::BatchView<Spanify::RSRAdjacencyView<K>> graph,
    std::span<std::array<T,2>> layout_2d,
    std::span<std::array<T,3>> xyz_3d,
    batch::BatchStateView state)
{
    TEMPLATE_TYPEDEFS(T,K);
    constexpr real_t scalerad = 4.0;

    auto [adj_flat, deg_flat, twin_flat] = graph.spans();
    (void)deg_flat; (void)twin_flat;

    std::span<std::array<K,3>> A_cubic(
        reinterpret_cast<std::array<K,3>*>(adj_flat.data()),
        adj_flat.size() / 3);
    std::span<coord2d> layout_cd(
        reinterpret_cast<coord2d*>(layout_2d.data()),
        layout_2d.size());

    auto statuses   = state.status;
    const int N        = graph.N();
    const int capacity = graph.size();

    SyclEventImpl projection_done = Q->submit([&](handler& h) {
        local_accessor<node_t, 1>  work_queue_memory(N*2, h);
        local_accessor<int, 1>     smem(N, h);
        local_accessor<coord2d, 1> atomic_coordinate_memory(N, h);
        local_accessor<coord3d, 1> xyz_smem(N, h);

        h.parallel_for(sycl::nd_range(sycl::range(N*capacity), sycl::range(N)),
        [=](nd_item<1> nditem) {
            auto cta        = nditem.get_group();
            auto tid        = nditem.get_local_linear_id();
            auto isomer_idx = nditem.get_group_linear_id();

            if (!statuses[isomer_idx].is_set(StatusEnum::FULLERENEGRAPH_PREPARED)) return;

            const auto cubic_neighbours = A_cubic.subspan(isomer_idx * N, N);
            const auto xys_acc = layout_cd.subspan(isomer_idx * N, N);
            auto xyz_acc = xyz_3d.subspan(isomer_idx * N, N);

            atomic_coordinate_memory[tid] = {T(0), T(0)};
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
            real_t Ravg = real_t(0);
            xyz_smem[tid] = xyz;
            sycl::group_barrier(cta);
            real_t local_Ravg = real_t(0);
            for (size_t i = 0; i < 3; i++) { local_Ravg += norm(xyz_smem[tid] - xyz_smem[neighbours[i]]); }
            Ravg = sycl::reduce_over_group(cta, local_Ravg, sycl::plus<real_t>{})/real_t(3*N);
            xyz *= scalerad*real_t(1.5)/Ravg;
            xyz_acc[tid] = xyz;
            if (tid == 0) statuses[isomer_idx] |= StatusEnum::NOT_CONVERGED;
        });
    });
    return SyclEvent(std::move(projection_done));
}

template <typename T, typename K>
SyclEvent SphericalProjectionFunctor<T,K>::compute(
    SyclQueue& Q,
    batch::BatchView<Spanify::RSRAdjacencyView<K>> graph,
    std::span<std::array<T,2>> layout_2d,
    std::span<std::array<T,3>> xyz_3d,
    batch::BatchStateView state) {
    return spherical_projection_view_batch_impl<T,K>(Q, graph, layout_2d, xyz_3d, state);
}

template struct SphericalProjectionFunctor<float,uint16_t>;
template struct SphericalProjectionFunctor<float,uint32_t>;
template struct SphericalProjectionFunctor<double,uint16_t>;
template struct SphericalProjectionFunctor<double,uint32_t>;