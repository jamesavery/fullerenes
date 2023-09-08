#include <CL/sycl.hpp>
#include "numeric"
#include <vector>
#include <tuple>
#include <iterator>
#include <type_traits>
#include <fullerenes/sycl-isomer-batch.hh>
namespace spp {
    #include "forcefield-includes.cc"
}
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
CuDeque(local_accessor<T,1> memory, const int capacity): array(memory), front(0), back(0), q_size(0), capacity(capacity) {}
    
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
        return q_size == 0;
    }

    /**
     * @brief  Returns true if the queue is full and false otherwise. 
     * @param  None
     * @retval True if the queue is full and false otherwise.
     */
    bool full(){
        return q_size == capacity;
    }
    
    /**
     * @brief  This function is used to pop the first element of the queue.
     * @param  None
     * @retval First element of the queue
     */
    T pop_front(){
        if (!empty()){
            T return_val = array[front];
            front = (front + 1) % capacity ;
            q_size--;
            return return_val;
        }
        assert(false);
        return T(); //Compiler wants a return statement
    }

    /**
     * @brief Returns the last element of the queue and removes it from the queue.
     * @param None
     * @return The last element of the queue
     */
    T pop_back(){
        if (!empty())
        {
            T return_val = array[back];
            back = back > 0 ? back-1 : capacity-1;
            q_size--;
            return return_val;
        }
        assert(false);
        return T(); //Compiler wants a return statement
    }

    /** @brief Insert a value into the back of the queue
     *  @param val the value to insert
     */
    void push_back(T val){
        assert(!full());
        back = (back + 1) % capacity;
        array[back] = val;
        q_size++;
    }

    /** @brief Insert a value into the front of the queue
     *  @param val the value to insert
     */
    void push_front(T val){
        assert(!full());
        front = front > 0 ? front-1 : capacity-1;
        array[front] = val;
        q_size++;
    }
};

template <typename K>
K multiple_source_shortest_paths(const sycl::group<1>& cta, const accessor<K,1>& cubic_neighbours, local_accessor<K,1>& distances, local_accessor<K,1>& smem){
    INT_TYPEDEFS(K);
    using namespace spp;
    auto N = cta.get_local_linear_range();
    auto tid = cta.get_local_linear_id();
    auto isomer_idx = cta.get_group_linear_id();
    DeviceCubicGraph FG(cubic_neighbours, isomer_idx*3*N);
    node_t outer_face[6]; memset(outer_face, std::numeric_limits<node_t>::max(), 6*sizeof(node_t));
    uint8_t Nface = FG.get_face_oriented(0,FG.cubic_neighbours[0], outer_face);
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
                auto w = FG.cubic_neighbours[v*3 + i];
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
void spherical_projection(sycl::queue& Q, IsomerBatch<T,K>& batch, const LaunchPolicy policy){
    using namespace spp;
    TEMPLATE_TYPEDEFS(T,K);
    constexpr real_t scalerad = 4.0;
    if(policy == LaunchPolicy::SYNC) Q.wait();
    Q.submit([&](handler& h) {
        auto N = batch.N();
        auto Nf = batch.Nf();
        auto max_iter = N * 10;
        accessor<real_t, 1> xys_acc(batch.xys, h, write_only); 
        accessor<node_t, 1> cubic_neighbours(batch.cubic_neighbours, h, read_only);
        accessor<IsomerStatus, 1> statuses(batch.statuses, h, read_only);
        accessor<coord3d, 1> xyz_acc(batch.X, h, write_only);

        local_accessor<node_t, 1>   work_queue_memory(N*2, h);
        local_accessor<node_t, 1>   smem(N, h);
        local_accessor<coord2d, 1>  atomic_coordinate_memory(N, h);
        local_accessor<coord3d, 1>  xyz_smem(N, h);


        h.parallel_for<class tutte>(sycl::nd_range(sycl::range(N*batch.capacity()), sycl::range(N)), [=](nd_item<1> nditem) {
            auto cta = nditem.get_group();
            auto tid = nditem.get_local_linear_id();
            auto isomer_idx = nditem.get_group_linear_id();
            if (statuses[isomer_idx] != IsomerStatus::EMPTY){
            NodeNeighbours node_graph(cubic_neighbours, cta);
            node3 neighbours = node_graph.cubic_neighbours;
            node_t distance = multiple_source_shortest_paths(cta, cubic_neighbours, smem, work_queue_memory);
            node_t d_max = reduce_over_group(cta, distance, maximum<node_t>{});
            smem[tid] = 0;
            sycl::group_barrier(cta);
            sycl::atomic_fetch_add(&smem[distance], 1);
            sycl::group_barrier(cta);
            node_t num_same_dist = smem[distance];
            sycl::group_barrier(cta);
            coord2d xys = xys_acc[isomer_idx*N + tid];
            sycl::atomic_fetch_add(&atomic_coordinate_memory[distance][0], xys[0]);
            sycl::atomic_fetch_add(&atomic_coordinate_memory[distance][1], xys[1]);
            sycl::group_barrier(cta);
            coord2d centroid = atomic_coordinate_memory[distance]/num_same_dist;
            coord2d xy = xys - centroid;
            real_t dtheta = real_t(M_PI)/real_t(d_max+1);
            real_t phi = dtheta*(distance+0.5);
            real_t theta = sycl::atan2(xy[0],xy[1]);
            coord2d spherical_coords = {theta, phi};
            coord3d xyz = {cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)};
            
            coord3d cm = reduce_over_group(cta, xyz, plus<coord3d>{})/real_t(N);
            xyz -= cm;
            real_t Ravg = 0.0;
            xyz_smem[tid] = xyz;
            sycl::group_barrier(cta);
            real_t local_Ravg = 0.0;
            for (size_t i = 0; i < 3; i++){ local_Ravg += norm(xyz_smem[tid] - xyz_smem[neighbours[i]]); }
            Ravg = reduce_over_group(cta, local_Ravg, plus<real_t>{})/real_t(3*N);
            xyz *= scalerad/Ravg;
            xyz_acc[isomer_idx*N + tid] = xyz;}
        });
    });
    if(policy == LaunchPolicy::SYNC) Q.wait();
}