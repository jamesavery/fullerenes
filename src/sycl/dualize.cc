#include <CL/sycl.hpp>
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>
#include <oneapi/dpl/iterator>
#include <oneapi/dpl/numeric>
#include "numeric"
#include <vector>
#include <iostream>
#include <fstream>
#include <tuple>
#include <iterator>
#include <type_traits>
#include <fullerenes/sycl-isomer-batch.hh>
#include "forcefield-includes.cc"

//Template specialisation for dualize


#define UINT_TYPE_MAX std::numeric_limits<UINT_TYPE>::max()

template<int MaxDegree, typename K>
struct DeviceDualGraph{
    //Check that K is integral
    INT_TYPEDEFS(K);

    const K* dual_neighbours;                          //(Nf x MaxDegree)
    const K* face_degrees;                            //(Nf x 1)
    
    DeviceDualGraph(const K* dual_neighbours, const K* face_degrees) : dual_neighbours(dual_neighbours), face_degrees(face_degrees) {}

    K arc_ix(const K u, const K v) const{
        for (uint8_t j = 0; j < face_degrees[u]; j++){
            if (dual_neighbours[u*MaxDegree + j] == v) return j;
        }

        assert(false);
	return -1;		// Make compiler happy
    }
    K arc_ix(const node2& e){ return arc_ix(e[0], e[1]); }

    /**
     * @brief returns the next node in the clockwise order around u
     * @param v the current node around u
     * @param u the node around which the search is performed
     * @return the next node in the clockwise order around u
     */
    K next(const K u, const K v) const{
        K j = arc_ix(u,v);
        return dual_neighbours[u*MaxDegree + ((j+1)%face_degrees[u])];
    }
    
    /**
     * @brief returns the prev node in the clockwise order around u
     * @param v the current node around u
     * @param u the node around which the search is performed
     * @return the previous node in the clockwise order around u
     */
    K prev(const K u, const K v) const{
        K j = arc_ix(u,v);
        return dual_neighbours[u*MaxDegree + ((j-1+face_degrees[u])%face_degrees[u])];
    }

    /**
     * @brief Find the node that comes next on the face. given by the edge (u,v)
     * @param u Source of the edge.
     * @param v Destination node.
     * @return The node that comes next on the face.
     */
    K next_on_face(const K u, const K v) const{
        return prev(v,u);
    }

    /**
     * @brief Find the node that comes next on the face. given by the edge (u,v)
     * @param u Source of the edge.
     * @param v Destination node.
     * @return The node that comes next on the face.
     */
    K prev_on_face(const K u, const K v) const{
        return next(v,u);
    }

    /**
     * @brief Finds the canonical triangle arc of the triangle (u,v,w)
     * 
     * @param u source node
     * @param v target node
     * @return canonical triangle arc 
     */
    node2 get_canonical_triangle_arc(const K u, const K v) const{
        //In a triangle u, v, w there are only 3 possible representative arcs, the canonical arc is chosen as the one with the smalles source node.
        node2 min_edge = {u,v};
        K w = next(u,v);
        if (v < u && v < w) min_edge = {v, w};
        if (w < u && w < v) min_edge = {w, u};
        return min_edge;
    }

    node2 get_canonical_arc(K u, K v) const{
        int i = 0;
        auto start_node = u;
        node2 min_edge = {u,v};
        while (v!= start_node){
            K w = next_on_face(u, v);
            u = v; v = w;
            if(u < min_edge[0]) min_edge = {u,v};
            ++i;
        }
        //assert(next_on_face(u,v) == start_node);
        return min_edge;
    }
};
template <typename T, typename K>
void dualize(sycl::queue&Q, IsomerBatch<T,K>& batch, const LaunchPolicy policy){
    INT_TYPEDEFS(K);
    constexpr int     MaxDegree = 6;
    constexpr node_t EMPTY_NODE = std::numeric_limits<node_t>::max();
    
    if(policy == LaunchPolicy::SYNC) Q.wait();
    Q.submit([&](handler &h) {
        auto N = batch.N();
        auto Nf = batch.Nf();
        auto capacity = batch.capacity();
	//        std::cout << "Entered dualize" << std::endl;

	// auto d = Q.get_device();
	// std::cout << "Running on " << d.get_info<info::device::name>() << "\n";
	
        //Create local accessors to shared memory
        local_accessor<node_t, 1>   triangle_numbers(Nf*MaxDegree, h);
        local_accessor<node_t, 1>   cached_neighbours(Nf*MaxDegree, h);
        local_accessor<node_t, 1>   cached_degrees(Nf, h);
        local_accessor<node2, 1>    arc_list(N, h);

        //Create device accessors
        accessor     cubic_neighbours_dev(batch.cubic_neighbours, h, write_only);
        accessor     face_degrees_dev    (batch.face_degrees, h, read_only);
        accessor     dual_neighbours_dev (batch.dual_neighbours, h, read_only);
        accessor     statuses_dev        (batch.statuses, h, read_only);
        /* 
        std::cout << N * capacity << std::endl; */
        h.parallel_for(sycl::nd_range(sycl::range{N*capacity}, sycl::range{N}), [=](nd_item<1> nditem) {
            auto cta = nditem.get_group();
            node_t f = nditem.get_local_linear_id();    // Face-node index
            auto isomer = nditem.get_group_linear_id(); // Isomer    index
            if(statuses_dev[isomer] != IsomerStatus::NOT_CONVERGED) return;
	    
            cta.async_work_group_copy(cached_neighbours.get_pointer(), dual_neighbours_dev.get_pointer() + isomer*Nf*MaxDegree, Nf*MaxDegree);
            cta.async_work_group_copy(cached_degrees.get_pointer(),    face_degrees_dev.get_pointer()    + isomer*Nf, Nf);

            DeviceDualGraph<MaxDegree, node_t> FD(cached_neighbours.get_pointer(), cached_degrees.get_pointer());
            node_t canon_arcs[MaxDegree]; for(size_t i=0;i<MaxDegree;i++) canon_arcs[i] = EMPTY_NODE; // JA: memset was bug: it writes byte-values, but node_t is 16/32bit.

            node_t rep_count  = 0;
            sycl::group_barrier(cta);

            if (f < Nf){
                for (node_t i = 0; i < FD.face_degrees[f]; i++){
                    auto canon_arc = FD.get_canonical_triangle_arc(f, FD.dual_neighbours[f*MaxDegree + i]);
                    if (canon_arc[0] == f){
                        canon_arcs[i] = canon_arc[1];
                        rep_count++;
                    }
                }
            }
            sycl::group_barrier(cta);

            node_t scan_result = exclusive_scan_over_group(cta, rep_count, sycl::plus<node_t>{});

            if (f < Nf){
                node_t arc_count = 0;
                for (node_t i = 0; i < FD.face_degrees[f]; i++){
                    if(canon_arcs[i] != std::numeric_limits<node_t>::max()){
                        triangle_numbers[f*MaxDegree + i] = scan_result + arc_count;
                        ++arc_count;
                    }    
                }
            }
            sycl::group_barrier(cta);

            if (f < Nf){
                for (node_t i = 0; i < FD.face_degrees[f]; i++){
                    if(canon_arcs[i] != std::numeric_limits<node_t>::max()){
                        node_t u = triangle_numbers[f*MaxDegree + i];
                        arc_list[u] = {f, canon_arcs[i]};
                    }
                }
            }
            sycl::group_barrier(cta);
//
            auto [u, v] = arc_list[f];
            auto w = FD.next(u,v);
//
            node2 edge_b = FD.get_canonical_triangle_arc(v, u); cubic_neighbours_dev[isomer*N*3 + f*3 + 0] = triangle_numbers[edge_b[0]*MaxDegree + FD.arc_ix(edge_b)];
            node2 edge_c = FD.get_canonical_triangle_arc(w, v); cubic_neighbours_dev[isomer*N*3 + f*3 + 1] = triangle_numbers[edge_c[0]*MaxDegree + FD.arc_ix(edge_c)];
            node2 edge_d = FD.get_canonical_triangle_arc(u, w); cubic_neighbours_dev[isomer*N*3 + f*3 + 2] = triangle_numbers[edge_d[0]*MaxDegree + FD.arc_ix(edge_d)];
        });
    });

    if(policy == LaunchPolicy::SYNC) Q.wait();
}

int roundUp(int numToRound, int multiple) 
{
    assert(multiple);
    return ((numToRound + multiple - 1) / multiple) * multiple;
}

template <typename T, typename K>
void dualize_V1(sycl::queue&Q, IsomerBatch<T,K>& batch, const LaunchPolicy policy){
    INT_TYPEDEFS(K);
    constexpr int     MaxDegree = 6;
    constexpr node_t EMPTY_NODE = std::numeric_limits<node_t>::max();
    
    if(policy == LaunchPolicy::SYNC) Q.wait();
    auto subgroup_size = Q.get_device().get_info<info::device::sub_group_sizes>()[0];
    size_t lcm = roundUp(batch.Nf(), subgroup_size);
    Q.submit([&](handler &h) {
        auto N = batch.N();
        auto Nf = batch.Nf();
        auto capacity = batch.capacity();
	//        std::cout << "Entered dualize" << std::endl;

	// auto d = Q.get_device();
	// std::cout << "Running on " << d.get_info<info::device::name>() << "\n";
	
        //Create local accessors to shared memory
        local_accessor<node_t, 1>   triangle_numbers(Nf*MaxDegree, h);
        local_accessor<node_t, 1>   cached_neighbours(Nf*MaxDegree, h);
        local_accessor<node_t, 1>   cached_degrees(Nf, h);
        local_accessor<node2, 1>    arc_list(N, h);

        //Create device accessors
        accessor     cubic_neighbours_dev(batch.cubic_neighbours, h, write_only);
        accessor     face_degrees_dev    (batch.face_degrees, h, read_only);
        accessor     dual_neighbours_dev (batch.dual_neighbours, h, read_only);
        accessor     statuses_dev        (batch.statuses, h, read_only);
        /* 
        std::cout << N * capacity << std::endl; */
        h.parallel_for(sycl::nd_range(sycl::range{lcm*capacity}, sycl::range{lcm}), [=](nd_item<1> nditem) {
            auto cta = nditem.get_group();
            node_t f = nditem.get_local_linear_id();    // Face-node index
            auto isomer = nditem.get_group_linear_id(); // Isomer    index
            if(statuses_dev[isomer] != IsomerStatus::NOT_CONVERGED) return;
	    
            cta.async_work_group_copy(cached_neighbours.get_pointer(), dual_neighbours_dev.get_pointer() + isomer*Nf*MaxDegree, Nf*MaxDegree);
            cta.async_work_group_copy(cached_degrees.get_pointer(),    face_degrees_dev.get_pointer()    + isomer*Nf, Nf);

            DeviceDualGraph<MaxDegree, node_t> FD(cached_neighbours.get_pointer(), cached_degrees.get_pointer());
            node_t canon_arcs[MaxDegree]; for(size_t i=0;i<MaxDegree;i++) canon_arcs[i] = EMPTY_NODE; // JA: memset was bug: it writes byte-values, but node_t is 16/32bit.

            node_t rep_count  = 0;
            sycl::group_barrier(cta);

            if (f < Nf){
                for (node_t i = 0; i < FD.face_degrees[f]; i++){
                    auto canon_arc = FD.get_canonical_triangle_arc(f, FD.dual_neighbours[f*MaxDegree + i]);
                    if (canon_arc[0] == f){
                        canon_arcs[i] = canon_arc[1];
                        rep_count++;
                    }
                }
            }
            sycl::group_barrier(cta);

            node_t scan_result = exclusive_scan_over_group(cta, rep_count, sycl::plus<node_t>{});

            if (f < Nf){
                node_t arc_count = 0;
                for (node_t i = 0; i < FD.face_degrees[f]; i++){
                    if(canon_arcs[i] != std::numeric_limits<node_t>::max()){
                        triangle_numbers[f*MaxDegree + i] = scan_result + arc_count;
                        ++arc_count;
                    }    
                }
            }
            sycl::group_barrier(cta);

            if (f < Nf){
                for (node_t i = 0; i < FD.face_degrees[f]; i++){
                    if(canon_arcs[i] != std::numeric_limits<node_t>::max()){
                        node_t u = triangle_numbers[f*MaxDegree + i];
                        arc_list[u] = {f, canon_arcs[i]};
                    }
                }
            }
            sycl::group_barrier(cta);
//          
            for(auto tix = f; tix < N; tix += lcm){
                auto [u, v] = arc_list[tix];
                auto w = FD.next(u,v);
//
                node2 edge_b = FD.get_canonical_triangle_arc(v, u); cubic_neighbours_dev[isomer*N*3 + tix*3 + 0] = triangle_numbers[edge_b[0]*MaxDegree + FD.arc_ix(edge_b)];
                node2 edge_c = FD.get_canonical_triangle_arc(w, v); cubic_neighbours_dev[isomer*N*3 + tix*3 + 1] = triangle_numbers[edge_c[0]*MaxDegree + FD.arc_ix(edge_c)];
                node2 edge_d = FD.get_canonical_triangle_arc(u, w); cubic_neighbours_dev[isomer*N*3 + tix*3 + 2] = triangle_numbers[edge_d[0]*MaxDegree + FD.arc_ix(edge_d)];
            }
        });
    });

    if(policy == LaunchPolicy::SYNC) Q.wait();
}

template <typename K>
struct DualWorkingBuffers{
    std::vector<sycl::device> devices;

    std::vector<sycl::buffer<K,1>> cannon_ixs;
    std::vector<sycl::buffer<K,1>> rep_count;
    std::vector<sycl::buffer<K,1>> scan_array;
    std::vector<sycl::buffer<K,1>> triangle_numbers;
    std::vector<sycl::buffer<K,1>> arc_list;

    std::vector<std::array<size_t, 3>> dims; //

    //SYCL lacks a find for retrieving a unique identifier for a given queue, platform, or device so we have to do this manually
    int get_device_index(const sycl::device& device){
        for(int i = 0; i < devices.size(); i++){
            if(devices[i] == device){
                return i;
            }
        }
        return -1;
    }

    DualWorkingBuffers(size_t Nin, size_t Nout, size_t MaxDegree, const sycl::device& device){
        devices = {device};
        cannon_ixs       = std::vector<sycl::buffer<K,1>>(1, sycl::buffer<K,1>(Nin*MaxDegree));
        rep_count        = std::vector<sycl::buffer<K,1>>(1, sycl::buffer<K,1>(Nin));
        scan_array       = std::vector<sycl::buffer<K,1>>(1, sycl::buffer<K,1>(Nin));
        triangle_numbers = std::vector<sycl::buffer<K,1>>(1, sycl::buffer<K,1>(Nin*MaxDegree));
        arc_list         = std::vector<sycl::buffer<K,1>>(1, sycl::buffer<K,1>(Nout*2));
        dims = std::vector<std::array<size_t, 3>>(1, {Nin, Nout, MaxDegree});
    }

    void resize(size_t Nin, size_t Nout, size_t MaxDegree, size_t idx){
        cannon_ixs[idx]         = sycl::buffer<K,1>(Nin*MaxDegree);
        rep_count[idx]          = sycl::buffer<K,1>(Nin);
        scan_array[idx]         = sycl::buffer<K,1>(Nin);
        triangle_numbers[idx]   = sycl::buffer<K,1>(Nin*MaxDegree);
        arc_list[idx]           = sycl::buffer<K,1>(Nout*2);
        dims[idx]               = {Nin, Nout, MaxDegree};
    }
    
    void append_buffers(size_t Nin, size_t Nout, size_t MaxDegree, const sycl::device& device){
        cannon_ixs.push_back        (sycl::buffer<K,1>(Nin*MaxDegree));
        rep_count.push_back         (sycl::buffer<K,1>(Nin));
        scan_array.push_back        (sycl::buffer<K,1>(Nin));
        triangle_numbers.push_back  (sycl::buffer<K,1>(Nin*MaxDegree));
        arc_list.push_back          (sycl::buffer<K,1>(Nout*2));
        devices.push_back           (device);
    }

    void update(size_t Nin, size_t Nout, size_t MaxDegree, const sycl::device& device){
        auto idx = get_device_index(device);
        if(idx == -1){
            append_buffers(Nin, Nout, MaxDegree, device);
            //std::cout << "Appended buffers for device " << device.get_info<sycl::info::device::name>() << " ID: " << devices.size() - 1 << "\n";
        }
        else if (Nin != dims[idx][0] || Nout != dims[idx][1] || MaxDegree != dims[idx][2]){
            resize(Nin, Nout, MaxDegree, idx);
            //std::cout << "Resized buffers for device " << device.get_info<sycl::info::device::name>() << " ID: " << idx << "\n";
        }
    }
};
template <typename K, int MaxDegIn, int MaxDegOut> class DualizeGeneralStep1 {};
template <typename K, int MaxDegIn, int MaxDegOut> class DualizeGeneralStep2 {};
template <typename K, int MaxDegIn, int MaxDegOut> class DualizeGeneralStep3 {};
template <typename K, int MaxDegIn, int MaxDegOut> class DualizeGeneralStep4 {};
template <typename K, int MaxDegIn, int MaxDegOut> class VisDualizeGeneralStep1 {};
template <typename K, int MaxDegIn, int MaxDegOut> class VisDualizeGeneralStep2 {};
template <typename K, int MaxDegIn, int MaxDegOut> class VisDualizeGeneralStep3 {};
template <typename K, int MaxDegIn, int MaxDegOut> class VisDualizeGeneralStep4 {};

template<int MaxDegIn, int MaxDegOut, typename K>
void dualize_general_for_visualization(sycl::queue& Q, sycl::buffer<K,1>& G_in, sycl::buffer<K,1>& Deg_in, sycl::buffer<K,1>& G_out, sycl::buffer<K,1>& Deg_out, int Nin, int Nout, LaunchPolicy policy, bool output_intermediate){
    INT_TYPEDEFS(K);
    if(policy == LaunchPolicy::SYNC) Q.wait();
    sycl::device d = Q.get_device();
    auto exec_pol = oneapi::dpl::execution::make_device_policy(d);

    auto output_buf_to_file = [&](sycl::buffer<K,1>& buf, std::string filename){
        std::vector<K> data(buf.get_count());
        Q.wait();
        std::ofstream file(filename);
        host_accessor acc(buf, read_only);
        for(int i = 0; i < buf.get_count(); i++){
            data[i] = acc[i];
        }
        file.write((char*)data.data(), data.size()*sizeof(K));
    };

    if(output_intermediate) output_buf_to_file(G_in, "G_in_N=" + std::to_string(Nin) + "_dims_" + std::to_string(Nin) + "_X_" + std::to_string(MaxDegIn) + "_.uint16");
    //Find maximum workgroup size
    auto max_workgroup_size = d.get_info<sycl::info::device::max_work_group_size>();
    static DualWorkingBuffers<K> buffers(Nin, Nout, MaxDegIn, d);
    buffers.update(Nin, Nout, MaxDegIn, d);
    sycl::buffer<K,1> face_list(Nout*MaxDegOut); //List of original vertices comprising each face
    sycl::buffer<K,1> neighbouring_output_faces(MaxDegIn*Nin); //List of faces that each vertex is a part of

    auto idx = buffers.get_device_index(d);
    size_t workgroup_size1 = std::min((int)max_workgroup_size, Nin);
    size_t workgroup_size2 = std::min((int)max_workgroup_size, Nout);
    size_t grid_size1 = roundUp(Nin, workgroup_size1);
    size_t grid_size2 = roundUp(Nout, workgroup_size2);
    Q.submit([&](handler &h) {
        accessor G_in_acc   (G_in,   h, read_only);
        accessor Deg_in_acc (Deg_in, h, read_only);
        accessor G_out_acc  (G_out,  h, write_only);
        accessor Deg_out_acc(Deg_out, h, write_only);
        accessor cannon_ixs_acc         (buffers.cannon_ixs[idx], h, read_write);
        accessor rep_count_acc          (buffers.rep_count[idx], h, read_write);
        accessor scan_array_acc         (buffers.scan_array[idx], h, read_write);
        accessor triangle_numbers_acc   (buffers.triangle_numbers[idx], h, read_write);
        accessor arc_list_acc           (buffers.arc_list[idx], h, read_write);

        
        h.parallel_for<VisDualizeGeneralStep1<K,MaxDegIn, MaxDegOut>>(nd_range(range{grid_size1}, range{workgroup_size1}), [=](nd_item<1> nditem) {
            auto thid = nditem.get_global_linear_id();
            DeviceDualGraph<MaxDegIn, K> FD(G_in_acc.get_pointer(), Deg_in_acc.get_pointer());
            K rep_count = 0;
            if (thid < Nin){
                for (int i = 0; i < FD.face_degrees[thid]; i++){
                    auto canon_arc = FD.get_canonical_arc(thid, FD.dual_neighbours[thid*MaxDegIn + i]);
                    if (canon_arc[0] == thid){
                        cannon_ixs_acc[thid*MaxDegIn + rep_count] = i;
                        ++rep_count;
                    }
                }
                rep_count_acc[thid] = rep_count;
                scan_array_acc[thid] = rep_count;
            }
        });
    });
    if (output_intermediate){
        output_buf_to_file(buffers.cannon_ixs[idx], "cannon_ixs_N=" + std::to_string(Nout) + "_dims_" + std::to_string(Nin) + "_X_" + std::to_string(MaxDegIn) + "_.uint16");
        output_buf_to_file(buffers.rep_count[idx], "rep_count_N=" + std::to_string(Nout) + "_dim_" + std::to_string(Nin) + "_.uint16");
    }


    oneapi::dpl::exclusive_scan(
            exec_pol,
            oneapi::dpl::begin(buffers.scan_array[idx]),
            oneapi::dpl::end(buffers.scan_array[idx]),
            oneapi::dpl::begin(buffers.scan_array[idx]),
            0);

    if(output_intermediate) output_buf_to_file(buffers.scan_array[idx], "scan_array_N=" + std::to_string(Nout) + "_dim_" + std::to_string(Nin) + "_.uint16");

    Q.submit([&](handler &h) {
        accessor G_in_acc   (G_in,   h, read_only);
        accessor Deg_in_acc (Deg_in, h, read_only);
        accessor G_out_acc  (G_out,  h, write_only);
        accessor Deg_out_acc(Deg_out, h, write_only);
        accessor cannon_ixs_acc         (buffers.cannon_ixs[idx], h, read_only);
        accessor rep_count_acc          (buffers.rep_count[idx], h, read_only);
        accessor scan_array_acc         (buffers.scan_array[idx], h, read_only);
        accessor triangle_numbers_acc   (buffers.triangle_numbers[idx], h, read_write);
        accessor arc_list_acc           (buffers.arc_list[idx], h, read_write);

        h.parallel_for<VisDualizeGeneralStep2<K,MaxDegIn, MaxDegOut>>(nd_range(range{grid_size1}, range{workgroup_size1}), [=](nd_item<1> nditem) {
            auto idx = nditem.get_global_linear_id();
            DeviceDualGraph<MaxDegIn, K> FD(G_in_acc.get_pointer(), Deg_in_acc.get_pointer());
            if (idx < Nin){
            K rep_count = rep_count_acc[idx];
            K scan_result = scan_array_acc[idx];
            for (int ii = 0; ii < rep_count; ii++){
                K i = cannon_ixs_acc[idx*MaxDegIn + ii];
                K triangle_id = scan_result + ii;
                triangle_numbers_acc[idx*MaxDegIn + i] = triangle_id;
                arc_list_acc[triangle_id*2 + 0] = idx;
                arc_list_acc[triangle_id*2 + 1] = FD.dual_neighbours[idx*MaxDegIn + i];
            }
            }
        });
    });
    if(output_intermediate){
        output_buf_to_file(buffers.triangle_numbers[idx], "triangle_numbers_N=" + std::to_string(Nout) + "_dims_" + std::to_string(Nin) + "_X_" + std::to_string(MaxDegIn) + "_.uint16");
        output_buf_to_file(buffers.arc_list[idx], "arc_list_N=" + std::to_string(Nout) + "_dims_" + std::to_string(Nin) + "_X_" + std::to_string(2) + "_.uint16");
    }

    
    Q.submit([&](handler &h) {
        accessor G_in_acc   (G_in,   h, read_only);
        accessor Deg_in_acc (Deg_in, h, read_only);
        accessor G_out_acc  (G_out,  h, write_only);
        accessor Deg_out_acc(Deg_out, h, write_only);
        accessor triangle_numbers_acc   (buffers.triangle_numbers[idx], h, read_write);
        accessor arc_list_acc           (buffers.arc_list[idx], h, read_write);
        accessor face_list_acc (face_list, h, write_only);

        h.parallel_for<VisDualizeGeneralStep3<K,MaxDegIn, MaxDegOut>>(nd_range(range{grid_size2}, range{workgroup_size2}), [=](nd_item<1> nditem) {
            auto tidx = nditem.get_global_linear_id();
            DeviceDualGraph<MaxDegIn, K> FD(G_in_acc.get_pointer(), Deg_in_acc.get_pointer());
            if (tidx < Nout){
            K u = arc_list_acc[tidx*2 + 0]; K v = arc_list_acc[tidx*2 + 1];
            auto n_count = 0;
            auto u0 = u;
            node2 edge = FD.get_canonical_arc(v, u); G_out_acc[tidx*MaxDegOut] = triangle_numbers_acc[edge[0]*MaxDegIn + FD.arc_ix(edge)];
            face_list_acc[tidx*MaxDegOut] = u;
            face_list_acc[tidx*MaxDegOut + 1] = v;
            while(v != u0 && n_count < MaxDegOut){
                ++n_count;
                auto w = FD.next_on_face(u,v);
                if(w != u0) face_list_acc[tidx*MaxDegOut + n_count + 1] = w;
                u = v; v = w;
                edge = FD.get_canonical_arc(v, u); G_out_acc[tidx*MaxDegOut + n_count] = triangle_numbers_acc[edge[0]*MaxDegIn + FD.arc_ix(edge)];
            }
            Deg_out_acc[tidx] = n_count+1;
            }
        });
    });
    if(output_intermediate){ 
    Q.submit([&](handler &h) {
        accessor G_in_acc   (G_in,   h, read_only);
        accessor Deg_in_acc (Deg_in, h, read_only);
        accessor G_out_acc  (G_out,  h, write_only);
        accessor Deg_out_acc(Deg_out, h, write_only);
        accessor triangle_numbers_acc           (buffers.triangle_numbers[idx], h, read_write);
        accessor arc_list_acc                   (buffers.arc_list[idx], h, read_write);
        accessor neighbouring_output_faces_acc  (neighbouring_output_faces, h, write_only);

        h.parallel_for<VisDualizeGeneralStep4<K,MaxDegIn, MaxDegOut>>(nd_range(range{grid_size1}, range{workgroup_size1}), [=](nd_item<1> nditem) {
            auto tidx = nditem.get_global_linear_id();
            DeviceDualGraph<MaxDegIn, K> FD(G_in_acc.get_pointer(), Deg_in_acc.get_pointer());
            if (tidx < Nin){
                K u = tidx;
                for (int i = 0; i < FD.face_degrees[u]; i++){
                    auto v = FD.dual_neighbours[u*MaxDegIn + i];
                    auto arc = FD.get_canonical_arc(u,v);
                    auto triangle_id = triangle_numbers_acc[arc[0]*MaxDegIn + FD.arc_ix(arc)];
                    neighbouring_output_faces_acc[u*MaxDegIn + i] = triangle_id;
                }
            }
        });
    });
    }

    if(output_intermediate) {
        output_buf_to_file(G_out, "G_out_N=" + std::to_string(Nout) + "_dims_" + std::to_string(Nin) + "_X_" + std::to_string(MaxDegOut) + "_.uint16");
        output_buf_to_file(Deg_out, "Deg_out_N=" + std::to_string(Nout) + "_dim_" + std::to_string(Nout) + "_.uint16");
        output_buf_to_file(face_list, "face_list_N=" + std::to_string(Nout) + "_dims_" + std::to_string(Nout) + "_X_" + std::to_string(MaxDegOut) + "_.uint16");
        output_buf_to_file(neighbouring_output_faces, "neighbouring_output_faces_N=" + std::to_string(Nin) + "_dims_" + std::to_string(Nin) + "_X_" + std::to_string(MaxDegIn) + "_.uint16");
    }
    
    if(policy == LaunchPolicy::SYNC) Q.wait();
}

template<int MaxDegIn, int MaxDegOut, typename K>
void dualize_general(sycl::queue& Q, sycl::buffer<K,1>& G_in, sycl::buffer<K,1>& Deg_in, sycl::buffer<K,1>& G_out, sycl::buffer<K,1>& Deg_out, int Nin, int Nout, LaunchPolicy policy){
    INT_TYPEDEFS(K);
    if(policy == LaunchPolicy::SYNC) Q.wait();
    sycl::device d = Q.get_device();
    auto exec_pol = oneapi::dpl::execution::make_device_policy(d);

    //Find maximum workgroup size
    auto max_workgroup_size = d.get_info<sycl::info::device::max_work_group_size>();
    static DualWorkingBuffers<K> buffers(Nin, Nout, MaxDegIn, d);
    buffers.update(Nin, Nout, MaxDegIn, d);

    auto idx = buffers.get_device_index(d);
    size_t workgroup_size1 = std::min((int)max_workgroup_size, Nin);
    size_t workgroup_size2 = std::min((int)max_workgroup_size, Nout);
    size_t grid_size1 = roundUp(Nin, workgroup_size1);
    size_t grid_size2 = roundUp(Nout, workgroup_size2);
    Q.submit([&](handler &h) {
        accessor G_in_acc   (G_in,   h, read_only);
        accessor Deg_in_acc (Deg_in, h, read_only);
        accessor G_out_acc  (G_out,  h, write_only);
        accessor Deg_out_acc(Deg_out, h, write_only);
        accessor cannon_ixs_acc         (buffers.cannon_ixs[idx], h, read_write);
        accessor rep_count_acc          (buffers.rep_count[idx], h, read_write);
        accessor scan_array_acc         (buffers.scan_array[idx], h, read_write);
        accessor triangle_numbers_acc   (buffers.triangle_numbers[idx], h, read_write);
        accessor arc_list_acc           (buffers.arc_list[idx], h, read_write);

        
        h.parallel_for<DualizeGeneralStep1<K,MaxDegIn, MaxDegOut>>(nd_range(range{grid_size1}, range{workgroup_size1}), [=](nd_item<1> nditem) {
            auto thid = nditem.get_global_linear_id();
            DeviceDualGraph<MaxDegIn, K> FD(G_in_acc.get_pointer(), Deg_in_acc.get_pointer());
            K rep_count = 0;
            if (thid < Nin){
                for (int i = 0; i < FD.face_degrees[thid]; i++){
                    auto canon_arc = FD.get_canonical_arc(thid, FD.dual_neighbours[thid*MaxDegIn + i]);
                    if (canon_arc[0] == thid){
                        cannon_ixs_acc[thid*MaxDegIn + rep_count] = i;
                        ++rep_count;
                    }
                }
                rep_count_acc[thid] = rep_count;
                scan_array_acc[thid] = rep_count;
            }
        });
    });
    oneapi::dpl::exclusive_scan(
            exec_pol,
            oneapi::dpl::begin(buffers.scan_array[idx]),
            oneapi::dpl::end(buffers.scan_array[idx]),
            oneapi::dpl::begin(buffers.scan_array[idx]),
            0);

    Q.submit([&](handler &h) {
        accessor G_in_acc   (G_in,   h, read_only);
        accessor Deg_in_acc (Deg_in, h, read_only);
        accessor G_out_acc  (G_out,  h, write_only);
        accessor Deg_out_acc(Deg_out, h, write_only);
        accessor cannon_ixs_acc         (buffers.cannon_ixs[idx], h, read_only);
        accessor rep_count_acc          (buffers.rep_count[idx], h, read_only);
        accessor scan_array_acc         (buffers.scan_array[idx], h, read_only);
        accessor triangle_numbers_acc   (buffers.triangle_numbers[idx], h, read_write);
        accessor arc_list_acc           (buffers.arc_list[idx], h, read_write);

        h.parallel_for<DualizeGeneralStep2<K,MaxDegIn, MaxDegOut>>(nd_range(range{grid_size1}, range{workgroup_size1}), [=](nd_item<1> nditem) {
            auto idx = nditem.get_global_linear_id();
            DeviceDualGraph<MaxDegIn, K> FD(G_in_acc.get_pointer(), Deg_in_acc.get_pointer());
            if (idx < Nin){
            K rep_count = rep_count_acc[idx];
            K scan_result = scan_array_acc[idx];
            for (int ii = 0; ii < rep_count; ii++){
                K i = cannon_ixs_acc[idx*MaxDegIn + ii];
                K triangle_id = scan_result + ii;
                triangle_numbers_acc[idx*MaxDegIn + i] = triangle_id;
                arc_list_acc[triangle_id*2 + 0] = idx;
                arc_list_acc[triangle_id*2 + 1] = FD.dual_neighbours[idx*MaxDegIn + i];
            }
            }
        });
    });
    
    Q.submit([&](handler &h) {
        accessor G_in_acc   (G_in,   h, read_only);
        accessor Deg_in_acc (Deg_in, h, read_only);
        accessor G_out_acc  (G_out,  h, write_only);
        accessor Deg_out_acc(Deg_out, h, write_only);
        accessor triangle_numbers_acc   (buffers.triangle_numbers[idx], h, read_write);
        accessor arc_list_acc           (buffers.arc_list[idx], h, read_write);

        h.parallel_for<DualizeGeneralStep3<K,MaxDegIn, MaxDegOut>>(nd_range(range{grid_size2}, range{workgroup_size2}), [=](nd_item<1> nditem) {
            auto tidx = nditem.get_global_linear_id();
            DeviceDualGraph<MaxDegIn, K> FD(G_in_acc.get_pointer(), Deg_in_acc.get_pointer());
            if (tidx < Nout){
            K u = arc_list_acc[tidx*2 + 0]; K v = arc_list_acc[tidx*2 + 1];
            auto n_count = 0;
            auto u0 = u;
            node2 edge = FD.get_canonical_arc(v, u); G_out_acc[tidx*MaxDegOut] = triangle_numbers_acc[edge[0]*MaxDegIn + FD.arc_ix(edge)];
            while(v != u0 && n_count < MaxDegOut){
                ++n_count;
                auto w = FD.next_on_face(u,v);
                u = v; v = w;
                edge = FD.get_canonical_arc(v, u); G_out_acc[tidx*MaxDegOut + n_count] = triangle_numbers_acc[edge[0]*MaxDegIn + FD.arc_ix(edge)];
            }
            Deg_out_acc[tidx] = n_count+1;
            }
        });
    });
    
    if(policy == LaunchPolicy::SYNC) Q.wait();

}

template void dualize<float,uint16_t>(sycl::queue&Q, IsomerBatch<float,uint16_t>& batch, const LaunchPolicy policy);
template void dualize_V1<float,uint16_t>(sycl::queue&Q, IsomerBatch<float,uint16_t>& batch, const LaunchPolicy policy);
template void dualize<double,uint16_t>(sycl::queue&Q, IsomerBatch<double,uint16_t>& batch, const LaunchPolicy policy);
template void dualize_V1<double,uint16_t>(sycl::queue&Q, IsomerBatch<double,uint16_t>& batch, const LaunchPolicy policy);
template void dualize_general<6, 3, uint16_t>(sycl::queue&Q, sycl::buffer<uint16_t,1>& G_in, sycl::buffer<uint16_t,1>& Deg_in, sycl::buffer<uint16_t,1>& G_out, sycl::buffer<uint16_t,1>& Deg_out, int Nin, int Nout, LaunchPolicy policy);
template void dualize_general<3, 6, uint16_t>(sycl::queue&Q, sycl::buffer<uint16_t,1>& G_in, sycl::buffer<uint16_t,1>& Deg_in, sycl::buffer<uint16_t,1>& G_out, sycl::buffer<uint16_t,1>& Deg_out, int Nin, int Nout, LaunchPolicy policy);
template void dualize_general_for_visualization<6, 3, uint16_t>(sycl::queue&Q, sycl::buffer<uint16_t,1>& G_in, sycl::buffer<uint16_t,1>& Deg_in, sycl::buffer<uint16_t,1>& G_out, sycl::buffer<uint16_t,1>& Deg_out, int Nin, int Nout, LaunchPolicy policy, bool output_intermediate);
template void dualize_general_for_visualization<3, 6, uint16_t>(sycl::queue&Q, sycl::buffer<uint16_t,1>& G_in, sycl::buffer<uint16_t,1>& Deg_in, sycl::buffer<uint16_t,1>& G_out, sycl::buffer<uint16_t,1>& Deg_out, int Nin, int Nout, LaunchPolicy policy, bool output_intermediate);
