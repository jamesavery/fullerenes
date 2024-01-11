#include <CL/sycl.hpp>
#include "numeric"
#include <vector>
#include <tuple>
#include <iterator>
#include <type_traits>
#include <fullerenes/sycl-isomer-batch.hh>
#include "forcefield-includes.cc"

template <typename T, typename K>
void nop_kernel(sycl::queue&Q, IsomerBatch<T,K>& batch, const LaunchPolicy policy){
    TEMPLATE_TYPEDEFS(T,K);
    
    if(policy == LaunchPolicy::SYNC) Q.wait();
    Q.submit([&](handler &h) {
        auto N = batch.N();
        auto Nf = batch.Nf();
        auto capacity = batch.capacity();
    

        auto num_compute_units = Q.get_device().get_info<info::device::max_compute_units>();
        auto n_blocks_strided = num_compute_units;
        //Create device accessors
        accessor     cubic_neighbours_dev(batch.cubic_neighbours, h, read_write);
        accessor     face_degrees_dev(batch.face_degrees, h, read_write);
        accessor     dual_neighbours_dev(batch.dual_neighbours, h, read_write);
        accessor     xys_dev(batch.xys, h, read_write);
        accessor     X_dev(batch.X, h, read_write);
        accessor     statuses_dev(batch.statuses, h, read_write);
        accessor     iterations_dev(batch.iterations, h, read_write);

        h.parallel_for<class dualize>(sycl::nd_range(sycl::range{N*capacity}, sycl::range{N}), [=](nd_item<1> nditem) {
            auto i = nditem.get_global_id(0);
            auto bid = nditem.get_group_linear_id();
            auto cn = cubic_neighbours_dev[i];
            auto fd = face_degrees_dev[i];
            auto dn = dual_neighbours_dev[i];
            auto xy = xys_dev[i];
            auto X = X_dev[i];
            auto s = statuses_dev[i];
            auto it = iterations_dev[i];
            //Do nothing.
            
        });
    });
    if(policy == LaunchPolicy::SYNC) Q.wait();
}

template void nop_kernel<float,uint16_t>(sycl::queue&, IsomerBatch<float,uint16_t>&, const LaunchPolicy);