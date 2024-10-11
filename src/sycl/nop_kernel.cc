#include <CL/sycl.hpp>
#include "numeric"
#include <vector>
#include <tuple>
#include <iterator>
#include <type_traits>
#include "forcefield-includes.cc"
#include <fullerenes/sycl-headers/misc-enums.hh>
#include <fullerenes/sycl-headers/sycl-fullerene-structs.hh>

template <typename T, typename K>
void nop_kernel(sycl::queue&Q, FullereneBatch<T,K>& batch, const LaunchPolicy policy){
    TEMPLATE_TYPEDEFS(T,K);
    
    auto batch_view = FullereneBatchView<T,K>(batch);

    if(policy == LaunchPolicy::SYNC) Q.wait();
    Q.submit([&](handler &h) {
        auto N = batch.N_;
        auto Nf = batch.Nf_;
        auto capacity = batch.capacity();
    

        auto num_compute_units = Q.get_device().get_info<info::device::max_compute_units>();
        auto n_blocks_strided = num_compute_units;
        //Create device accessors

        h.parallel_for<class dualize>(sycl::nd_range(sycl::range{N*capacity}, sycl::range{N}), [=](nd_item<1> nditem) {
            auto i = nditem.get_global_id(0);
            auto bid = nditem.get_group_linear_id();
            auto cn = batch_view.d_.A_cubic_[i];
            auto fd = batch_view.d_.deg_[i];
            auto dn = batch_view.d_.A_dual_[i];
            auto xy = batch_view.d_.X_cubic_[i];
            auto X = batch_view.d_.X_dual_[i];
            auto s = batch_view.m_.flags_[i];
            auto it = batch_view.m_.iterations_[i];
            //Do nothing.
            
        });
    });
    if(policy == LaunchPolicy::SYNC) Q.wait();
}

template void nop_kernel<float,uint16_t>(sycl::queue&, FullereneBatch<float,uint16_t>&, const LaunchPolicy);