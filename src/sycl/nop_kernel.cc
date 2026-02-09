#include <sycl/sycl.hpp>
#include "numeric"
#include <vector>
#include <tuple>
#include <iterator>
#include <type_traits>
#include <fullerenes/kernel-headers/base-functor.hh>
#include "forcefield-includes.cc"

template <typename T, typename K>
void nop_kernel(sycl::queue&Q, FullereneBatch<T,K>& batch, const LaunchPolicy policy){
    TEMPLATE_TYPEDEFS(T,K);
    
    auto batch_view = FullereneBatchView<T,K>(batch);

    if(policy == LaunchPolicy::SYNC) Q.wait();
    Q.submit([&](handler &h) {
        auto N = batch.N_;
        auto capacity = batch.capacity();

        h.parallel_for<class dualize>(sycl::nd_range(sycl::range{N*capacity}, sycl::range{N}), [=](nd_item<1> nditem) {
            auto i = nditem.get_global_id(0);
            [[maybe_unused]] auto cn = batch_view.d_.A_cubic_[i];
            [[maybe_unused]] auto fd = batch_view.d_.deg_[i];
            [[maybe_unused]] auto dn = batch_view.d_.A_dual_[i];
            [[maybe_unused]] auto xy = batch_view.d_.X_cubic_[i];
            [[maybe_unused]] auto X = batch_view.d_.X_dual_[i];
            [[maybe_unused]] auto s = batch_view.m_.flags_[i];
            [[maybe_unused]] auto it = batch_view.m_.iterations_[i];
            //Do nothing.
            
        });
    });
    if(policy == LaunchPolicy::SYNC) Q.wait();
}

template void nop_kernel<float,uint16_t>(sycl::queue&, FullereneBatch<float,uint16_t>&, const LaunchPolicy);