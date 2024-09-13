#pragma once
#include "queue-impl.cc"

template <typename KernelName>
struct LaunchConfig{
    using int_t = size_t;
    int_t preferred_wg_multiple;
    int_t max_wg_size;
    int_t compile_wg_size;
    int_t num_regs;
    sycl::device dev;
    
    LaunchConfig(SyclQueue& Q){
        dev = Q->get_device();
        auto context = Q->get_context();
        auto kernel_bundle = sycl::get_kernel_bundle<sycl::bundle_state::executable>(context, {dev});
        try{
            auto kernel = kernel_bundle.get_kernel(sycl::get_kernel_id<KernelName>());
            preferred_wg_multiple = kernel.template get_info<sycl::info::kernel_device_specific::preferred_work_group_size_multiple>(dev);
            max_wg_size = kernel.template get_info<sycl::info::kernel_device_specific::work_group_size>(dev);
            compile_wg_size = kernel.template get_info<sycl::info::kernel_device_specific::compile_work_group_size>(dev)[0];
            num_regs = kernel.template get_info<sycl::info::kernel_device_specific::ext_codeplay_num_regs>(dev);
        } catch (sycl::exception& e){
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }

    bool is_batch_execution_possible(int_t N){
        return N <= max_wg_size;
    }

    sycl::nd_range<1> batch_nd_range(int_t N, int_t N_tasks){
        if(!is_batch_execution_possible(N)) assert(!"N must be less than or equal to the maximum work group size");
        return sycl::nd_range<1>(sycl::range<1>(N*N_tasks), sycl::range<1>(N));
    }

    sycl::nd_range<1> isomer_nd_range(int_t N){
        auto n_work_groups  = (N + max_wg_size - 1) / max_wg_size;
        auto best_wg_size   = (N + n_work_groups - 1) / n_work_groups;
        return sycl::nd_range<1>(sycl::range<1>(best_wg_size*n_work_groups), sycl::range<1>(best_wg_size));
    }

    size_t max_concurrent_launches(int_t N){
        if(dev.is_cpu()) return dev.get_info<sycl::info::device::max_compute_units>();
        auto max_comp_units = dev.get_info<sycl::info::device::max_compute_units>();
        auto sub_group_size = dev.get_info<sycl::info::device::sub_group_sizes>()[0];
        auto ideal_nd_range = isomer_nd_range(N);
        auto round_to_nearest_multiple = [](int_t x, int_t y){return y*((x + y - 1) / y);};
        auto work_group_size = round_to_nearest_multiple(ideal_nd_range.get_local_range().size(), sub_group_size);
        auto n_work_groups = ideal_nd_range.get_global_range().size() / ideal_nd_range.get_local_range().size();
        auto n_work_groups_per_cu = max_wg_size / work_group_size;
        auto n_cus_required = n_work_groups / n_work_groups_per_cu;
        auto max_concurrent_launches = max_comp_units / n_cus_required;
        return std::min((int_t)max_concurrent_launches, (int_t)max_comp_units);
    }

    inline friend std::ostream& operator<<(std::ostream& os, const LaunchConfig& config){
        os << "Kernel launch configuration for kernel: " << typeid(KernelName).name() << std::endl;
        os << "Preferred work group size multiple: " << config.preferred_wg_multiple << std::endl;
        os << "Max work group size: " << config.max_wg_size << std::endl;
        os << "Compile work group size: " << config.compile_wg_size << std::endl;
        os << "Number of registers: " << config.num_regs << std::endl;
        return os;
    }
};