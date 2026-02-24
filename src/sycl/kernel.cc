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

        // Backend-agnostic defaults that are available across SYCL implementations.
        preferred_wg_multiple = 1;
        max_wg_size = dev.get_info<sycl::info::device::max_work_group_size>();
        compile_wg_size = max_wg_size;
        num_regs = 0;

        auto sub_groups = dev.get_info<sycl::info::device::sub_group_sizes>();
        if(!sub_groups.empty()) preferred_wg_multiple = sub_groups[0];
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
        auto sub_groups = dev.get_info<sycl::info::device::sub_group_sizes>();
        auto sub_group_size = sub_groups.empty() ? size_t(1) : sub_groups[0];
        auto ideal_nd_range = isomer_nd_range(N);
        auto round_to_nearest_multiple = [](int_t x, int_t y){return y*((x + y - 1) / y);};
        auto local_size = ideal_nd_range.get_local_range().size();
        auto work_group_size = std::max<size_t>(1, round_to_nearest_multiple(local_size, sub_group_size));
        auto n_work_groups = std::max<size_t>(1, ideal_nd_range.get_global_range().size() / std::max<size_t>(1, local_size));
        auto n_work_groups_per_cu = std::max<size_t>(1, max_wg_size / work_group_size);
        auto n_cus_required = std::max<size_t>(1, n_work_groups / n_work_groups_per_cu);
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