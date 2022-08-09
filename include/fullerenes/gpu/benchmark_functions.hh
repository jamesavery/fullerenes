#include "isomerspace_kernel.hh"
#include "chrono"
namespace cuda_benchmark {
    bool test_global_scan(const int n_elements, const int n_times);
    
    std::chrono::microseconds benchmark_scan(const int n_elements, const int n_times, const int scan_version);

    std::chrono::nanoseconds benchmark_reduction(const int n_elements, const int n_times, const int reduction_version);

    size_t n_blocks(const int n_elements);

}
