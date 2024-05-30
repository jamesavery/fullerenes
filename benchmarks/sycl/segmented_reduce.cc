#include <CL/sycl.hpp>
using namespace cl::sycl;
//Typename T, is a floating point type, e.g. float, double
//N is the number of elements in each segment (Regular segment size)
//EPT is the number of elements processed by each thread in a segment
template <typename T, int N>
void segmented_reduce_V0(const sycl::nd_item<1>& nditem, const sycl::local_accessor<T,1>& in, const sycl::local_accessor<T,1>& out, const int nseg){
    auto cta = nditem.get_group();
    int tid = cta.get_local_id(0);
    int bdim = cta.get_local_range(0);
    auto simd = nditem.get_sub_group();
    auto nsimds = bdim / simd.get_local_range()[0];
    auto simdid = simd.get_group_id()[0];
    auto lane = simd.get_local_linear_id();
    for (int seg = simdid; seg < nseg; seg += nsimds){
        T sum = sycl::joint_reduce(simd, in.get_pointer() + seg * N, in.get_pointer() + (seg + 1) * N, T(0), std::plus<T>());
    }
}

template <typename T, int N>
void segmented_reduce_V1(const sycl::nd_item<1>& nditem, const sycl::local_accessor<T,1>& in, const sycl::local_accessor<T,1>& out, const int nseg){
    auto cta = nditem.get_group();
    int tid = cta.get_local_id(0);
    int bdim = cta.get_local_range(0);
    auto simd = nditem.get_sub_group();
    auto nsimds = bdim / simd.get_local_range()[0];
    auto simdid = simd.get_group_id()[0];
    auto lane = simd.get_local_linear_id();
    for (int seg = simdid; seg < nseg; seg += nsimds){
        T sum = T(0);
        for (int i = lane; i < N; i += simd.get_local_range()[0]){
            sum += in[seg * N + i];
        }
        out[seg] = sum;
    }
}

constexpr int pow2(int n){
    return n == 0 ? 1 : 2 * pow2(n - 1);
}

int main(int argc, char** argv){
    constexpr const int N = pow2(10);
    const int nseg = pow2(13); //1024
    const int n = N * nseg;
    const int bdim = pow2(8);
    const int gdim = pow2(9); //512
    std::vector<float> in(n);
    std::vector<float> out(nseg);
    for (int i = 0; i < n; i++) in[i] = 1.0f;
    {
        buffer<float, 1> in_buf(in.data(), n);
        buffer<float, 1> out_buf(out.data(), nseg);
        queue myQueue;
        myQueue.submit([&](handler& h){
            auto in = in_buf.get_access<access::mode::read>(h);
            auto out = out_buf.get_access<access::mode::write>(h);
            local_accessor<float, 1> local_in(range<1>(bdim), h);
            local_accessor<float, 1> local_out(range<1>(bdim), h);
            h.parallel_for<class segmented_reduce>(nd_range<1>(gdim * bdim, bdim), [=](nd_item<1> nditem){
                auto cta = nditem.get_group();
                auto bid = nditem.get_group_linear_id();
                auto event = cta.async_work_group_copy(local_in.get_pointer(), in.get_pointer() + bid * bdim, bdim);
                event.wait();


                segmented_reduce_V0<float, N>(nditem, local_in, local_out, nseg);
            });
        });
    }
    return 0;
}