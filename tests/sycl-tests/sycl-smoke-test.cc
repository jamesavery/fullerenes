#include <gtest/gtest.h>
#include <sycl/sycl.hpp>
#include <vector>
#include <numeric>

TEST(SyclSmoke, DeviceDiscovery) {
    auto platforms = sycl::platform::get_platforms();
    ASSERT_FALSE(platforms.empty());
    bool found = false;
    for (auto& p : platforms)
        for (auto& d : p.get_devices())
            if (d.is_cpu() || d.is_gpu()) found = true;
    EXPECT_TRUE(found) << "No CPU or GPU SYCL device found";
}

TEST(SyclSmoke, ParallelAdd) {
    constexpr int N = 1024;
    sycl::queue Q;

    std::vector<float> a(N), b(N), c(N, 0.0f);
    std::iota(a.begin(), a.end(), 0.0f);
    std::iota(b.begin(), b.end(), 1000.0f);

    {
        sycl::buffer<float> buf_a(a.data(), sycl::range<1>(N));
        sycl::buffer<float> buf_b(b.data(), sycl::range<1>(N));
        sycl::buffer<float> buf_c(c.data(), sycl::range<1>(N));

        Q.submit([&](sycl::handler& h) {
            auto acc_a = buf_a.get_access<sycl::access::mode::read>(h);
            auto acc_b = buf_b.get_access<sycl::access::mode::read>(h);
            auto acc_c = buf_c.get_access<sycl::access::mode::write>(h);
            h.parallel_for(sycl::range<1>(N), [=](sycl::id<1> i) {
                acc_c[i] = acc_a[i] + acc_b[i];
            });
        });
    }

    for (int i = 0; i < N; i++)
        EXPECT_FLOAT_EQ(c[i], float(i) + float(1000 + i));
}

TEST(SyclSmoke, Reduction) {
    constexpr int N = 256;
    sycl::queue Q;

    std::vector<float> data(N);
    std::iota(data.begin(), data.end(), 1.0f);
    float expected = N * (N + 1) / 2.0f;

    // USM-pointer reduction: the buffer-based sycl::reduction(buffer, handler, op)
    // overload is not provided by all AdaptiveCpp backends, so use the
    // device-pointer reduction form, which is portable.
    float* d_data = sycl::malloc_shared<float>(N, Q);
    float* d_sum  = sycl::malloc_shared<float>(1, Q);
    std::copy(data.begin(), data.end(), d_data);
    *d_sum = 0.0f;

    Q.submit([&](sycl::handler& h) {
        h.parallel_for(sycl::range<1>(N), sycl::reduction(d_sum, sycl::plus<float>()),
            [=](sycl::id<1> i, auto& s) { s += d_data[i]; });
    }).wait();

    float result = *d_sum;
    sycl::free(d_data, Q);
    sycl::free(d_sum, Q);

    EXPECT_NEAR(result, expected, 0.1f);
}
