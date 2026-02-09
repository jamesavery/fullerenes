#include <sycl-headers/sycl-vector.hh>
#include <sycl-headers/sycl-span.hh>
#include <linalg.hh>
#include <numeric>
#include <algorithm>
#include <gtest/gtest.h>


class LapackTest : public ::testing::TestWithParam<int> {
    protected:
        int N = GetParam();
        SyclVector<float> A (N*N);
        std::fill(A.begin(), A.end(), rand() % 100);
        void SetUp() override {
            //BuckyGen::next_fullerene(BQ, G);
        }
    
        void TearDown() override {
            //BuckyGen::stop(BQ);
        }
    };

TEST_P(LapackTest, CholeskyQR) {
    SyclQueue ctx;
    linalg::DenseMatHandle<float, linalg::BatchType::Batched> handleA(A.data(), N, N, N, N*N, 1);
    linalg::ortho_buffer_size<Backend::CUDA>(ctx, handleA, linalg::Transpose::NoTrans, linalg::OrthoAlgorithm::Cholesky);
}