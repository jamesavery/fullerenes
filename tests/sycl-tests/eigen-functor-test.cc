#include <gtest/gtest.h>
#include <fullerenes/kernel-headers/all-kernels.hh>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/polyhedron.hh>
#include <fullerenes/planargraph.hh>
#include <fullerenes/spiral.hh>
#include <fullerenes/kernel-headers/sycl-parallel-primitives.hh>
#include <iostream>
#include <algorithm>
#include <thread>
#include <chrono>
#include <future>

class MinimumProblem : public ::testing::TestWithParam<int> {
protected:
    int N = GetParam();
    using T = float;
    using K = uint16_t;
    BuckyGen::buckygen_queue BQ = BuckyGen::start(N, false, false);
    FullereneBatch<T, K> batch_1;
    FullereneBatch<T, K> batch_2;
    Graph G(N/2 + 2, GRAPH_DMAX);
    
    
    void SetUp() override {
        BuckyGen::next_fullerene(BQ, G);
        BuckyGen::next_fullerene(BQ, G);
        BuckyGen::next_fullerene(BQ, G);
        BuckyGen::next_fullerene(BQ, G);
        batch_1.push_back(G, 0);
        BuckyGen::next_fullerene(BQ, G);
        batch_1.push_back(G, 1);
        BuckyGen::next_fullerene(BQ, G);
        batch_2.push_back(G, 0);
        BuckyGen::next_fullerene(BQ, G);
        batch_2.push_back(G, 1);

    }

    void TearDown() override {
        BuckyGen::stop(BQ);
    }
};


TEST_P(MinimumProblem, AllTestsInOne) {
    {
        auto batch1 = batch_1;
        auto batch2 = batch_2;
        DualizeFunctor<T, uint16_t> dualize;
        TutteFunctor<T, uint16_t> tutte;
        SphericalProjectionFunctor<T, uint16_t> spherical_projection;
        ForcefieldOptimizeFunctor<PEDERSEN, T, uint16_t> forcefield_optimize;
        HessianFunctor<PEDERSEN, T, uint16_t> compute_hessian;
        EigenFunctor<EigensolveMode::ENDS, T, uint16_t> spectrum_ends;

        SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);

        dualize(Q, batch1, LaunchPolicy::SYNC);
        tutte(Q, batch1, LaunchPolicy::SYNC);
        spherical_projection(Q, batch1, LaunchPolicy::SYNC);
        forcefield_optimize(Q, batch1, LaunchPolicy::SYNC, 5*N, 5*N);

        FullereneBatch<T, uint16_t> batch1_copy = batch1;
        SyclVector<T> fullerene1_hessian(90*N);
        SyclVector<K> fullerene1_cols(90*N);
        SyclVector<T> fullerene1_eigs(2);
        SyclVector<T> fullerene1_eigvecs(2*N);

        SyclVector<T> fullerene2_hessian(90*N);
        SyclVector<K> fullerene2_cols(90*N);
        SyclVector<T> fullerene2_eigs(2);
        SyclVector<T> fullerene2_eigvecs(2*N);

        SyclVector<T> batch1_hessians(90*N * batch1.size());
        SyclVector<K> batch1_cols(90*N * batch1.size());
        SyclVector<T> batch1_eigenvalues(2 * batch1.size());
        SyclVector<T> batch1_eigenvectors(2*N * batch1.size());

        compute_hessian(Q, batch1, LaunchPolicy::SYNC , batch1_hessians, batch1_cols);
        compute_hessian(Q, FullereneBatchView(batch1_copy, 0, 1), LaunchPolicy::SYNC , fullerene1_hessian, fullerene1_cols);
        compute_hessian(Q, FullereneBatchView(batch1_copy, 1, 1), LaunchPolicy::SYNC , fullerene2_hessian, fullerene2_cols);
        ASSERT_EQ(batch1_hessians.subspan(0, 90*N), fullerene1_hessian.to_span());
        ASSERT_EQ(batch1_hessians.subspan(90*N, 90*N), fullerene2_hessian.to_span());

        spectrum_ends(Q, batch1, LaunchPolicy::SYNC, batch1_hessians, batch1_cols, 50, batch1_eigenvalues, batch1_eigenvectors);
        spectrum_ends(Q, FullereneBatchView(batch1_copy, 0, 1), LaunchPolicy::SYNC, fullerene1_hessian, fullerene1_cols, 50, fullerene1_eigs, fullerene1_eigvecs);
        spectrum_ends(Q, FullereneBatchView(batch1_copy, 1, 1), LaunchPolicy::SYNC, fullerene2_hessian, fullerene2_cols, 50, fullerene2_eigs, fullerene2_eigvecs);

        ASSERT_EQ(batch1_eigenvalues.subspan(0, 2), fullerene1_eigs.to_span());
        ASSERT_EQ(batch1_eigenvalues.subspan(2, 2), fullerene2_eigs.to_span());

        EigenFunctor<EigensolveMode::FULL_SPECTRUM, T, uint16_t> spectrum_full;
        SyclVector<T> fullerene1_eigs_full(N*3);
        SyclVector<T> fullerene1_eigvecs_full(N*3*N*3);
        SyclVector<T> fullerene2_eigs_full(N*3);
        SyclVector<T> fullerene2_eigvecs_full(N*3*N*3);


        SyclVector<T> batch1_eigenvalues_full(3 * N * batch1.size());
        SyclVector<T> batch1_eigenvectors_full(3*N * 3 * N * batch1.size());
        spectrum_full(Q, batch1, LaunchPolicy::SYNC, batch1_hessians, batch1_cols, N*3-6, batch1_eigenvalues_full, batch1_eigenvectors_full);
        spectrum_full(Q, FullereneBatchView(batch1_copy, 0, 1), LaunchPolicy::SYNC, fullerene1_hessian, fullerene1_cols, N*3-6, fullerene1_eigs_full, fullerene1_eigvecs_full);
        spectrum_full(Q, FullereneBatchView(batch1_copy, 1, 1), LaunchPolicy::SYNC, fullerene2_hessian, fullerene2_cols, N*3-6, fullerene2_eigs_full, fullerene2_eigvecs_full);

        ASSERT_EQ(batch1_eigenvalues_full.subspan(0, 3*N), fullerene1_eigs_full.to_span());
        ASSERT_EQ(batch1_eigenvalues_full.subspan(3*N, 3*N), fullerene2_eigs_full.to_span());

        std::sort(fullerene1_eigs_full.begin(), fullerene1_eigs_full.end());
        std::sort(fullerene2_eigs_full.begin(), fullerene2_eigs_full.end());
        std::cout << "Eigenvalues: " << fullerene1_eigs_full << std::endl;
        std::cout << "Eigenvalues: " << fullerene2_eigs_full << std::endl;
    }
}

INSTANTIATE_TEST_SUITE_P(, MinimumProblem, ::testing::Range(60, 61, 20)); 

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}