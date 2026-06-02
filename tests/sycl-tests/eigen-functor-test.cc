// Eigen functor test: verify that running the view-based pipeline on a
// batch of 2 produces per-slot hessian/eigenvalue blocks that match the
// results of running the same kernels on each isomer individually
// (via BatchView::slice and matching xyz subspans).

#include <gtest/gtest.h>
#include <fullerenes/kernel-headers/all-kernels.hh>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/graph.hh>
#include <fullerenes/dense_graph.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>

namespace {

class EigenTest : public ::testing::TestWithParam<int> {};

TEST_P(EigenTest, BatchSliceConsistency) {
    using T = float;
    using K = uint16_t;
    using RSR = Spanify::RSRAdjacencyView<K>;
    const int N  = GetParam();
    const int Nf = N/2 + 2;
    const int B  = 2;

    SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);

    batch::Batch<RSR>             src_dual (Nf, B, 6);
    batch::Batch<RSR>             dst_cubic(N,  B, 3);
    batch::BatchState             st(B);
    SyclVector<std::array<K,6>>   faces_cubic(B * Nf);
    SyclVector<std::array<K,3>>   faces_dual (B * N);
    SyclVector<std::array<T,2>>   layout2d(B * N);
    SyclVector<std::array<T,3>>   xyz(B * N);

    auto BQ = BuckyGen::start(N, false, false);
    Graph G(Nf, GRAPH_DMAX);
    // Skip a few to get distinct isomers across slots.
    for (int i = 0; i < 4; ++i) BuckyGen::next_fullerene(BQ, G);
    for (int b = 0; b < B; ++b) {
        ASSERT_TRUE(BuckyGen::next_fullerene(BQ, G));
        auto vc  = src_dual.view_capacity();
        auto& adj = std::get<0>(vc.spans());
        auto& deg = std::get<1>(vc.spans());
        for (int u = 0; u < Nf; ++u) {
            int du = G.deg[u];
            for (int k = 0; k < 6; ++k)
                adj[b*Nf*6 + u*6 + k] =
                    K(k < du ? G.neighbours[u*G.dmax + k]
                             : std::numeric_limits<K>::max());
            deg[b*Nf + u] = uint8_t(du);
        }
        st.push_back(uint64_t(b),
                     StatusFlag(int(StatusEnum::DUAL_INITIALIZED)));
    }
    BuckyGen::stop(BQ);
    src_dual.resize(B);

    dualize<T,K>(Q, src_dual.view(),
                    dst_cubic.view_capacity(), st.view(),
                    std::span<std::array<K,6>>(faces_cubic.data(), faces_cubic.size()),
                    std::span<std::array<K,3>>(faces_dual.data(),  faces_dual.size())).wait();
    dst_cubic.resize(B);
    tutte_layout<T,K>(Q, dst_cubic.view(),
                  std::span<std::array<T,2>>(layout2d.data(), layout2d.size()),
                  st.view()).wait();
    spherical_projection<T,K>(Q, dst_cubic.view(),
                std::span<std::array<T,2>>(layout2d.data(), layout2d.size()),
                std::span<std::array<T,3>>(xyz.data(), xyz.size()),
                st.view()).wait();
    forcefield_optimize<PEDERSEN,T,K>(Q, dst_cubic.view(),
               std::span<std::array<T,3>>(xyz.data(), xyz.size()),
               st.view(), 5*N, 5*N).wait();

    // Full-batch hessian + eigen.
    SyclVector<T> hess_b(B * 90 * N);
    SyclVector<K> cols_b(B * 90 * N);
    SyclVector<T> eigs_b(B * 2);
    SyclVector<T> evec_b(B * 2 * N);
    constexpr int n_lanczos = 50;
    const int n_lanczos_full = N*3 - 6;
    const int n_lanczos_max  = std::max(n_lanczos, n_lanczos_full);
    SyclVector<T> off_b (B * n_lanczos_max);
    SyclVector<T> qm_b  (B * n_lanczos_max * n_lanczos_max);
    SyclVector<T> lan_b (B * n_lanczos_max * N * 3);
    SyclVector<T> diag_b(B * n_lanczos_max);
    SyclVector<K> ends_b(B * 2);

    compute_hessian<PEDERSEN,T,K>(Q, dst_cubic.view(),
                    std::span<std::array<T,3>>(xyz.data(), xyz.size()),
                    st.view(),
                    std::span<T>(hess_b.data(), hess_b.size()),
                    std::span<K>(cols_b.data(), cols_b.size())).wait();
    eigensolve<EigensolveMode::ENDS,T,K>(Q, std::span<std::array<T,3>>(xyz.data(), xyz.size()),
                 N, B,
                 std::span<T>(hess_b.data(), hess_b.size()),
                 std::span<K>(cols_b.data(), cols_b.size()),
                 n_lanczos,
                 std::span<T>(eigs_b.data(), eigs_b.size()),
                 std::span<T>(evec_b.data(), evec_b.size()),
                 std::span<T>(off_b.data(), off_b.size()),
                 std::span<T>(qm_b.data(), qm_b.size()),
                 std::span<T>(lan_b.data(), lan_b.size()),
                 std::span<T>(diag_b.data(), diag_b.size()),
                 std::span<K>(ends_b.data(), ends_b.size())).wait();

    SyclVector<T> eigs_full_b(B * 3 * N);
    SyclVector<T> evec_full_b(B * 3 * N * 3 * N);
    eigensolve<EigensolveMode::FULL_SPECTRUM,T,K>(Q, std::span<std::array<T,3>>(xyz.data(), xyz.size()),
                 N, B,
                 std::span<T>(hess_b.data(), hess_b.size()),
                 std::span<K>(cols_b.data(), cols_b.size()),
                 N*3 - 6,
                 std::span<T>(eigs_full_b.data(), eigs_full_b.size()),
                 std::span<T>(evec_full_b.data(), evec_full_b.size()),
                 std::span<T>(off_b.data(), off_b.size()),
                 std::span<T>(qm_b.data(), qm_b.size()),
                 std::span<T>(lan_b.data(), lan_b.size()),
                 std::span<T>(diag_b.data(), diag_b.size()),
                 std::span<K>(ends_b.data(), ends_b.size())).wait();

    // Per-slot reference: rerun on each slice individually.
    auto run_single = [&](int b,
                          SyclVector<T>& hess1,
                          SyclVector<K>& cols1,
                          SyclVector<T>& eigs1,
                          SyclVector<T>& evec1,
                          SyclVector<T>& eigs_f1,
                          SyclVector<T>& evec_f1) {
        // For slice(b,1) the kernels iterate over view.size()==1, so we pass
        // single-isomer-sized output buffers + xyz subspan starting at b*N.
        SyclVector<T> off1 (n_lanczos_max);
        SyclVector<T> qm1  (n_lanczos_max * n_lanczos_max);
        SyclVector<T> lan1 (n_lanczos_max * N * 3);
        SyclVector<T> diag1(n_lanczos_max);
        SyclVector<K> end1 (2);

        std::span<std::array<T,3>> xyz1(xyz.data() + std::size_t(b) * N, N);
        auto bv1 = dst_cubic.view().slice(std::size_t(b), 1);
        auto st1 = st.slice(std::size_t(b), 1);

        compute_hessian<PEDERSEN,T,K>(Q, bv1, xyz1, st1,
                        std::span<T>(hess1.data(), hess1.size()),
                        std::span<K>(cols1.data(), cols1.size())).wait();
        eigensolve<EigensolveMode::ENDS,T,K>(Q, xyz1, N, 1,
                     std::span<T>(hess1.data(), hess1.size()),
                     std::span<K>(cols1.data(), cols1.size()),
                     n_lanczos,
                     std::span<T>(eigs1.data(), eigs1.size()),
                     std::span<T>(evec1.data(), evec1.size()),
                     std::span<T>(off1.data(), off1.size()),
                     std::span<T>(qm1.data(), qm1.size()),
                     std::span<T>(lan1.data(), lan1.size()),
                     std::span<T>(diag1.data(), diag1.size()),
                     std::span<K>(end1.data(), end1.size())).wait();
        eigensolve<EigensolveMode::FULL_SPECTRUM,T,K>(Q, xyz1, N, 1,
                     std::span<T>(hess1.data(), hess1.size()),
                     std::span<K>(cols1.data(), cols1.size()),
                     N*3 - 6,
                     std::span<T>(eigs_f1.data(), eigs_f1.size()),
                     std::span<T>(evec_f1.data(), evec_f1.size()),
                     std::span<T>(off1.data(), off1.size()),
                     std::span<T>(qm1.data(), qm1.size()),
                     std::span<T>(lan1.data(), lan1.size()),
                     std::span<T>(diag1.data(), diag1.size()),
                     std::span<K>(end1.data(), end1.size())).wait();
    };

    SyclVector<T> hess0 (90 * N), hess1 (90 * N);
    SyclVector<K> cols0 (90 * N), cols1_ (90 * N);
    SyclVector<T> eigs0 (2),       eigs1__(2);
    SyclVector<T> evec0 (2 * N),   evec1__(2 * N);
    SyclVector<T> eigsf0(3 * N),   eigsf1 (3 * N);
    SyclVector<T> evecf0(3 * N * 3 * N), evecf1(3 * N * 3 * N);
    run_single(0, hess0, cols0, eigs0, evec0, eigsf0, evecf0);
    run_single(1, hess1, cols1_, eigs1__, evec1__, eigsf1, evecf1);

    auto eq = [](const SyclVector<T>& a, std::span<const T> b, T tol) {
        if (a.size() != b.size()) return false;
        for (std::size_t i = 0; i < a.size(); ++i) {
            if (!(std::abs(a[i] - b[i]) <= tol * (std::abs(a[i]) + std::abs(b[i]) + T(1))))
                return false;
        }
        return true;
    };
    EXPECT_TRUE(eq(hess0,  std::span<const T>(hess_b.data(),         90*N), T(1e-3)));
    EXPECT_TRUE(eq(hess1,  std::span<const T>(hess_b.data() + 90*N,  90*N), T(1e-3)));
    EXPECT_TRUE(eq(eigs0,  std::span<const T>(eigs_b.data(),         2),    T(1e-3)));
    EXPECT_TRUE(eq(eigs1__, std::span<const T>(eigs_b.data() + 2,    2),    T(1e-3)));
    EXPECT_TRUE(eq(eigsf0, std::span<const T>(eigs_full_b.data(),         3*N), T(1e-3)));
    EXPECT_TRUE(eq(eigsf1, std::span<const T>(eigs_full_b.data() + 3*N,   3*N), T(1e-3)));
}

INSTANTIATE_TEST_SUITE_P(, EigenTest, ::testing::Values(60));

} // namespace

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
