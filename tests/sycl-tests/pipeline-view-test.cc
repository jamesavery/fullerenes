// Phase 7 pipeline test: verify the view-based functor overloads can be chained
// as a full geometry pipeline:
//   DualizeFunctor  → cubic adjacency in BatchView
//   TutteFunctor    → 2D layout
//   SphericalProjectionFunctor -> 3D coordinates
//
// This test intentionally stops after spherical projection. Forcefield
// optimisation behaviour is validated separately.

#include <gtest/gtest.h>
#include <fullerenes/kernel-headers/all-kernels.hh>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/graph.hh>
#include <fullerenes/dense_graph.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>

#include <array>
#include <cstdint>
#include <cmath>
#include <numeric>
#include <iostream>

namespace {

TEST(PipelineViewTest, C60ViewPipelineToSphericalProjection) {
    using T = float;
    using K = uint16_t;

    constexpr int N        = 60;
    constexpr int Nf       = N/2 + 2;   // 32
    constexpr int capacity = 5;

    SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);

    // -----------------------------------------------------------------
    // Step 0: generate dual graphs from BuckyGen
    // -----------------------------------------------------------------
    using RSR = Spanify::RSRAdjacencyView<K>;
    batch::Batch<RSR> src_dual(Nf, capacity, /*dmax*/6);
    batch::Batch<RSR> dst_cubic(N, capacity, /*dmax*/3);
    batch::BatchState st(capacity);
    SyclVector<std::array<K,6>> faces_cubic_buf(capacity * Nf);
    SyclVector<std::array<K,3>> faces_dual_buf (capacity * N);

    {
        // Generate dual graphs directly into src_dual.
        auto BQ = BuckyGen::start(N, false, false);
        Graph G(Nf, GRAPH_DMAX);
        auto src_spans = src_dual.view_capacity().spans();
        auto& src_adj  = std::get<0>(src_spans);
        auto& src_deg  = std::get<1>(src_spans);
        for (int b = 0; b < capacity; ++b) {
            ASSERT_TRUE(BuckyGen::next_fullerene(BQ, G));
            for (int u = 0; u < Nf; ++u) {
                int du = G.deg[u];
                for (int k = 0; k < 6; ++k)
                    src_adj[b*Nf*6 + u*6 + k] =
                        K(k < du ? G.neighbours[u*G.dmax + k]
                                 : std::numeric_limits<K>::max());
                src_deg[b*Nf + u] = uint8_t(du);
            }
            st.push_back(uint64_t(b), StatusFlag(int(StatusEnum::DUAL_INITIALIZED)));
        }
        BuckyGen::stop(BQ);
    }

    // -----------------------------------------------------------------
    // Step 1: dualize — dual → cubic adjacency
    // -----------------------------------------------------------------
    dualize<T,K>(Q,
        src_dual.view_capacity(), dst_cubic.view_capacity(), st.view(),
        std::span<std::array<K,6>>(faces_cubic_buf.data(), faces_cubic_buf.size()),
        std::span<std::array<K,3>>(faces_dual_buf.data(),  faces_dual_buf.size())
    ).wait();

    // -----------------------------------------------------------------
    // Step 2: tutte_layout — cubic adjacency → 2D layout
    // -----------------------------------------------------------------
    SyclVector<std::array<T,2>> layout2d(capacity * N);
    tutte_layout<T,K>(Q,
        dst_cubic.view_capacity(),
        std::span<std::array<T,2>>(layout2d.data(), layout2d.size()),
        st.view()
    ).wait();

    // Verify 2D layout is non-zero for each isomer.
    for (int b = 0; b < capacity; ++b) {
        double sum = 0.0;
        for (int u = 0; u < N; ++u) {
            sum += std::abs(layout2d[b*N + u][0]);
            sum += std::abs(layout2d[b*N + u][1]);
        }
        EXPECT_GT(sum, 0.0) << "Tutte layout all-zero for isomer " << b;
    }

    // -----------------------------------------------------------------
    // Step 3: spherical_projection — 2D layout → 3D coords
    // -----------------------------------------------------------------
    SyclVector<std::array<T,3>> xyz(capacity * N);
    spherical_projection<T,K>(Q,
        dst_cubic.view_capacity(),
        std::span<std::array<T,2>>(layout2d.data(), layout2d.size()),
        std::span<std::array<T,3>>(xyz.data(),      xyz.size()),
        st.view()
    ).wait();

    // After spherical projection coords are scaled (not unit sphere) — just
    // verify they are finite and non-zero.
    for (int b = 0; b < capacity; ++b) {
        double sum = 0.0;
        for (int u = 0; u < N; ++u) {
            const auto& p = xyz[b*N + u];
            EXPECT_TRUE(std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]))
                << "Non-finite coord after SphericalProjection: isomer " << b << " node " << u;
            sum += std::abs(p[0]) + std::abs(p[1]) + std::abs(p[2]);
        }
        EXPECT_GT(sum, 0.0) << "All-zero xyz after SphericalProjection: isomer " << b;
    }

    // Pipeline coverage for Phase 7 ends at spherical projection.
}

} // namespace
