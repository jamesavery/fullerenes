// Phase 7 commit 1: verify that the new view-based DualizeFunctor overload
// produces the same cubic adjacency as the legacy FullereneBatchView path.

#include <gtest/gtest.h>
#include <fullerenes/kernel-headers/all-kernels.hh>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/graph.hh>
#include <fullerenes/dense_graph.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>

#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

template<typename K>
class DualizeViewTest : public ::testing::Test {};

TEST(DualizeViewTest, MatchesFullereneBatchPathC60) {
    using T = float;
    using K = uint16_t;

    constexpr int N        = 60;
    constexpr int Nf       = N/2 + 2;      // 32
    constexpr int capacity = 3;

    // --- Reference path: FullereneBatch + legacy compute() ------------------
    FullereneBatch<T,K> ref_batch(N, capacity);
    auto BQ = BuckyGen::start(N, false, false);
    Graph G(Nf, GRAPH_DMAX);
    for (int i = 0; i < capacity; ++i) {
        ASSERT_TRUE(BuckyGen::next_fullerene(BQ, G));
        ref_batch.push_back(G, uint64_t(i));
    }
    BuckyGen::stop(BQ);

    // --- Allocate view-path buffers -----------------------------------------
    using RSR = Spanify::RSRAdjacencyView<K>;
    batch::Batch<RSR> src(Nf, capacity, /*dmax*/6);
    batch::Batch<RSR> dst(N,  capacity, /*dmax*/3);
    batch::BatchState st(capacity);
    SyclVector<std::array<K,6>> faces_cubic_buf(capacity * Nf);
    SyclVector<std::array<K,3>> faces_dual_buf (capacity * N);

    SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);

    DualizeFunctor<T,K> dualize;
    dualize(Q, ref_batch, LaunchPolicy::SYNC);
    Q.wait();

    // Populate src from reference dual graph.
    auto src_spans = src.view_capacity().spans();
    auto& src_adj  = std::get<0>(src_spans); // std::span<K>
    auto& src_deg  = std::get<1>(src_spans); // std::span<uint8_t>

    for (int isomer = 0; isomer < capacity; ++isomer) {
        for (int u = 0; u < Nf; ++u) {
            const auto& nbrs = ref_batch.d_.A_dual_[isomer*Nf + u];
            for (int k = 0; k < 6; ++k) {
                src_adj[isomer*Nf*6 + u*6 + k] = K(nbrs[k]);
            }
            src_deg[isomer*Nf + u] = uint8_t(ref_batch.d_.deg_[isomer*Nf + u]);
        }
        st.push_back(uint64_t(isomer),
                     StatusFlag(int(StatusEnum::DUAL_INITIALIZED)));
    }

    auto src_view = src.view_capacity();
    auto dst_view = dst.view_capacity();
    ASSERT_EQ(src_view.size(), capacity);
    ASSERT_EQ(dst_view.size(), capacity);
    ASSERT_EQ(st.size(),        capacity);

    auto ev = dualize.compute(Q, src_view, dst_view, st.view(),
                              Span<std::array<K,6>>(faces_cubic_buf.data(), faces_cubic_buf.size()),
                              Span<std::array<K,3>>(faces_dual_buf.data(),  faces_dual_buf.size()));
    ev.wait();
    Q.wait();

    // --- Compare cubic adjacency -------------------------------------------
    auto dst_spans = dst.view_capacity().spans();
    auto& dst_adj  = std::get<0>(dst_spans); // std::span<K>

    int mismatches = 0;
    for (int isomer = 0; isomer < capacity; ++isomer) {
        for (int u = 0; u < N; ++u) {
            for (int k = 0; k < 3; ++k) {
                K got = dst_adj[isomer*N*3 + u*3 + k];
                K ref = ref_batch.d_.A_cubic_[isomer*N + u][k];
                if (got != ref) ++mismatches;
            }
        }
        // Status should have FULLERENEGRAPH_PREPARED set.
        EXPECT_TRUE(st.status()[isomer].is_set(StatusEnum::FULLERENEGRAPH_PREPARED))
            << "isomer " << isomer << " status not advanced";
    }
    EXPECT_EQ(mismatches, 0);
}

} // namespace
