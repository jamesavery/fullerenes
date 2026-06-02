// Phase-9 smoke test: drive dualize -> tutte -> spherical_projection
// through a single Batch<RSRPolyhedronView> (graph + xyz carried
// atomically). Adjacency-only kernels receive Spanify::as_adjacency_view
// and geometry kernels receive Spanify::points_span. This verifies the
// Poly view slicing helpers work end-to-end without an external xyz
// scratch buffer.
//
// Forcefield / Hessian / Eigen / shape-property kernels are exercised
// by forcefield-functor-test and eigen-functor-test; we do not re-run
// them here.

#include <gtest/gtest.h>
#include <fullerenes/kernel-headers/all-kernels.hh>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/graph.hh>
#include <fullerenes/dense_graph.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>

#include <array>
#include <cmath>
#include <cstdint>

namespace {

class FunctorTests : public ::testing::TestWithParam<int> {};

TEST_P(FunctorTests, DualizeTutteSphThroughPolyView) {
    using T    = float;
    using K    = uint16_t;
    using Adj  = Spanify::RSRAdjacencyView<K>;
    using Poly = Spanify::RSRPolyhedronView<T,K>;

    const int N  = GetParam();
    const int Nf = N/2 + 2;

    SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);
    auto device = Device::get_devices(DeviceType::GPU)[0];
    int  cu     = device.get_property(DeviceProperty::MAX_COMPUTE_UNITS);
    const int capacity = std::max(2, cu * 2);

    batch::Batch<Adj>             src_dual (Nf, capacity, /*dmax*/6);
    batch::Batch<Poly>            dst_poly (N,  capacity, /*dmax*/3);
    batch::BatchState             st(capacity);
    SyclVector<std::array<K,6>>   faces_cubic(capacity * Nf);
    SyclVector<std::array<K,3>>   faces_dual (capacity * N);
    SyclVector<std::array<T,2>>   layout2d   (capacity * N);

    // Fill src_dual from BuckyGen (one distinct isomer per slot).
    auto BQ = BuckyGen::start(N, false, false);
    Graph G(Nf, GRAPH_DMAX);
    {
        auto vc   = src_dual.view_capacity();
        auto& adj = std::get<0>(vc.spans());
        auto& deg = std::get<1>(vc.spans());
        for (int b = 0; b < capacity; ++b) {
            ASSERT_TRUE(BuckyGen::next_fullerene(BQ, G));
            for (int u = 0; u < Nf; ++u) {
                int du = G.degree(u);
                for (int k = 0; k < 6; ++k)
                    adj[b*Nf*6 + u*6 + k] =
                        K(k < du ? G.nbrs(u)[k]
                                  : std::numeric_limits<K>::max());
                deg[b*Nf + u] = uint8_t(du);
            }
            st.push_back(uint64_t(b),
                         StatusFlag(int(StatusEnum::DUAL_INITIALIZED)));
        }
        src_dual.resize(capacity);
    }
    BuckyGen::stop(BQ);

    // Dualize writes into the adjacency prefix of dst_poly.
    {
        auto dst_adj_cap = Spanify::as_adjacency_view(dst_poly.view_capacity());
        dualize<T,K>(Q, src_dual.view_capacity(), dst_adj_cap, st.view(),
                        std::span<std::array<K,6>>(faces_cubic.data(), faces_cubic.size()),
                        std::span<std::array<K,3>>(faces_dual.data(),  faces_dual.size())).wait();
    }
    dst_poly.resize(capacity);

    auto dst_adj = Spanify::as_adjacency_view(dst_poly.view());
    auto dst_xyz = Spanify::points_span  (dst_poly.view());

    tutte_layout<T,K>(Q, dst_adj,
                  std::span<std::array<T,2>>(layout2d.data(), layout2d.size()),
                  st.view()).wait();

    // After tutte the 2D layout must be finite for every vertex.
    for (int b = 0; b < capacity; ++b)
        for (int u = 0; u < N; ++u)
            for (int c = 0; c < 2; ++c)
                ASSERT_TRUE(std::isfinite(layout2d[b*N + u][c]))
                    << "non-finite tutte b=" << b << " u=" << u << " c=" << c;

    spherical_projection<T,K>(Q, dst_adj,
                std::span<std::array<T,2>>(layout2d.data(), layout2d.size()),
                dst_xyz,
                st.view()).wait();

    // After sph the xyz stored in the Poly view must be finite and
    // non-trivial on every slot.
    for (int b = 0; b < capacity; ++b) {
        T r2_max = 0;
        for (int u = 0; u < N; ++u) {
            for (int c = 0; c < 3; ++c)
                ASSERT_TRUE(std::isfinite(dst_xyz[b*N + u][c]))
                    << "non-finite sph b=" << b << " u=" << u << " c=" << c;
            T r2 = 0;
            for (int c = 0; c < 3; ++c) r2 += dst_xyz[b*N + u][c] * dst_xyz[b*N + u][c];
            if (r2 > r2_max) r2_max = r2;
        }
        EXPECT_GT(r2_max, T(0)) << "zero-magnitude xyz on slot " << b;
    }
}

INSTANTIATE_TEST_SUITE_P(, FunctorTests, ::testing::Values(60));

} // namespace

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
