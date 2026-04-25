// Forcefield determinism test: N independent BatchQueues carrying
// (graph + xyz) atomically via RSRPolyhedronView must produce identical
// xyz output after running the optimizer in lockstep over multiple rounds.

#include <gtest/gtest.h>
#include <fullerenes/kernel-headers/all-kernels.hh>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/graph.hh>
#include <fullerenes/dense_graph.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/batch/batch_queue.hh>
#include <fullerenes/sycl-headers/sycl-vector.hh>

#include <array>
#include <cstdint>
#include <cmath>
#include <numeric>
#include <vector>

namespace {

class ForceFieldTest : public ::testing::TestWithParam<int> {};

TEST_P(ForceFieldTest, DeterministicAcrossParallelQueues) {
    using T    = float;
    using K    = uint16_t;
    using Adj  = Spanify::RSRAdjacencyView<K>;
    using Poly = Spanify::RSRPolyhedronView<T,K>;

    const int N  = GetParam();
    const int Nf = N/2 + 2;

    SyclQueue Q(Device::get_devices(DeviceType::GPU).at(0), true);

    auto BQ = BuckyGen::start(N, false, false);
    auto device = Device::get_devices(DeviceType::GPU)[0];
    int  cu     = device.get_property(DeviceProperty::MAX_COMPUTE_UNITS);
    int  M      = std::max(2, cu * 2);  // entries per queue

    // -----------------------------------------------------------------
    // Seed: dualize + tutte + sph on M isomers to produce (cubic graph,
    // xyz) pairs that populate the queues.
    // -----------------------------------------------------------------
    batch::Batch<Adj>             src_dual (Nf, M, 6);
    batch::Batch<Poly>            seed_poly(N,  M, 3);
    batch::BatchState             st(M);
    SyclVector<std::array<K,6>>   faces_cubic(M * Nf);
    SyclVector<std::array<K,3>>   faces_dual (M * N);
    SyclVector<std::array<T,2>>   layout2d(M * N);

    Graph G(Nf, GRAPH_DMAX);
    for (int i = 0; i < M; ++i) {
        ASSERT_TRUE(BuckyGen::next_fullerene(BQ, G));
        auto vc  = src_dual.view_capacity();
        auto& adj = std::get<0>(vc.spans());
        auto& deg = std::get<1>(vc.spans());
        for (int u = 0; u < Nf; ++u) {
            int du = G.degree(u);
            for (int k = 0; k < 6; ++k)
                adj[i*Nf*6 + u*6 + k] =
                    K(k < du ? G.nbrs(u)[k]
                              : std::numeric_limits<K>::max());
            deg[i*Nf + u] = uint8_t(du);
        }
        st.push_back(uint64_t(i),
                     StatusFlag(int(StatusEnum::DUAL_INITIALIZED)));
    }
    BuckyGen::stop(BQ);
    src_dual.resize(M);

    DualizeFunctor<T,K>             dualize;
    TutteFunctor<T,K>               tutte;
    SphericalProjectionFunctor<T,K> sph;
    ForcefieldOptimizeFunctor<PEDERSEN, T, K> ff;

    // Dualize into the adjacency prefix of seed_poly; points stay unused
    // until sph writes them.
    auto seed_adj = Spanify::as_adjacency_view(seed_poly.view_capacity());
    dualize.compute(Q, src_dual.view(), seed_adj, st.view(),
                    std::span<std::array<K,6>>(faces_cubic.data(), faces_cubic.size()),
                    std::span<std::array<K,3>>(faces_dual.data(),  faces_dual.size())).wait();
    seed_poly.resize(M);
    auto seed_adj_sz = Spanify::as_adjacency_view(seed_poly.view());
    auto seed_xyz    = Spanify::points_span  (seed_poly.view());
    tutte.compute(Q, seed_adj_sz,
                  std::span<std::array<T,2>>(layout2d.data(), layout2d.size()),
                  st.view()).wait();
    sph.compute(Q, seed_adj_sz,
                std::span<std::array<T,2>>(layout2d.data(), layout2d.size()),
                seed_xyz,
                st.view()).wait();

    // -----------------------------------------------------------------
    // Build n_queues input queues, each holding M copies of the seed.
    // BatchQueue<Poly> stores graph + xyz atomically.
    // -----------------------------------------------------------------
    constexpr int n_queues = 10;
    std::vector<batch::BatchQueue<Poly>> in_q;
    std::vector<batch::BatchQueue<Poly>> out_q;
    in_q.reserve(n_queues); out_q.reserve(n_queues);
    for (int j = 0; j < n_queues; ++j) {
        in_q.emplace_back(N, M, 3);
        out_q.emplace_back(N, M, 3);
        auto sv = seed_poly.view();
        for (int i = 0; i < M; ++i) {
            in_q[j].push_back(sv[std::size_t(i)],
                              uint64_t(i),
                              StatusFlag(int(StatusEnum::FULLERENEGRAPH_PREPARED)),
                              0);
        }
    }

    // Per-round scratch batch + state (one per queue).
    std::vector<batch::Batch<Poly>>  opt_b;
    std::vector<batch::BatchState>   opt_s;
    opt_b.reserve(n_queues); opt_s.reserve(n_queues);
    for (int j = 0; j < n_queues; ++j) {
        opt_b.emplace_back(N, M, 3);
        opt_s.emplace_back(M);
    }

    for (int round = 0; round < 5; ++round) {
        for (int j = 0; j < n_queues; ++j) {
            opt_b[j].clear(); opt_s[j].clear();

            // Queue -> Batch: atomic (graph + xyz) transfer.
            batch::queue_push(opt_b[j], opt_s[j], in_q[j]);

            // Run FF on the adjacency prefix + xyz points of the batch.
            auto bv_adj = Spanify::as_adjacency_view(opt_b[j].view());
            auto bv_xyz = Spanify::points_span  (opt_b[j].view());
            ff.compute(Q, bv_adj, bv_xyz, opt_s[j].view(), N, 10*N).wait();

            // Batch -> Queue: atomic (graph + xyz) transfer back.
            batch::queue_push(out_q[j], opt_b[j], opt_s[j]);
        }

        // All queues must remain identical.
        for (int j = 1; j < n_queues; ++j) {
            ASSERT_EQ(out_q[j].size(), out_q[0].size())
                << "queue size diverged at round " << round << " queue " << j;
            for (int i = 0; i < out_q[0].size(); ++i) {
                auto a = out_q[0].at(i);
                auto b = out_q[j].at(i);
                for (int u = 0; u < N; ++u)
                    for (int c = 0; c < 3; ++c) {
                        ASSERT_FLOAT_EQ(b.points[u][c], a.points[u][c])
                            << "xyz mismatch round=" << round
                            << " q=" << j << " u=" << u << " c=" << c;
                    }
            }
        }
    }
}

INSTANTIATE_TEST_SUITE_P(, ForceFieldTest, ::testing::Values(60));

} // namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
