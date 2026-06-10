// Tests for the new low-level batch types: batch::Batch<V>,
// batch::BatchState, batch::BatchQueue<V> and the queue_push transfer
// helpers. Replaces the legacy FullereneBatch/FullereneQueue test suite.

#include <gtest/gtest.h>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/graph.hh>
#include <fullerenes/dense_graph.hh>
#include <fullerenes/batch/batch.hh>
#include <fullerenes/batch/batch_state.hh>
#include <fullerenes/batch/batch_queue.hh>
#include <fullerenes/batch/utilities.hh>
#include <fullerenes/sycl-headers/sycl-status-enum.hh>

#include <array>
#include <cstdint>
#include <cstdlib>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>

// Initialize SYCL runtime early so that Batch/Queue (which uses SyclVector
// storage when FULLERENES_ENABLE_SYCL is defined) lazy-init does not leave
// subsystems half-brought-up, which otherwise causes a hang at global
// destructor teardown in this test binary.
#include <fullerenes/kernel-headers/all-kernels.hh>

namespace {

using K = uint16_t;
using RSR = Spanify::RSRAdjacencyView<K>;

class BatchFixture : public ::testing::TestWithParam<int> {
protected:
    int N  = GetParam();
    int Nf = N/2 + 2;
    static Graph& shared_graph(int N) {
        static Graph G(N/2 + 2, GRAPH_DMAX);
        static bool initialized = false;
        if (!initialized) {
            auto BQ = BuckyGen::start(N, false, false);
            BuckyGen::next_fullerene(BQ, G);
            BuckyGen::stop(BQ);
            initialized = true;
        }
        return G;
    }
    Graph& G = shared_graph(N);
};

TEST_P(BatchFixture, BatchConstructAndPush) {
    batch::Batch<RSR> b(Nf, /*capacity*/3, /*dmax*/6);
    EXPECT_EQ(b.size(), 0);
    EXPECT_EQ(b.capacity(), 3);
    EXPECT_EQ(b.N(), Nf);
    EXPECT_EQ(b.dmax(), 6);

    batch::push_back(b, G);
    EXPECT_EQ(b.size(), 1);

    auto v = b.view()[0];
    auto& adj = std::get<0>(v.to_tuple());
    EXPECT_EQ(int(adj[0]), G.nbrs(0)[0]);
}

TEST_P(BatchFixture, BatchClear) {
    batch::Batch<RSR> b(Nf, 5, 6);
    batch::push_back(b, G);
    batch::push_back(b, G);
    EXPECT_EQ(b.size(), 2);
    b.clear();
    EXPECT_EQ(b.size(), 0);
    EXPECT_EQ(b.capacity(), 5);
}

TEST_P(BatchFixture, BatchStateBasics) {
    batch::BatchState s(4);
    EXPECT_EQ(s.size(), 0);
    EXPECT_EQ(s.capacity(), 4);
    s.push_back(7, StatusFlag(int(StatusEnum::DUAL_INITIALIZED)), 3, -1);
    EXPECT_EQ(s.size(), 1);
    auto sv = s.view();
    EXPECT_EQ(int(sv.id[0]), 7);
    EXPECT_EQ(int(sv.iteration[0]), 3);
    EXPECT_TRUE(sv.status[0].is_set(StatusEnum::DUAL_INITIALIZED));

    s.write_slot(2, 99,
                 StatusFlag(int(StatusEnum::FULLERENEGRAPH_PREPARED)), 11);
    uint64_t id; StatusFlag f; int32_t it;
    s.read_slot(2, id, f, it);
    EXPECT_EQ(int(id), 99);
    EXPECT_EQ(int(it), 11);
    EXPECT_TRUE(f.is_set(StatusEnum::FULLERENEGRAPH_PREPARED));
}

TEST_P(BatchFixture, BatchQueuePushPopOrder) {
    batch::BatchQueue<RSR> q(Nf, /*capacity*/4, /*dmax*/6);
    EXPECT_TRUE(q.empty());

    for (int i = 0; i < 4; ++i)
        batch::push_back(q, G, uint64_t(100 + i),
                         StatusFlag(int(StatusEnum::DUAL_INITIALIZED)),
                         int32_t(i));
    EXPECT_EQ(q.size(), 4);
    EXPECT_EQ(q.front_index(), 0);
    EXPECT_EQ(q.back_index(),  3);

    // Pop two and check ordering.
    for (int i = 0; i < 2; ++i) {
        RSR out_v = q.storage().view_capacity()[0];  // dummy stamp
        // For pop, we need an addressable RSR; use the storage to give us
        // properly-shaped spans (pop_front does an internal copy_fields).
        // Simpler: just discard the first two and check IDs of remaining.
        uint64_t id; StatusFlag f; int32_t it;
        q.state().read_slot(q.front_index(), id, f, it);
        EXPECT_EQ(int(id), 100 + i);
        q.discard_front();
    }
    EXPECT_EQ(q.size(), 2);
    EXPECT_EQ(q.front_index(), 2);
}

TEST_P(BatchFixture, BatchQueueGrowsCircularly) {
    batch::BatchQueue<RSR> q(Nf, /*capacity*/0, /*dmax*/6);

    for (int i = 0; i < 5; ++i)
        batch::push_back(q, G, uint64_t(i),
                         StatusFlag(int(StatusEnum::DUAL_INITIALIZED)));
    EXPECT_EQ(q.size(), 5);
    EXPECT_GE(q.capacity(), 5);
    // After growth, logical front == 0.
    EXPECT_EQ(q.front_index(), 0);
    for (int i = 0; i < 5; ++i) {
        uint64_t id; StatusFlag f; int32_t it;
        q.state().read_slot(q.slot_of(i), id, f, it);
        EXPECT_EQ(int(id), i);
    }
}

TEST_P(BatchFixture, QueuePushBatchToQueue) {
    batch::Batch<RSR> b(Nf, 3, 6);
    batch::BatchState s(3);
    for (int i = 0; i < 3; ++i) {
        batch::push_back(b, G);
        s.push_back(uint64_t(i),
                    StatusFlag(int(StatusEnum::DUAL_INITIALIZED)));
    }
    batch::BatchQueue<RSR> q(Nf, 3, 6);
    int n = batch::queue_push(q, b, s,
        batch::StatusPredicate{StatusEnum::DUAL_INITIALIZED});
    EXPECT_EQ(n, 3);
    EXPECT_EQ(q.size(), 3);
}

TEST_P(BatchFixture, QueuePushQueueToBatchWithPredicate) {
    batch::BatchQueue<RSR> q(Nf, 5, 6);
    for (int i = 0; i < 5; ++i)
        batch::push_back(q, G, uint64_t(i),
                         StatusFlag(i % 2 == 0 ? int(StatusEnum::DUAL_INITIALIZED)
                                               : int(StatusEnum::EMPTY)));
    batch::Batch<RSR> dst_b(Nf, 5, 6);
    batch::BatchState dst_s(5);
    int n = batch::queue_push(dst_b, dst_s, q,
        batch::StatusPredicate{StatusEnum::DUAL_INITIALIZED, StatusEnum::EMPTY});
    EXPECT_EQ(n, 1);    // only the front (id=0) matches before stopping.
    EXPECT_EQ(dst_b.size(), 1);
    EXPECT_EQ(int(dst_s.view().id[0]), 0);
    EXPECT_EQ(q.size(), 4);
}

TEST_P(BatchFixture, QueuePushSetsConsumedFlag) {
    batch::Batch<RSR> b(Nf, 2, 6);
    batch::BatchState s(2);
    for (int i = 0; i < 2; ++i) {
        batch::push_back(b, G);
        s.push_back(uint64_t(i),
                    StatusFlag(int(StatusEnum::DUAL_INITIALIZED)));
    }
    batch::BatchQueue<RSR> q(Nf, 2, 6);
    batch::queue_push(q, b, s, batch::MatchAll{}, StatusEnum::EMPTY);
    auto sv = s.view();
    EXPECT_TRUE(sv.status[0].is_set(StatusEnum::EMPTY));
    EXPECT_TRUE(sv.status[1].is_set(StatusEnum::EMPTY));
}

INSTANTIATE_TEST_SUITE_P(, BatchFixture, ::testing::Values(60));

} // namespace

int main(int argc, char** argv) {
    // Put ourselves into our own process group so that any child processes
    // spawned (e.g. AdaptiveCpp helper daemons, BuckyGen workers) can be
    // reaped together at exit. Without this, child processes inherit the
    // parent's stdout/stderr pipe and block the parent's controlling
    // process (e.g. ctest) from observing EOF until they exit on their own.
    setpgid(0, 0);

    // Force SYCL runtime initialization up front.
    { SyclQueue Q{Device::get_devices(DeviceType::GPU).at(0), true}; (void)Q; }
    testing::InitGoogleTest(&argc, argv);
    int rc = RUN_ALL_TESTS();

    // Flush before we forcibly tear down.
    std::cout.flush();
    std::cerr.flush();
    // Close stdio so any descriptors inherited by child processes that we
    // cannot directly reap do not keep the parent's stdout pipe alive.
    ::close(0); ::close(1); ::close(2);
    // Ignore SIGTERM so we survive the broadcast below, then kill the rest
    // of our process group.
    ::signal(SIGTERM, SIG_IGN);
    ::killpg(0, SIGTERM);
    std::_Exit(rc);
}
