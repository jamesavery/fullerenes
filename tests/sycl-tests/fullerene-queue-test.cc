#include <fullerenes/sycl-headers/all-kernels.hh>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/isomerdb.hh>
#include <numeric>
#include <gtest/gtest.h>


class PushBackTest : public ::testing::TestWithParam<int> {
protected:
    using Fullerene = Fullerene<float, uint16_t>;
    using FullereneBatch = FullereneBatch<float, uint16_t>;
    using FullereneQueue = FullereneQueue<float, uint16_t>;
    using FullereneBatchView = FullereneBatchView<float, uint16_t>;
    
    int N = GetParam();
    int n_max_isomers = IsomerDB::number_isomers(N);
    Graph G = Graph(neighbours_t(N/2 + 2), true);
    FullereneBatch batch = FullereneBatch(N, std::min(1,n_max_isomers));
    FullereneQueue queue = FullereneQueue(N, std::min(1,n_max_isomers));

    BuckyGen::buckygen_queue BQ = BuckyGen::start(N, false, false);
    SyclQueue Q = SyclQueue(Device::get_devices(DeviceType::GPU)[0]);
    void SetUp() override {
        BuckyGen::next_fullerene(BQ, G);
    }

    void TearDown() override {
        //BuckyGen::stop(BQ);
    }
};

class FullereneTest : public PushBackTest {
    void SetUp() override {
        PushBackTest::SetUp();
        batch.push_back(G, 0);
        queue.push_back(G, 0);
    }
};

class FullereneBatchTest : public PushBackTest {};

class FullereneQueueTest : public PushBackTest {};

class FullereneBatchViewTest : public PushBackTest {

    void SetUp() override {
        BuckyGen::buckygen_queue BQ;
        BQ = BuckyGen::start(N, false, false);
        auto n_resize = std::min(40, n_max_isomers);
        for (int i = 0; i < n_resize; i++) {
            BuckyGen::next_fullerene(BQ, G);
            batch.push_back(G, i);
        }
    }

    void TearDown() override {
        BuckyGen::stop(BQ);
    }
};

class FullereneFixture : public ::testing::Test {
protected:
    using Fullerene = Fullerene<float, uint16_t>;
    using FullereneBatch = FullereneBatch<float, uint16_t>;
    using FullereneQueue = FullereneQueue<float, uint16_t>;
    using FullereneBatchView = FullereneBatchView<float, uint16_t>;
};

TEST_P(PushBackTest, PushGraphFromBuckygenToFullereneBatch) {
    batch.push_back(G, 0);
    Graph G_batch(batch[0]);
    // Assert that the data is correctly transferred from Buckygen to batch and back
    EXPECT_EQ(G_batch.neighbours, G.neighbours);
}

TEST_P(PushBackTest, PushGraphFromBuckygenToFullereneQueue) {
    queue.push_back(G, 0);
    Graph G_queue(queue.front());
    // Assert that the data is correctly transferred from Buckygen to queue and back
    EXPECT_EQ(G_queue.neighbours, G.neighbours);
}

TEST_P(PushBackTest, PushMoreGraphsThanOriginalCapacityOfBatch) {
    batch.push_back(G, 0);
    batch.push_back(G, 1);
    batch.push_back(G, 2);
    // Assert that the batch is full
    EXPECT_EQ(batch.size(), 3);
}

TEST_P(FullereneTest, CopyAssignment) {
    auto fullerene = batch[0];
    EXPECT_EQ(fullerene.N_, N);
    EXPECT_EQ(fullerene.Nf_, N/2 + 2);
    auto view1 = FullereneBatchView(batch, 0, 1);
    EXPECT_EQ(fullerene.d_, view1.d_);
    EXPECT_EQ(fullerene.m_.flags_, view1.m_.flags_[0]);
    EXPECT_EQ(fullerene.m_.iterations_, view1.m_.iterations_[0]);
    EXPECT_EQ(fullerene.m_.ID_, view1.m_.ID_[0]);
    EXPECT_EQ(fullerene.m_.iterations_, view1.m_.iterations_[0]);
}

TEST_P(FullereneTest, EqualityOperator) {
    auto fullerene1 = batch[0];
    auto fullerene2 = queue[0];
    EXPECT_TRUE(fullerene1 == fullerene2);
}

TEST_P(FullereneTest, CopyMethod) {
    auto fullerene = batch[0];
    auto fullerene_copy = fullerene;
    EXPECT_EQ(fullerene_copy, fullerene);
}

TEST_P(FullereneTest, CopyConstructor) {
    Fullerene fullerene_copy(batch[0]);
    EXPECT_EQ(fullerene_copy, batch[0]);
}

TEST_P(FullereneTest, CopyAssignmentOperator) {
    Fullerene fullerene = batch[0];
    EXPECT_EQ(fullerene, batch[0]);
}

TEST_P(FullereneTest, CopyAssignmentOperatorFromGraph) {
    Fullerene fullerene = batch[0];
    EXPECT_EQ(fullerene, batch[0]);
    fullerene = G;
    EXPECT_EQ(fullerene, batch[0]);
}

TEST_P(FullereneTest, CopyAssignmentOperatorFromPolyhedron) {
    FullereneDual dual(G);
    dual.update();
    PlanarGraph PG = dual.dual_graph();
    PG.layout2d = PG.tutte_layout();
    Polyhedron P(PG);
    P.points = P.zero_order_geometry();
    P.optimize();
    FullereneBatch batch1(N, 1);
    batch1.push_back(P);
    EXPECT_FLOAT_EQ(((Polyhedron)batch1[0]).volume_divergence(), P.volume_divergence());
}

// Test FullereneBatch class
TEST_P(FullereneBatchTest, Constructor) {
    FullereneBatch test_batch(N, 1);
    EXPECT_EQ(test_batch.size(), 0);
    EXPECT_EQ(test_batch.capacity(), 1);
}

TEST_P(FullereneBatchTest, PushBack) {
    FullereneBatch batch(N, 5);
    batch.push_back(G);
    EXPECT_EQ(batch.size(), 1);
    EXPECT_EQ(batch.capacity(), 5);
}

TEST_P(FullereneBatchTest, Resize) {
    FullereneBatch batch(N, 5);
    batch.resize(10);
    EXPECT_EQ(batch.size(), 0);
    EXPECT_EQ(batch.capacity(), 10);
}

TEST_P(FullereneBatchTest, IndexOperator) {
    FullereneBatch batch(N, 5);
    batch.push_back(G);
    Fullerene fullerene = batch[0];
    Graph G_batch(fullerene);
    EXPECT_EQ(fullerene.N_, N);
    EXPECT_EQ(G_batch.neighbours, G.neighbours);
}

TEST_P(FullereneBatchTest, EqualityOperator) {
    FullereneBatch batch1(N, 1);
    FullereneBatch batch2(N, 1);
    batch1.push_back(G);
    batch2.push_back(G);
    auto batch1_view = (FullereneBatchView)batch1;
    auto batch2_view = (FullereneBatchView)batch2;

    EXPECT_EQ(batch1, batch2);
    EXPECT_EQ(batch1_view, batch2_view);
}

// Test FullereneQueue class
TEST_P(FullereneQueueTest, Constructor) {
    FullereneQueue queue;
    EXPECT_EQ(queue.size(), 0);
    EXPECT_EQ(queue.capacity(), 0);
}

TEST_P(FullereneQueueTest, AssignmentOperator) {
    FullereneQueue queue1(N, 1);
    queue1.push_back(G);
    FullereneQueue queue2 = queue1;
    EXPECT_EQ(queue1, queue2);
}

TEST_P(FullereneQueueTest, PushBack) {
    FullereneQueue queue(N, 5);
    queue.push_back(G);
    EXPECT_EQ(queue.size(), 1);
    EXPECT_EQ(queue.capacity(), 5);
}

TEST_P(FullereneQueueTest, PushQueueToBatch) {
    FullereneQueue queue2;
    FullereneBatch batch(N, 3);
    for (int i = 0; i < 5; i++) {
        queue2.push_back(G);
    }
    SyclQueue Q;
    QueueUtil::push(Q, batch, queue2, ConditionFunctor(0, StatusEnum::DUAL_INITIALIZED));
    EXPECT_EQ(queue2.size(), 2);
    EXPECT_EQ(batch.size(), 3);
    for (int i = 0; i < 3; i++) {
        EXPECT_EQ(Graph(batch[i]).neighbours, G.neighbours);
    }
}


TEST_P(FullereneQueueTest, QueueManipulation) {
    FullereneQueue queue;
    FullereneBatch batch(N, 3);
    EXPECT_EQ(queue.size(), 0);
    EXPECT_EQ(queue.capacity(), 0);
    auto full_ix = 0;
    for(int i = 0; i < 5; i++){
        queue.push_back(G, full_ix++); 
    } //Queue starts at 0 capacity and resizes to 1, 2, 4, 8
    EXPECT_EQ(queue.capacity(), 8);

    for (int i = 0; i < 5; i++) {
        EXPECT_EQ(queue[i].m_.ID_.get(), i);    
    }

    EXPECT_EQ(queue.size(), 5);
    EXPECT_EQ(queue.front_index(), 0);
    EXPECT_EQ(queue.back_index(), 4);

    SyclQueue Q;
    QueueUtil::push(Q, queue, batch, ConditionFunctor(StatusEnum::DUAL_INITIALIZED));
    EXPECT_EQ(queue.size(), 5);
    EXPECT_EQ(queue.front_index(), 0);
    EXPECT_EQ(queue.back_index(), 4);

    QueueUtil::push(Q, batch, queue, ConditionFunctor(0, StatusEnum::DUAL_INITIALIZED));
    EXPECT_EQ(queue.size(), 2);
    EXPECT_EQ(queue.front_index(), 3);
    EXPECT_EQ(queue.back_index(), 4);

    for (int i = 0; i < 5; i++) {
        queue.push_back(G, full_ix++);
    }

    EXPECT_EQ(queue.size(), 7);
    EXPECT_EQ(queue.front_index(), 3);
    EXPECT_EQ(queue.back_index(), 1);
    std::array expected_batch_IDs = {0, 1, 2};
    std::array expected_queue_IDs = {3, 4, 5, 6, 7, 8, 9};


    for (int i = 0; i < 3; i++) {
        EXPECT_EQ(batch[i].m_.ID_.get(), expected_batch_IDs[i]);
    }

    for (int i = 0; i < 7; i++) {
        EXPECT_EQ(queue[i].m_.ID_.get(), expected_queue_IDs[i]);
    }

    queue.push_back(G, full_ix++);

    EXPECT_EQ(queue.size(), 8); //Queue is full
    EXPECT_EQ(queue.capacity(), 8);
    EXPECT_EQ(queue.front_index(), 3);
    EXPECT_EQ(queue.back_index(), 2);

    queue.push_back(G, full_ix++);

    EXPECT_EQ(queue.size(), 9); //Queue will now have resized to 16
    EXPECT_EQ(queue.capacity(), 16);
    EXPECT_EQ(queue.front_index(), 0);
    EXPECT_EQ(queue.back_index(), 8);

    std::array expected_queue_IDs_2 = {3, 4, 5, 6, 7, 8, 9, 10, 11};

    for (int i = 0; i < 9; i++) {
        EXPECT_EQ(queue[i].m_.ID_.get(), expected_queue_IDs_2[i]);
    }

    QueueUtil::push(Q, queue, batch, ConditionFunctor(StatusEnum::DUAL_INITIALIZED));

    EXPECT_EQ(queue.size(), 12);
    EXPECT_EQ(queue.front_index(), 0);
    EXPECT_EQ(queue.back_index(), 11);

    std::array expected_queue_IDs_3 = {3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2};
    
    for (int i = 0; i < 12; i++) {
        EXPECT_EQ(queue[i].m_.ID_.get(), expected_queue_IDs_3[i]);
        EXPECT_EQ(Graph(queue[i]).neighbours, G.neighbours);
    }

}

// Test That None of the non-owning containers have default constructors (This is to prevent accidental misuse)

TEST_F(FullereneFixture, NoDefaultViewConstructor) {
    EXPECT_FALSE(std::is_default_constructible_v<FullereneBatchView>);
}

TEST_F(FullereneFixture, NoDefaultFullereneConstructor) {
    EXPECT_FALSE(std::is_default_constructible_v<Fullerene>);
}

TEST_P(FullereneTest, DeletingFullereneDoesNotInvalidateBatch) {
    EXPECT_EQ(batch.size(), 1);
    auto ptr = batch.d_.A_cubic_.data();
    {
        Fullerene fullerene = batch[0];
        fullerene.m_.flags_ = StatusEnum::DUAL_INITIALIZED;
    }
    EXPECT_EQ(batch.size(), 1);
    EXPECT_EQ(batch.d_.A_cubic_.data(), ptr);
}

/* TEST_P(FullereneTest, AccessingFullereneAfterDeletingBatchShouldSegfault) {
    auto return_fullerene = [&](){
        FullereneBatch batch(N, 1);
        batch.push_back(G);
        return Fullerene(batch[0]);
    };
    auto fullerene = return_fullerene();
    EXPECT_TRUE(fullerene.N_ == N);
    ASSERT_EXIT((std::cout << fullerene.d_.A_cubic_[0], exit(0)), ::testing::KilledBySignal(SIGSEGV), ".*");
} */



TEST_P(FullereneBatchViewTest, ConstructFromBatch) {
    FullereneBatchView view(batch, 0, 1);
    EXPECT_EQ(view.size(), 1);
    EXPECT_EQ(view.capacity(), 1);
    EXPECT_TRUE(std::equal(view.begin(), view.end(), batch.begin())   );
}

TEST_P(FullereneBatchViewTest, ConstructFromBatchWithRange) {
    FullereneBatchView view(batch, 0, 5);
    EXPECT_EQ(view.size(), 5);
    EXPECT_EQ(view.capacity(), 5);
    EXPECT_TRUE(std::equal(view.begin(), view.end(), batch.begin()));
}

TEST_P(FullereneBatchViewTest, ConstructFromView) {

    FullereneBatchView view1(batch, 0, 5);
    FullereneBatchView view2(view1, 2, 3);
    EXPECT_EQ(view2.size(), 3);
    EXPECT_EQ(view2.capacity(), 3);
    EXPECT_TRUE(std::equal(view2.begin(), view2.end(), view1.begin() + 2));
    EXPECT_TRUE(std::equal(view2.begin(), view2.end(), batch.begin() + 2));
}

TEST_P(FullereneBatchViewTest, ConstructFromIteratorPair) {
    FullereneBatchView view(batch.begin() + 2, batch.end());
    FullereneBatchView view2(batch, 2);
    EXPECT_EQ(view.size(), batch.size() - 2);
    EXPECT_EQ(view.capacity(), batch.size() - 2);
    EXPECT_EQ(view, view2);
    EXPECT_TRUE(std::equal(view.begin(), view.end(), batch.begin() + 2));
}


TEST_P(FullereneBatchViewTest, EqualityOperator) {
    FullereneBatchView view1(batch, 3, 7);
    FullereneBatchView view2(batch, 3, 7);
    EXPECT_EQ(view1, view2);
}

TEST_P(FullereneBatchViewTest, CopyAssignment) {
    FullereneBatchView view1(batch, 0, 1);
    FullereneBatchView view2 = view1;
    EXPECT_EQ(view1, view2);
}

TEST_P(FullereneBatchViewTest, CopyConstruction) {
    FullereneBatchView view1(batch, 0, 1);
    FullereneBatchView view2(view1);
    EXPECT_EQ(view1, view2);
}

TEST_P(FullereneBatchViewTest, ViewIndexOperator) {
    FullereneBatchView view(batch, 5, 10);
    Fullerene fullerene = view[0];
    Fullerene fullerene2 = batch[5];
    EXPECT_EQ(fullerene, fullerene2);
}

TEST_P(FullereneQueueTest, Resize) {
    FullereneQueue queue(N, 5);
    queue.resize(10);
    EXPECT_EQ(queue.size(), 0);
    EXPECT_EQ(queue.capacity(), 10);
}

TEST_P(FullereneQueueTest, IndexOperator) {
    FullereneQueue queue(N, 5);
    queue.push_back(G);
    Fullerene fullerene = queue[0];
    Graph G_queue(fullerene);
    EXPECT_EQ(fullerene.N_, N);
    EXPECT_EQ(G_queue.neighbours, G.neighbours);
}

TEST_P(FullereneQueueTest, EqualityOperator) {
    FullereneQueue queue1(N, 1);
    FullereneQueue queue2(N, 1);
    queue1.push_back(G);
    queue2.push_back(G);
    EXPECT_EQ(queue1, queue2);
}

    
INSTANTIATE_TEST_SUITE_P(, FullereneQueueTest, ::testing::Range(60, 61, 20));
INSTANTIATE_TEST_SUITE_P(, FullereneBatchTest, ::testing::Range(60, 61, 20));
INSTANTIATE_TEST_SUITE_P(, PushBackTest, ::testing::Range(60, 61, 20));
INSTANTIATE_TEST_SUITE_P(, FullereneTest, ::testing::Range(60, 61, 20));
INSTANTIATE_TEST_SUITE_P(, FullereneBatchViewTest, ::testing::Range(60, 61, 20));

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}