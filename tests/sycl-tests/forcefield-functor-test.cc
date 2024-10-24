#include <fullerenes/sycl-headers/all-kernels.hh>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/isomerdb.hh>
#include <numeric>
#include <execution>
#include <gtest/gtest.h>
#include <iomanip>


class ForceFieldTest : public ::testing::TestWithParam<int> {
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

    DualizeFunctor<float, uint16_t> dualize;
    TutteFunctor<float, uint16_t> tutte;
    SphericalProjectionFunctor<float, uint16_t> spherical_projection;
    ForcefieldOptimizeFunctor<ForcefieldType::PEDERSEN, float, uint16_t> FF;

    BuckyGen::buckygen_queue BQ = BuckyGen::start(N, false, false);
    SyclQueue Q = SyclQueue(Device::get_devices(DeviceType::GPU)[0]);
    void SetUp() override {
        BuckyGen::next_fullerene(BQ, G);
    }

    void TearDown() override {
        //BuckyGen::stop(BQ);
    }
};

TEST_P(ForceFieldTest, TestForceFieldOptimizeFunctor) {
    auto device = Device::get_devices(DeviceType::GPU)[0];
    for (int i = 0; i < device.get_property(DeviceProperty::MAX_COMPUTE_UNITS)*5; i++) {
        BuckyGen::next_fullerene(BQ, G);
        queue.push_back(G);
    }
    dualize(Q, queue, LaunchPolicy::SYNC);
    tutte(Q, queue, LaunchPolicy::SYNC);
    spherical_projection(Q, queue, LaunchPolicy::SYNC); //Isomers are now prepared for forcefield optimization

    FullereneBatch opt_batch(N, device.get_property(DeviceProperty::MAX_COMPUTE_UNITS));
    FullereneQueue out_queue(N, 1);

    std::vector<FullereneQueue> input_queues(10, FullereneQueue(N, device.get_property(DeviceProperty::MAX_COMPUTE_UNITS)*5));
    for (int i = 0; i < input_queues[0].capacity(); i++) {
        BuckyGen::next_fullerene(BQ, G);
        input_queues[0].push_back(G);
    }
    dualize(Q, input_queues[0], LaunchPolicy::SYNC);
    tutte(Q, input_queues[0], LaunchPolicy::SYNC);
    spherical_projection(Q, input_queues[0], LaunchPolicy::SYNC);
    
    for (int i = 1 ; i < 10; i++) {
        input_queues[i] = input_queues[0];
    }
    

    std::vector<FullereneBatch> opt_batches(10, FullereneBatch(N, device.get_property(DeviceProperty::MAX_COMPUTE_UNITS)));
    std::vector<FullereneQueue> out_queues(10, FullereneQueue(N, 1));
    
    for (int j = 0; j < 10; j++)    ASSERT_EQ(input_queues[j],input_queues[0]); //Sanity check
    
    for (int i = 0; i < 50; i++) {
        for (int j = 0; j < 10; j++) {
            QueueUtil::push(Q, opt_batches[j], input_queues[j], ConditionFunctor(0, 0, StatusEnum::EMPTY | StatusEnum::CONVERGED_3D | StatusEnum::FAILED_3D));
            FF(Q, opt_batches[j], LaunchPolicy::SYNC, N, 10*N);
            QueueUtil::push(Q, out_queues[j], opt_batches[j], ConditionFunctor(0, 0, StatusEnum::CONVERGED_3D | StatusEnum::FAILED_3D), StatusEnum::EMPTY);
        }
        bool all_equal = true;
        for (int j = 0; j < 10; j++) {
            all_equal &= out_queues[j] == out_queues[0];
        }
        ASSERT_TRUE(all_equal);
    } 

    std::cout << out_queues[0].m_.flags_ << std::endl;
    std::cout << out_queues[0].m_.iterations_ << std::endl;
    std::cout << out_queues[0].size() << std::endl;
    std::cout << input_queues[0].size() << std::endl;
    std::cout << std::transform_reduce(opt_batches[0].m_.flags_.begin(), opt_batches[0].m_.flags_.end(), 0, std::plus<int>(), [](auto x) {return x.is_set(StatusEnum::NOT_CONVERGED) ? 1 : 0;}) << std::endl;

    for (int i = 0; i < 21; i++) {
        /* QueueUtil::push(Q, opt_batch, queue, ConditionFunctor(0, 0, StatusEnum::EMPTY | StatusEnum::CONVERGED_3D | StatusEnum::FAILED_3D));
        FF(Q, opt_batch, LaunchPolicy::SYNC, N, N*20);
        std::cout << "Flags: " << opt_batch.m_.flags_ << std::endl;
        std::cout << "Iterations: " << opt_batch.m_.iterations_ << std::endl;

        auto qsize_before = out_queue.size();
        QueueUtil::push(Q, out_queue, opt_batch, ConditionFunctor(0, 0, StatusEnum::CONVERGED_3D | StatusEnum::FAILED_3D), StatusEnum::EMPTY);
        std::cout << "Pushed " << out_queue.size() - qsize_before << " to out_queue" << std::endl; */
    }


}

INSTANTIATE_TEST_SUITE_P(, ForceFieldTest, ::testing::Range(60, 61, 20));

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}