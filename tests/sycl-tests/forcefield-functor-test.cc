#include <fullerenes/sycl-headers/all-kernels.hh>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/isomerdb.hh>
#include <numeric>
#include <gtest/gtest.h>


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

    for (int i = 0; i < 21; i++) {
        QueueUtil::push(Q, opt_batch, queue, ConditionFunctor(0, 0, StatusEnum::EMPTY | StatusEnum::CONVERGED_3D | StatusEnum::FAILED_3D));
        FF(Q, opt_batch, LaunchPolicy::SYNC, N, N*20);
        std::cout << "Flags: " << opt_batch.m_.flags_ << std::endl;
        std::cout << "Iterations: " << opt_batch.m_.iterations_ << std::endl;

        auto Qsize_before = out_queue.size();
        QueueUtil::push(Q, out_queue, opt_batch, ConditionFunctor(0, 0, StatusEnum::CONVERGED_3D | StatusEnum::FAILED_3D), StatusEnum::EMPTY);
        std::cout << "Pushed " << out_queue.size() - Qsize_before << " to out_queue" << std::endl;
    }


}

INSTANTIATE_TEST_SUITE_P(, ForceFieldTest, ::testing::Range(60, 61, 20));

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}