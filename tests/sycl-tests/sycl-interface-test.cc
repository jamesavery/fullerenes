#include <iostream>
#include <numeric>
#include <algorithm>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/sycl-headers/sycl-wrappers.hh>
#include <fullerenes/sycl-headers/sycl-kernels.hh>
#include <fullerenes/isomerdb.hh>
#include <sys/types.h>


int main(int argc, char** argv) {
    constexpr int N = 30;
    BuckyGen::buckygen_queue BQ = BuckyGen::start(N, false, false);
    Graph G(neighbours_t(N/2 + 2));

    auto n_isomers = IsomerDB::number_isomers(N);
    FullereneBatch<float, uint16_t> batch(N, n_isomers);
    FullereneQueue<float, uint16_t> FQ(N, n_isomers);
    
    BuckyGen::next_fullerene(BQ, G);

    //std::cout << "Verify Equality of first isomer: " << std::boolalpha << (batch[0] == FQ[0]) << std::endl;

    
    std::cout << "Number of isomers: " << n_isomers << std::endl;
    batch.resize(n_isomers);
    std::for_each(batch.begin(), batch.end(), 
                [&, i = 0] (auto fullerene) mutable {
                BuckyGen::next_fullerene(BQ, G);
                batch.push_back(G, i++); });
   
    
    ConditionFunctor c1(StatusFlag::CONVERGED_2D, StatusFlag::EMPTY | StatusFlag::CONVERGED_3D);


    Device device = Device::get_devices(DeviceType::GPU)[0];
    std::cout << device.get_name() << std::endl;
    SyclQueue Q(device, true);
    
    DualizeFunctor<float, uint16_t> dualize;
    auto batchview = FullereneBatchView(batch);
    
    dualize(Q, batch[1], LaunchPolicy::ASYNC);
    //std::for_each_n(batchview.begin() + 1, 2, [&](auto fullerene){dualize(Q, fullerene, LaunchPolicy::ASYNC);} );

    push(Q, FQ, batch, ConditionFunctor(StatusFlag::CUBIC_INITIALIZED));
    push(Q, batch, FQ, ConditionFunctor(0, StatusFlag::CUBIC_INITIALIZED));
    //std::cout << "Device has " << device.get_property(DeviceProperty::MAX_COMPUTE_UNITS) << " compute units\n";
    //std::cout << "Device has " << device.get_property(DeviceProperty::MAX_WORK_GROUP_SIZE) << " max work group size\n";
    //std::cout << "Device has " << device.get_property(DeviceProperty::MAX_CLOCK_FREQUENCY) << " max clock frequency\n";

    
    return 0;
}