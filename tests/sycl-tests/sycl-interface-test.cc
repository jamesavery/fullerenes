#include <iostream>
#include <numeric>
#include <algorithm>
#include <fullerenes/sycl-wrappers.hh>
#include <fullerenes/buckygen-wrapper.hh>
#include <fullerenes/sycl-kernels.hh>

#include <unistd.h>
#include <sys/types.h>

int main(int argc, char** argv) {
    constexpr int N = 20;
    BuckyGen::buckygen_queue BQ = BuckyGen::start(N, false, false);
    Graph G(neighbours_t(N/2 + 2));

    BuckyGen::next_fullerene(BQ, G);

    FullereneBatch<float, uint16_t> batch(N, 1);
    batch.push_back(G);
    ConditionFunctor c1(StatusFlag::CONVERGED_2D, StatusFlag::EMPTY | StatusFlag::CONVERGED_3D);
    

    
    std::cout << "Batch Dual Graph: " << batch.d_.A_dual_ << std::endl;

    
    


    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <device_type>" << std::endl;
        return 1;
    }
    DeviceType type;
    if (std::string(argv[1]) == "cpu") { type = DeviceType::CPU; 
        std::cout << "Device type is CPU" << std::endl;
    } else if (std::string(argv[1]) == "gpu") { type = DeviceType::GPU;
        std::cout << "Device type is GPU" << std::endl;
    } else if (std::string(argv[1]) == "accelerator") { type = DeviceType::ACCELERATOR;
        std::cout << "Device type is Accelerator" << std::endl;
    } else {
        std::cerr << "Invalid device type" << std::endl;
        return 1;
    }

    SyclContext ctx = SyclContext();
    Device device(0, type);
    std::cout << ctx.device_get_name(device) << std::endl;
    SyclQueue Q(device, true);
    auto fullerene = batch[0];
    

    dualize(Q, fullerene, LaunchPolicy::SYNC);
    

    std::cout << "Batch Cubic Graph: " << batch.d_.A_cubic_ << std::endl;

    /* std::iota(vec.begin(), vec.end(), 0);
    std::transform(vec.begin(), vec.end(), vec.begin(), [](int x){return x*x;});
    for(auto val : vec){
        std::cout << val << std::endl;
    } */

    std::cout << "Device has " << ctx.device_get_property(device, DeviceProperty::MAX_COMPUTE_UNITS) << " compute units\n";
    std::cout << "Device has " << ctx.device_get_property(device, DeviceProperty::MAX_WORK_GROUP_SIZE) << " max work group size\n";
    std::cout << "Device has " << ctx.device_get_property(device, DeviceProperty::MAX_CLOCK_FREQUENCY) << " max clock frequency\n";

    
    return 0;
}