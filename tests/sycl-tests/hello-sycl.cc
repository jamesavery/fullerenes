#include <CL/sycl.hpp>
#include <iostream>

using namespace cl::sycl;
using namespace std;

int main(int, char**) {

  // This should work, but doesn't: It seems LUMI-G has a separate platform for every GPU-device!
  platform cpu_platform(cpu_selector_v);
  platform gpu_platform(gpu_selector_v);

  cout << "CPU-devices:\n";
  for(auto d: cpu_platform.get_devices())
    cout << "\t" << d.get_info<info::device::name>() << "\n";

  cout << "GPU-devices:\n";
  for(auto d: gpu_platform.get_devices())
    cout << "\t" << d.get_info<info::device::name>() << "\n";

  // Try again in a more pedestrian way:
  cout << "GPU-devices with pedestrian traversal:\n";  
  for(auto P: platform::get_platforms())
    for(auto d: P.get_devices())
      if(d.is_gpu()) cout << "\t" << d.get_info<info::device::name>() << "\n";

  return 0;
}

