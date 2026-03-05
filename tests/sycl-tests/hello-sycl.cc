#include <CL/sycl.hpp>
#include <iostream>

using namespace cl::sycl;
using namespace std;

int main(int, char**) {

  // This should work, but doesn't: It seems LUMI-G has a separate platform for every GPU-device!
  /* platform cpu_platform(cpu_selector_v);
  platform gpu_platform(gpu_selector_v);

  cout << "CPU-devices:\n";
  for(auto d: cpu_platform.get_devices())
    cout << "\t" << d.get_info<info::device::name>() << "\n";

  cout << "GPU-devices:\n";
  for(auto d: gpu_platform.get_devices())
    cout << "\t" << d.get_info<info::device::name>() << "\n";

  // Try again in a more pedestrian way:
  cout << "GPU-devices with pedestrian traversal:\n";
      */
  vector<queue> gpu_queues;
  for(auto P: platform::get_platforms())
    for(auto d: P.get_devices())
      if(d.is_gpu()){
	      gpu_queues.push_back(queue(d));
      } 

  // All right, now try to execute a kernel -- this seems to provoke the break! With a more informative error: 
  // /opt/cray/pe/cce/15.0.1/cce-clang/x86_64/bin/llvm-link: /tmp/cooltmpdir-GiXTzq/graph.cc-cce-openmp-amdgcn-amd-amdhsa.amdgpu: error: Unknown attribute kind (86) (Producer: 'LLVM16.0.6' Reader: 'LLVM 15.0.6')
  // /opt/cray/pe/cce/15.0.1/cce-clang/x86_64/bin/llvm-link: error:  loading file '/tmp/cooltmpdir-GiXTzq/graph.cc-cce-openmp-amdgcn-amd-amdhsa.amdgpu'
  // SO: looks like we found the culprit.		      
  for(auto Q: gpu_queues){
    Q.submit([&](handler &h) {
      h.parallel_for(1024, [=](id<1> idx) {
        
      });      
    });

    Q.wait();

    Q.submit([&](handler &h) {
      h.parallel_for(nd_range(range{1024*10}, range{1024}), [=](nd_item<1> idx) {
        //Print group size
        if (idx.get_local_id(0) == 0) {
          if(idx.get_group(0) == 0)
            sycl::_V1::ext::oneapi::experimental::printf("Group size: %d\n", idx.get_local_range(0));
        }
      });      
    });

    //Throws compilation error: (Must specify both a global and local range for nd_range)
    /* Q.submit([&](handler &h) {
      h.parallel_for(nd_range(range{1024*10}), [=](nd_item<1> idx) {
        //Print group size
        if (idx.get_local_id(0) == 0) {
          sycl::_V1::ext::oneapi::experimental::printf("Group size: %d\n", idx.get_local_range(0));
        }
      });      
    }); */

    Q.wait();
    Q.submit([&](handler &h) {
      h.parallel_for_work_group(range{1024*10}, [=](group<1> g) {
        //Print group size
        if (g.get_local_id(0) == 0) {
          if(g.get_group_id(0) == 0)
            sycl::_V1::ext::oneapi::experimental::printf("Group size: %d\n", g.get_local_range(0));
        }
      });
    });
    Q.wait();
    }

    //Enumerate all platforms and devices
    int device_count = 0;
    int platform_count = 0;
    std::vector<sycl::platform> platforms = sycl::platform::get_platforms();
    //Does comparison work?
    for(int i = 0; i < platforms.size(); i++){
      for(int j = i; j < platforms.size(); j++){
        if(platforms[i] == platforms[j])
          cout << "Platform " << i << " is equal to platform " << j << "\n";
        else
          cout << "Platform " << i << " is not equal to platform " << j << "\n";
      }
    }

    std::vector<sycl::device> devices;
    for(auto P: platform::get_platforms()){
      for(auto d: P.get_devices()){
        devices.push_back(d);
      }
    }

    std::vector<sycl::backend> backends;
    for(auto P: platform::get_platforms()){
      for(auto d: P.get_devices()){
        backends.push_back(d.get_backend());
      }
    }

    //Does comparison work?
    for(int i = 0; i < backends.size(); i++){
      for(int j = i; j < backends.size(); j++){
        if(backends[i] == backends[j])
          cout << "Backend " << i << " is equal to backend " << j << "\n";
        else
          cout << "Backend " << i << " is not equal to backend " << j << "\n";
      }
    }

    //Does comparison work?
    for(int i = 0; i < devices.size(); i++){
      for(int j = i; j < devices.size(); j++){
        if(devices[i] == devices[j])
          cout << "Device " << i << " is equal to device " << j << "\n";
        else
          cout << "Device " << i << " is not equal to device " << j << "\n";
      }
    }

    //What's different about device 2 and 3?
    devices[2].get_info<info::device::vendor_id>() == devices[3].get_info<info::device::vendor_id>()
      ? cout << "Device 2 and 3 have the same vendor id\n"
      : cout << "Device 2 and 3 have different vendor ids\n";

    devices[2].get_info<info::device::name>() == devices[3].get_info<info::device::name>()
      ? cout << "Device 2 and 3 have the same name\n"
      : cout << "Device 2 and 3 have different names\n";

    devices[2].get_info<info::device::driver_version>() == devices[3].get_info<info::device::driver_version>()
      ? cout << "Device 2 and 3 have the same driver version\n"
      : cout << "Device 2 and 3 have different driver versions\n";

    devices[2].get_info<info::device::profile>() == devices[3].get_info<info::device::profile>()
      ? cout << "Device 2 and 3 have the same profile\n"
      : cout << "Device 2 and 3 have different profiles\n";

    devices[2].get_info<info::device::version>() == devices[3].get_info<info::device::version>()
      ? cout << "Device 2 and 3 have the same version\n"
      : cout << "Device 2 and 3 have different versions\n";

    for(auto P: platform::get_platforms()){
      cout << "Platform has id: " << platform_count++ << " " ;
      for(auto d: P.get_devices()){
        d.has(aspect::emulated)
          ? cout << "Emulated device: "
          : cout << "Physical device: ";
          if (d.is_gpu())
            cout << "GPU: ";
          else if (d.is_cpu())
            cout << "CPU: ";
          else if (d.is_accelerator())
            cout << "Accelerator: ";
          else
            cout << "Unknown: ";
          cout << "\t" << d.get_info<info::device::name>() << "\n";
          device_count++;

        
      }

    }





  
  
  
  return 0;
}

