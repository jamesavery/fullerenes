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
  vector<queue> gpu_queues;
  for(auto P: platform::get_platforms())
    for(auto d: P.get_devices())
      if(d.is_gpu()){
	cout << "\t" << d.get_info<info::device::name>() << "\n";
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
  }
  
  
  return 0;
}

