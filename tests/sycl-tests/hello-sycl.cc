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
  
  // queue queue(selector, property::queue::in_order{});
  
  // std::cout << "Running on "
  // 	    << queue.get_device().get_info<info::device::name>()
  // 	    << "\n";

  // {
  //   { // start of scope, ensures data copied back to host
  //     buffer<float4, 1> a_sycl(&a, range<1>(1));
  //     buffer<float4, 1> b_sycl(&b, range<1>(1));
  //     buffer<float4, 1> c_sycl(&c, range<1>(1));

  //     queue.submit([&] (handler& h) {
  // 	auto a_acc = a_sycl.get_access<access::mode::read>(h);
  // 	auto b_acc = b_sycl.get_access<access::mode::read>(h);
  // 	auto c_acc = c_sycl.get_access<access::mode::discard_write>(h);

							
  // 	h.single_task<class vector_addition>([=] () {
  // 	  c_acc[0] = a_acc[0] + b_acc[0];
  // 	});

  // 	//	h.parallel_for(1024, [=](id<1> idx) { });
  //     });
  //   } // end of scope, ensures data copied back to host
    

  // }

  // std::cout << "  A { " << a.x() << ", " << a.y() << ", " << a.z() << ", " << a.w() << " }\n"
  // 	    << "+ B { " << b.x() << ", " << b.y() << ", " << b.z() << ", " << b.w() << " }\n"
  // 	    << "------------------\n"
  // 	    << "= C { " << c.x() << ", " << c.y() << ", " << c.z() << ", " << c.w() << " }"
  // 	    << std::endl;


  return 0;
}

