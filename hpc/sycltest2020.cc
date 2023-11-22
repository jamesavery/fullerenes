#include <iostream>

#include <CL/sycl.hpp>

using namespace cl::sycl;

int main(int, char**) {
  float4 a = { 1.0, 2.0, 3.0, 4.0 };
  float4 b = { 4.0, 3.0, 2.0, 1.0 };
  float4 c = { 0.0, 0.0, 0.0, 0.0 };

  auto selector =  gpu_selector_v;
  queue queue(selector, property::queue::in_order{});
  
  std::cout << "Running on "
	    << queue.get_device().get_info<info::device::name>()
	    << "\n";

  {
    { // start of scope, ensures data copied back to host
      buffer<float4, 1> a_sycl(&a, range<1>(1));
      buffer<float4, 1> b_sycl(&b, range<1>(1));
      buffer<float4, 1> c_sycl(&c, range<1>(1));

      queue.submit([&] (handler& h) {
	auto a_acc = a_sycl.get_access<access::mode::read>(h);
	auto b_acc = b_sycl.get_access<access::mode::read>(h);
	auto c_acc = c_sycl.get_access<access::mode::discard_write>(h);

							
	h.single_task<class vector_addition>([=] () {
	  c_acc[0] = a_acc[0] + b_acc[0];
	});

	//	h.parallel_for(1024, [=](id<1> idx) { });
      });
    } // end of scope, ensures data copied back to host
    

  }

  std::cout << "  A { " << a.x() << ", " << a.y() << ", " << a.z() << ", " << a.w() << " }\n"
	    << "+ B { " << b.x() << ", " << b.y() << ", " << b.z() << ", " << b.w() << " }\n"
	    << "------------------\n"
	    << "= C { " << c.x() << ", " << c.y() << ", " << c.z() << ", " << c.w() << " }"
	    << std::endl;


  return 0;
}

