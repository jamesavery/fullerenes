#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "forcefield2.cc" // Include CC-file
using namespace std;

namespace python_api { 
  namespace py = pybind11;
  template <typename T> using
  np_array = py::array_t<T, py::array::c_style | py::array::forcecast>; 

  template <typename T> T *ptr(np_array<T> &np_A)
  {
    auto A_info = np_A.request();
    return static_cast<T*>(A_info.ptr);
  }
  
  FullereneForcefield forcefield(np_array<int>&    neighbours,
				 np_array<real_t>& X0,
				 np_array<uint8_t>& face_right,
				 np_array<int>& next_on_face,
				 np_array<int>& prev_on_face)
  {
    auto X_info = X0.request();
    assert(X_info.shape.size() == 2 && X_info.shape[1] == 3);

    size_t   N     = X_info.shape[0];
    coord3d *X_ptr = static_cast<coord3d*>(X_info.ptr);
    
    return FullereneForcefield(N,
			       ptr(neighbours), X_ptr,
			       ptr(face_right), ptr(next_on_face),
			       ptr(prev_on_face));
  }
  
}


PYBIND11_MODULE(forcefield, m) {
    m.doc() = "Example Python Module in C++"; // optional module docstring

}
