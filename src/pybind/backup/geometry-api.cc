#include "geometry.cc"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>


namespace python_api { 
  namespace py = pybind11;

  typedef py::array_t<voxel_type, py::array::c_style | py::array::forcecast> np_voxelarray;
  typedef py::array_t<real_t, py::array::c_style | py::array::forcecast> np_realarray;  

  array<real_t,3> center_of_mass(const np_voxelarray &np_voxels){
    auto voxels_info    = np_voxels.request();

    return ::center_of_mass({voxels_info.ptr,voxels_info.shape});
  }


  array<real_t,9> inertia_matrix(const np_voxelarray &np_voxels, array<real_t,3>& cm){
    auto voxels_info    = np_voxels.request();
    
    return ::inertia_matrix({voxels_info.ptr,voxels_info.shape}, cm);
  }

  array<real_t,9> inertia_matrix_serial(const np_voxelarray &np_voxels, array<real_t,3>& cm){
    auto voxels_info    = np_voxels.request();
    
    return ::inertia_matrix_serial({voxels_info.ptr,voxels_info.shape}, cm);
  }  


  void sample_plane(const np_voxelarray &np_voxels,
		    const array<real_t,3> &cm,
		    const array<real_t,3> &u_axis,
		    const array<real_t,3> &v_axis,		    
		    np_voxelarray &np_plane_samples,
		    const array<real_t,3> &L)
  {
    auto voxels_info = np_voxels.request();
    auto plane_samples_info  = np_plane_samples.request();
    
    ::sample_plane({voxels_info.ptr, voxels_info.shape},
		   {cm,u_axis,v_axis},
		   {plane_samples_info.ptr, plane_samples_info.shape},
		   L);
  }


  void integrate_axes(const np_voxelarray &np_voxels,
		    const array<real_t,3> &x0,		    
		    const array<real_t,3> &v_axis,
		    const array<real_t,3> &w_axis,
		    const real_t v_min, const real_t w_min,
		    np_realarray &output)
  {
    auto voxels_info = np_voxels.request();
    auto output_info  = output.request();

    ::integrate_axes({voxels_info.ptr, voxels_info.shape},
		     x0,v_axis,w_axis,
		     v_min, w_min,
		     {output_info.ptr, output_info.shape});
  }

  void zero_outside_bbox(const array<real_t,9> &principal_axes,
			 const array<real_t,6> &parameter_ranges,
			 const array<real_t,3> &cm, // TOOD: Med eller uden voxelsize?
			 np_voxelarray &np_voxels)
  {
    auto voxels_info = np_voxels.request();
    
    ::zero_outside_bbox(principal_axes,
		      parameter_ranges,
		      cm, 
		      {voxels_info.ptr, voxels_info.shape});
  }
}



PYBIND11_MODULE(geometry, m) {
    m.doc() = "Voxel Geometry Module"; // optional module docstring

    m.def("center_of_mass",       &python_api::center_of_mass);
    m.def("inertia_matrix",       &python_api::inertia_matrix);
    m.def("inertia_matrix_serial",&python_api::inertia_matrix_serial);
    m.def("integrate_axes",       &python_api::integrate_axes);        
    m.def("sample_plane",         &python_api::sample_plane);
    m.def("zero_outside_bbox",    &python_api::zero_outside_bbox);    
}
