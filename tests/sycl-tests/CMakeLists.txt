
set(SYCL_EXECUTABLES
  sycl-pipeline-test
)

include_directories(/opt/intel/oneapi/compiler/2023.2.1/linux/include/sycl/CL)

foreach(EXECUTABLE ${SYCL_EXECUTABLES})
  set_source_files_properties(${EXECUTABLE}.cc PROPERTIES COMPILE_FLAGS "${ICPX_FLAGS}")    
endforeach()
  
foreach(EXECUTABLE ${SYCL_EXECUTABLES})
  add_executable(${EXECUTABLE} ${EXECUTABLE}.cc)
endforeach()

foreach(EXECUTABLE ${SYCL_EXECUTABLES})
  target_include_directories(${EXECUTABLE} PRIVATE "/opt/intel/oneapi/compiler/2023.2.1/linux/include/sycl" "/opt/intel/oneapi/compiler/2023.2.1/linux/include")
  target_link_libraries(${EXECUTABLE} PRIVATE sycl_fullerene_lib sycl fullerenes)
  target_link_options(${EXECUTABLE} PRIVATE 
    -O4
    -fsycl
    "SHELL:${ICPX_TARGET_FLAG}"
    "SHELL:-Xcuda-ptxas --verbose"
    "SHELL:-Xcuda-ptxas --maxrregcount=80"
    "SHELL:-Xcuda-ptxas --opt-level=3"
    "SHELL:-Xcuda-ptxas --generate-line-info"
    "SHELL:-Xcuda-ptxas --warn-on-double-precision-use"
    "SHELL:-Xcuda-ptxas --warn-on-spills"
    "SHELL:-Xcuda-ptxas --generate-line-info"
  )
endforeach()
