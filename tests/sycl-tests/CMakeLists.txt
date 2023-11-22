set(SYCL_EXECUTABLES
  hello-sycl
  sycl-pipeline-test
)

foreach(EXECUTABLE ${SYCL_EXECUTABLES})
  set_source_files_properties(${EXECUTABLE}.cc PROPERTIES COMPILE_FLAGS "${ICPX_FLAGS}")    
endforeach()
  
foreach(EXECUTABLE ${SYCL_EXECUTABLES})
  add_executable(${EXECUTABLE} ${EXECUTABLE}.cc)
endforeach()

foreach(EXECUTABLE ${SYCL_EXECUTABLES})
  target_include_directories(${EXECUTABLE} PRIVATE ${DPCPP_INCLUDE} ${DPCPP_INCLUDE}/sycl)
  target_link_libraries(${EXECUTABLE} PRIVATE sycl_fullerene_lib sycl fullerenes)
  target_link_options(${EXECUTABLE} PRIVATE 
    -O4
    -fsycl
    "SHELL:${DPCPP_TARGET_FLAGS}"
    #    "SHELL:-Xcuda-ptxas --verbose"
    #    "SHELL:-Xcuda-ptxas --maxrregcount=80"
    #    "SHELL:-Xcuda-ptxas --opt-level=3"
    #    "SHELL:-Xcuda-ptxas --generate-line-info"
    #    "SHELL:-Xcuda-ptxas --warn-on-double-precision-use"
    #    "SHELL:-Xcuda-ptxas --warn-on-spills"
    #    "SHELL:-Xcuda-ptxas --generate-line-info"
  )
endforeach()
