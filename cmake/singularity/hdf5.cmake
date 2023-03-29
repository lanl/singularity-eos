macro(singularity_enable_hdf5 target)
  # CMake's FindHDF5 package ignores any path in CMAKE_PREFIX_PATH that
  # doesn't contain the hdf5-config.cmake file in a certain location.
  # Some of our HDF5 builds don't have this.  So we have to give CMake
  # a little extra help...
  if(NOT HDF5_ROOT)
    find_path(HDF5_INCLUDE_DIR hdf5.h
      HINTS ENV HDF5_ROOT
      PATH_SUFFIXES include)
    get_filename_component(HDF5_ROOT ${HDF5_INCLUDE_DIR} DIRECTORY)
  endif()

  find_package(HDF5 COMPONENTS HL REQUIRED)


  target_include_directories(${target} SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(${target} PUBLIC ${HDF5_LIBRARIES})

  if(HDF5_IS_PARALLEL)
    find_package(MPI COMPONENTS CXX REQUIRED)
    target_link_libraries(${target} PUBLIC MPI::MPI_CXX)
  endif()

  target_compile_definitions(${target} PUBLIC SPINER_USE_HDF5)
  target_compile_definitions(${target} PUBLIC SINGULARITY_USE_HDF5)

endmacro()

