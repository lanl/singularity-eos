macro(singularity_enable_hdf5 target)
  if(NOT HDF5_ROOT)
    find_path(HDF5_INCLUDE_DIR hdf5.h
      HINTS ENV HDF5_ROOT
      PATH_SUFFIXES include)
    get_filename_component(HDF5_ROOT ${HDF5_INCLUDE_DIR} DIRECTORY)
  endif()

  if(NOT HDF5_LIBRARIES)
    find_package(HDF5 COMPONENTS C HL REQUIRED)
  endif()

  target_include_directories(${target} SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(${target} PUBLIC ${HDF5_LIBRARIES})

  if(HDF5_IS_PARALLEL)
    enable_language(C)
    find_package(MPI COMPONENTS C CXX REQUIRED)
    target_link_libraries(${target} PUBLIC MPI::MPI_CXX)
  endif()

  target_compile_definitions(${target} PUBLIC SPINER_USE_HDF5)
  target_compile_definitions(${target} PUBLIC SINGULARITY_USE_HDF5)

endmacro()

