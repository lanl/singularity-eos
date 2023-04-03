macro(singularity_import_kokkos)
  
  set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
  if(SINGULARITY_USE_CUDA)
    set(Kokkos_ENABLE_CUDA ON CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE ON CACHE BOOL "" FORCE)
  endif()

  if(NOT TARGET Kokkos::kokkos)
    add_subdirectory(utils/kokkos)
  endif()
endmacro()

macro(singularity_find_kokkos)
  find_package(Kokkos REQUIRED)
endmacro()

macro(singularity_enable_kokkos target)

  target_link_libraries(${target} PUBLIC Kokkos::kokkos)

  target_compile_definitions(${target} PUBLIC PORTABILITY_STRATEGY_KOKKOS)
  #TODO: shouldn't be needed
  target_compile_definitions(${target} PUBLIC SPINER_USE_KOKKOS)

  #  if(SINGULARITY_USE_KOKKOSKERNELS AND SINGULARITY_BUILD_CLOSURE)
  #  target_link_libraries(${target} PUBLIC Kokkos::kokkoskernels)
  #endif()
endmacro()

macro(singularity_import_kokkoskernels)
  set(KokkosKernels_ENABLE_TPL_BLAS OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_ENABLE_TPL_MKL OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_ENABLE_TPL_LAPACK OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_ENABLE_TPL_CUBLAS OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_ENABLE_TPL_CUSPARSE OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_ENABLE_TPL_MAGMA OFF CACHE BOOL "" FORCE)
  # Disable ETIs
  set(KokkosKernels_INST_COMPLEX_DOUBLE OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_COMPLEX_FLOAT OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_DOUBLE OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_EXECSPACE_OPENMP OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_EXECSPACE_SERIAL OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_EXECSPACE_THREADS OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_EXECSPACE_CUDA OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_FLOAT OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_LAYOUTLEFT OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_LAYOUTRIGHT OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_MEMSPACE_HOSTSPACE OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_OFFSET_INT OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_OFFSET_SIZE_T OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_ORDINAL_INT OFF CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_ORDINAL_INT OFF CACHE BOOL "" FORCE)
  if(NOT TARGET Kokkos::kokkoskernels)
    add_subdirectory(utils/kokkos-kernels)
  endif()
endmacro()

macro(singularity_find_kokkoskernels)
  find_package(KokkosKernels REQUIRED)
endmacro()

macro(singularity_enable_kokkoskernels target)
  target_compile_definitions(${target} PUBLIC SINGULARITY_USE_KOKKOSKERNELS)
  target_link_libraries(${target} PUBLIC Kokkos::kokkoskernels)
endmacro()

macro(_singularity_handle_custom_kokkos_tpl)
  set(options)
  set(one_value_args
    NAME
  )
  set(multi_value_args)

  cmake_parse_arguments(KTPL "${options}" "${one_value_args}"
    "${multi_value_args}" ${ARGN})

  if(NOT ("${KTPL_NAME}" STREQUAL "Kokkos" OR "${KTPL_NAME}" STREQUAL "KokkosKernels") )
    message(FATAL_ERROR "Custom TPL only supported for \"Kokkos\" and \"KokkosKernels\" (got \"${KTPL_NAME}\")")
  endif()

  string(TOUPPER ${KTPL_NAME} tplUp)
  string(TOLOWER ${KTPL_NAME} tplLo)
    
  # requires source dir
  if(NOT SINGULARITY_CUSTOM_${tplUp}_DIR)
    message(FATAL_ERROR "Selected custom ${KTPL_NAME} source, but SINGULARITY_CUSTOM_${tplUp}_DIR not set."
     "Set and reconfigure, or disable custom source.")
  else()
    if(NOT EXISTS ${SINGULARITY_CUSTOM_${tplUp}_DIR})
      message(FATAL_ERROR "Custom ${KTPL_NAME} source dir ${SINGULARITY_CUSTOM_${tplUP}_DIR} does not exist.")
    endif()
  endif()
  if(NOT SINGULARITY_CUSTOM_${tplUp}_CACHE)
    message(WARNING "Selected custom ${KTPL_NAME} source, but SINGULARITY_CUSTOM_${tplUp}_CACHE not set."
      " CMake will configure with default options."
      " See ${PROJECT_SOURCE_DIR}/utils/scripts/caches for sample caches.")
    set(SINGULARITY_CUSTOM_${tplUp}_CACHE "${PROJECT_SOURCE_DIR}/utils/scripts/caches/${tplLo}_default.cmake")
  endif() 
  if(NOT EXISTS ${SINGULARITY_CUSTOM_${tplUp}_CACHE})
    message(FATAL "Custom kokkos CMake cache file ${SINGULARITY_CUSTOM_${tplUp}_CACHE} does not exist.")
  endif()
  include(${SINGULARITY_CUSTOM_${tplUp}_CACHE})

  if(NOT TARGET Kokkos::${tplLo})
    add_subdirectory(${SINGULARITY_CUSTOM_${tplUp}_DIR} ${PROJECT_BINARY_DIR}/_deps/${tplLo}-build)
  else()
    message(WARNING "${KTPL_NAME} target already present, not using custom source directory.")
  endif()
 
endmacro()
