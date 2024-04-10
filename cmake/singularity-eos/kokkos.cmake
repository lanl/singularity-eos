macro(singularity_import_kokkos)

  set(Kokkos_ENABLE_SERIAL
      ON
      CACHE BOOL "" FORCE)
  if(SINGULARITY_USE_CUDA)
    set(Kokkos_ENABLE_CUDA
        ON
        CACHE BOOL "" FORCE)
    set(Kokkos_ENABLE_CUDA_LAMBDA
        ON
        CACHE BOOL "" FORCE)
  endif()

  if(Kokkos_ENABLE_CUDA)
    set(Kokkos_COMPILE_LANGUAGE CUDA)
  elseif(Kokkos_ENABLE_HIP)
    set(Kokkos_COMPILE_LANGUAGE HIP)
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
  # #TODO: shouldn't be needed target_compile_definitions(${target} PUBLIC
  # SPINER_USE_KOKKOS)

  # if(SINGULARITY_USE_KOKKOSKERNELS AND SINGULARITY_BUILD_CLOSURE)
  # target_link_libraries(${target} PUBLIC Kokkos::kokkoskernels) endif()
endmacro()

macro(singularity_import_kokkoskernels)
  set(KokkosKernels_ENABLE_TPL_BLAS
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_ENABLE_TPL_MKL
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_ENABLE_TPL_LAPACK
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_ENABLE_TPL_CUBLAS
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_ENABLE_TPL_CUSPARSE
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_ENABLE_TPL_MAGMA
      OFF
      CACHE BOOL "" FORCE)
  # Disable ETIs
  set(KokkosKernels_INST_COMPLEX_DOUBLE
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_COMPLEX_FLOAT
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_DOUBLE
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_EXECSPACE_OPENMP
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_EXECSPACE_SERIAL
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_EXECSPACE_THREADS
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_EXECSPACE_CUDA
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_FLOAT
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_LAYOUTLEFT
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_LAYOUTRIGHT
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_MEMSPACE_HOSTSPACE
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_OFFSET_INT
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_OFFSET_SIZE_T
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_ORDINAL_INT
      OFF
      CACHE BOOL "" FORCE)
  set(KokkosKernels_INST_ORDINAL_INT
      OFF
      CACHE BOOL "" FORCE)
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
