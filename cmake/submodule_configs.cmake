#------------------------------------------------------------------------------#
# © 2021-2023. Triad National Security, LLC. All rights reserved.  This
# program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by Triad
# National Security, LLC for the U.S.  Department of Energy/National
# Nuclear Security Administration. All rights in the program are
# reserved by Triad National Security, LLC, and the U.S. Department of
# Energy/National Nuclear Security Administration. The Government is
# granted for itself and others acting on its behalf a nonexclusive,
# paid-up, irrevocable worldwide license in this material to reproduce,
# prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.
#------------------------------------------------------------------------------#

macro(singularity_cmake_config pkg)
  # workaround
  # submodules will enable `BUILD_TESTING` by default
  # and this is used to disable it. 
  set(_storBT ${BUILD_TESTING})
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  if(${pkg} STREQUAL "Kokkos")
    set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
    if(SINGULARITY_USE_CUDA)
      set(Kokkos_ENABLE_CUDA ON CACHE BOOL "" FORCE)
      set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "" FORCE)
      set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE ON CACHE BOOL "" FORCE)
    endif()
  elseif(${pkg} STREQUAL "KokkosKernels")
    # Disable TPLs
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
  elseif(${pkg} STREQUAL "Eigen3")
    set(EIGEN_TEST_CXX11 OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
    #TODO: I don't like this as it can interfere with
    #     other projects, but I don't have a great solution
    #     to disable testing in here.
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  else()
    set(_PKG_DEFAULT_CONFIGURED ON)
  endif()
  if(_PKG_DEFAULT_CONFIGURED)
    singularity_msg(VERBOSE "[SUBMOD_CONFIG] \"${pkg}\" does not have a defined config, using defaults.")
    unset(_PKG_DEFAULT_CONFIGURED)
  else()
    #TODO: optional print of options
    singularity_msg(VERBOSE "[SUBMOD_CONFIG] \"${pkg}\" configured with defined config.")
  endif()
  set(BUILD_TESTING ${_storBT} CACHE BOOL "" FORCE)
  unset(_storBT)
endmacro()

