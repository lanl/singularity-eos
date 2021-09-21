#------------------------------------------------------------------------------#
# © 2021. Triad National Security, LLC. All rights reserved.  This
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
<<<<<<< HEAD
cmake_minimum_required(VERSION 3.19)
=======

cmake_minimum_required(VERSION 3.14)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# options
option (SINGULARITY_USE_HDF5 "Pull in hdf5" ON)
option (SINGULARITY_USE_FORTRAN "Enable fortran bindings" ON)
option (SINGULARITY_USE_KOKKOS "Use Kokkos for portability" OFF)
option (SINGULARITY_USE_EOSPAC "Pull in eospac" OFF)
option (SINGULARITY_USE_CUDA "Enable cuda support" OFF)
option (SINGULARITY_USE_KOKKOSKERNELS
  "Use kokkos-kernels for linear algebra" OFF)
option (SINGULARITY_INVERT_AT_SETUP "Use eospacs preinverted tables" OFF)
option (SINGULARITY_BUILD_TESTS "Compile tests" OFF)
option (SINGULARITY_BUILD_SESAME2SPINER "Compile sesame2spiner" OFF)
option (SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER "Compile stellarcollapse2spiner" OFF)
option (SINGULARITY_BETTER_DEBUG_FLAGS
  "Better debug flags for singularity" ON)
option (SINGULARITY_HIDE_MORE_WARNINGS "hide more warnings" OFF)
option (SINGULARITY_SUBMODULE_MODE "Submodule mode" OFF)
option (SINGULARITY_BUILD_CLOSURE "Mixed cell closure" ON)
option (SINGULARITY_TEST_SESAME "Test the Sesame table readers" OFF)
option (SINGULARITY_TEST_STELLAR_COLLAPSE "Test the stellar collapse table readers" OFF)
option (SINGULARITY_EOS_SKIP_EXTRAP "Turn off eospac extrap checks" OFF)
option (SINGULARITY_USE_SINGLE_LOGS "Use single precision logs. Can harm accuracy." OFF)
option (SINGULARITY_FMATH_USE_ORDER_4 "4th order interpolant for fast logs. Default is 7th order." OFF)
option (SINGULARITY_FMATH_USE_ORDER_5 "5th order interpolant for fast logs. Default is 7th order." OFF)

if (SINGULARITY_SUBMODULE_MODE)
  set(SINGULARITY_BETTER_DEBUG_FLAGS OFF CACHE BOOL "" FORCE)
  set(SINGULARITY_HIDE_MORE_WARNINGS ON CACHE BOOL "" FORCE)
endif()

# if cuda, need to set cxx compiler to nvcc_wrapper before compilers
# are checked
if(SINGULARITY_USE_KOKKOS AND SINGULARITY_USE_CUDA)
  # set nvcc wrapper default compiler
  if(NOT "$ENV{CXX}x" STREQUAL "x" AND
     "$ENV{NVCC_WRAPPER_DEFAULT_COMPILER}x" STREQUAL "x")
    set(ENV{NVCC_WRAPPER_DEFAULT_COMPILER} "$ENV{CXX}")
  endif()
  # set necessary kokkos build options if building inline
  if(NOT SINGULARITY_KOKKOS_INSTALL_DIR)
    set(NVCC_WRAPPER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/utils/kokkos/bin CACHE STRING "")
  else()
    set(NVCC_WRAPPER_DIR ${SINGULARITY_KOKKOS_INSTALL_DIR}/bin CACHE STRING "")
  endif()
  set(CMAKE_CXX_COMPILER ${NVCC_WRAPPER_DIR}/nvcc_wrapper CACHE STRING "")
endif()
>>>>>>> main

project(singularity-eos VERSION 1.5.0)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# Patches variant to be compatible with cuda
# Assumes "patch" is present on system
message(STATUS "Patching mpark::variant to support GPUs")
execute_process(COMMAND patch -N -s -V never
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/variant/include/mpark/variant.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/cuda_compatibility.patch)


add_library(${PROJECT_NAME}::flags INTERFACE IMPORTED GLOBAL)

include(SetupCustomMessage)
include(ProcessConfiguration)
include(SetupCompilers)
include(SetupFlags)

add_subdirectory(${PROJECT_SOURCE_DIR}/singularity-eos)

include(SetupInstall)

if(SINGULARITY_BUILD_TESTS)
  include(CTest)
  add_subdirectory(test)
endif()

set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")

include(CPack)

# feature_summary(WHAT ALL
#                 INCLUDE_QUIET_PACKAGES
#                 DESCRIPTION "Enabled Features:"
#                 VAR enabledFeaturesText)
# message(STATUS "${enabledFeaturesText}")
