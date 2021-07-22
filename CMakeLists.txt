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
if (SINGULARITY_SUBMODULE_MODE)
  set(SINGULARITY_BETTER_DEBUG_FLAGS OFF CACHE BOOL "" FORCE)
  set(SINGULARITY_HIDE_MORE_WARNINGS ON CACHE BOOL "" FORCE)
endif()

# if cuda, need to set cxx compiler to nvcc_wrapper before compilers
# are checked
if(SINGULARITY_USE_KOKKOS AND SINGULARITY_USE_CUDA)
  # json stuff won't compile with nvcc wrapper
  if(SINGULARITY_BUILD_SESAME2SPINER)
    message(FATAL_ERROR "sesame2spiner cannot be built with nvcc. If you wish to buld it, disable CUDA with -DSINGULARITY_USE_CUDA=OFF. Otherwise disable sesame2spiner with -DSINGULARITY_BUILD_SESAME2SPINER=OFF.")
  endif()
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

project(singularity-eos VERSION 1.5.0)

if(SINGULARITY_USE_FORTRAN)
  enable_language(CXX Fortran)
else()
  enable_language(CXX)
endif()
include(CMakeDetermineCXXCompiler)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Patches variant to be compatible with cuda
# Assumes "patch" is present on system
message(STATUS "Patching mpark::variant to support GPUs")
execute_process(COMMAND patch -N -s -V never
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/variant/include/mpark/variant.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/cuda_compatibility.patch)

add_library (singularity-eos::flags INTERFACE IMPORTED GLOBAL)
add_library (singularity-eos::libs INTERFACE IMPORTED)

# TODO: Is this generator expression correct/needed?
target_include_directories(singularity-eos::flags
                           INTERFACE
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utils>)

# xl fix
target_compile_options(singularity-eos::flags INTERFACE
                       $<$<COMPILE_LANG_AND_ID:CXX,XL>:-std=c++1y;-qxflag=disable__cplusplusOverride>)
target_link_options(singularity-eos::flags INTERFACE
                   $<$<COMPILE_LANG_AND_ID:CXX,XL>:-std=c++1y;-qxflag=disable__cplusplusOverride>)

#target_compile_features(singularity-eos::flags INTERFACE cxx_std_14)
#target_compile_features(singularity-eos::libs INTERFACE cxx_std_14)

# check for currently incompatible option
if(SINGULARITY_USE_CUDA AND NOT SINGULARITY_USE_KOKKOS)
  message(FATAL_ERROR "Cuda without kokkos is not currently supported")
endif()

if (SINGULARITY_USE_EOSPAC)
  if (NOT SINGULARITY_USE_HDF5)
    message(WARNING "EOSPAC enabled without hdf5.")
  endif()
  if (SINGULARITY_USE_CUDA)
    message(FATAL_ERROR "EOSPAC not supported with cuda")
  endif()
endif()

if (SINGULARITY_TEST_SESAME)
   target_compile_definitions(singularity-eos::flags INTERFACE SINGULARITY_TEST_SESAME)
endif()

if (SINGULARITY_USE_HDF5)
  target_compile_definitions(singularity-eos::flags INTERFACE
                             SPINER_USE_HDF)
  if(SINGULARITY_HDF5_INSTALL_DIR)
    list(APPEND CMAKE_PREFIX_PATH "${SINGULARITY_HDF5_INSTALL_DIR}")
  endif()
  find_package(HDF5 COMPONENTS C HL)
  if(HDF5_FOUND)
    add_library(singularity-eos::hdf5 INTERFACE IMPORTED GLOBAL)
    set_target_properties(singularity-eos::hdf5 PROPERTIES
      INTERFACE_LINK_LIBRARIES "${HDF5_LIBRARIES};${HDF5_HL_LIBRARIES}"
      INTERFACE_COMPILE_DEFINITIONS "SINGULARITY_USE_HDF5"
      INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIRS}")
    if(HDF5_IS_PARALLEL)
      if(SINGULARITY_MPI_INSTALL_DIR)
        list(APPEND CMAKE_PREFIX_PATH "${SINGULARITY_MPI_INSTALL_DIR}")
      endif()
      find_package(MPI COMPONENTS CXX)
      if(MPI_FOUND)
        target_include_directories(singularity-eos::libs INTERFACE "${MPI_CXX_INCLUDE_DIRS}")
      endif()
    endif()
  else()
    message(FATAL_ERROR "HDF5 was requested but not found. Can be disabled with -DSINGULARITY_USE_HDF5=OFF")
  endif()
  target_link_libraries(singularity-eos::libs INTERFACE singularity-eos::hdf5)
endif()

if (SINGULARITY_USE_KOKKOS)
  target_compile_definitions(singularity-eos::flags INTERFACE
                             PORTABILITY_STRATEGY_KOKKOS)
  # Import Kokkos if not already available as a build target
  if (NOT TARGET Kokkos::kokkos)
    if (SINGULARITY_KOKKOS_INSTALL_DIR)
      list(APPEND CMAKE_PREFIX_PATH "${SINGULARITY_KOKKOS_INSTALL_DIR}")
      set(Kokkos_ROOT "${SINGULARITY_KOKKOS_INSTALL_DIR}/lib64/cmake/Kokkos")
      find_package(Kokkos REQUIRED)
    else()
      if(SINGULARITY_USE_CUDA)
        set(Kokkos_ENABLE_CUDA ON CACHE BOOL "" FORCE)
        set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
        set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "" FORCE)
        set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE ON CACHE BOOL "" FORCE)
      endif()
      add_subdirectory(${PROJECT_SOURCE_DIR}/utils/kokkos)
      install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/utils/kokkos"
	DESTINATION "include"
	FILES_MATCHING PATTERN "*.hpp"
	)
    endif()
  endif()
  # cuda specific flags
  if(SINGULARITY_USE_CUDA)
    if(SINGULARITY_HIDE_MORE_WARNINGS)
      target_compile_options(
        singularity-eos::flags
         INTERFACE # Generator expression shamelessly copied from EAP
         "$<$<COMPILE_LANGUAGE:CXX>:--expt-relaxed-constexpr;>"
           )
    else()
      target_compile_options(
        singularity-eos::flags
         INTERFACE # Generator expression shamelessly copied from EAP
         "$<$<COMPILE_LANGUAGE:CXX>:--expt-relaxed-constexpr;-Xcudafe;--diag_suppress=esa_on_defaulted_function_ignored;>"
           )
    endif()
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
      target_compile_options(singularity-eos::flags INTERFACE
                             -use_fast_math)
    elseif(CMAKE_BUILD_TYPE STREQUAL "Debug" AND SINGULARITY_BETTER_DEBUG_FLAGS)
      target_compile_options(
        singularity-eos::flags
        INTERFACE
        "$<$<COMPILE_LANGUAGE:CXX>:-G;-lineinfo;>"
        )
    endif()
  endif()
  target_link_libraries(singularity-eos::libs INTERFACE Kokkos::kokkos)
endif ()

# linear algebra library
if(SINGULARITY_KOKKOSKERNELS_SUB_DIR)
  set(SINGULARITY_USE_KOKKOSKERNELS ON CACHE BOOL "" FORCE)
endif()
if(SINGULARITY_KOKKOSKERNELS_INSTALL_DIR)
  set(SINGULARITY_USE_KOKKOSKERNELS ON CACHE BOOL "" FORCE)
endif()
if(SINGULARITY_USE_CUDA AND NOT SINGULARITY_USE_KOKKOSKERNELS)
  message(WARNING "Kokkos Kernels is required if you want to use closures on GPU.")
endif()
if(SINGULARITY_USE_KOKKOSKERNELS)
  if(NOT DEFINED SINGULARITY_USE_KOKKOS)
    message(FATAL_ERROR "Kokkos is a required dependency for KokkosKernels")
  endif()
  # kokkos kernels
  if(SINGULARITY_KOKKOSKERNELS_SUB_DIR)
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
    add_subdirectory(${SINGULARITY_KOKKOSKERNELS_SUB_DIR}
                     ${PROJECT_BINARY_DIR}/utils/kokkos-kernels)
  else()
    if(SINGULARITY_KOKKOSKERNELS_INSTALL_DIR)
      set(KokkosKernels_ROOT "${SINGULARITY_KOKKOSKERNELS_INSTALL_DIR}")
    endif()
    find_package(KokkosKernels REQUIRED)
  endif()
  target_link_libraries(singularity-eos::flags INTERFACE
                        Kokkos::kokkoskernels)
  target_compile_definitions(singularity-eos::flags INTERFACE
                             SINGULARITY_USE_KOKKOSKERNELS)
else()
  # Eigen
  if (NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 REQUIRED)
    target_link_libraries(singularity-eos::flags INTERFACE Eigen3::Eigen)
    install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/utils/eigen"
      DESTINATION "include"
      FILES_MATCHING PATTERN "*.hpp"
      )
  endif()
endif()

if (SINGULARITY_USE_EOSPAC)
  target_compile_definitions(singularity-eos::flags INTERFACE
                             SINGULARITY_USE_EOSPAC
                             EOSPAC_WARN)
  if (SINGULARITY_INVERT_AT_SETUP)
    target_compile_definitions(singularity-eos::flags INTERFACE
                               SINGULARITY_INVERT_AT_SETUP)
  endif ()
  if (SINGULARITY_EOS_SKIP_EXTRAP)
    target_compile_definitions(singularity-eos::flags INTERFACE
                               SINGULARITY_EOS_SKIP_EXTRAP)
  endif ()
  if (SINGULARITY_EOSPAC_INSTALL_DIR)
    list(APPEND CMAKE_PREFIX_PATH "${SINGULARITY_EOSPAC_INSTALL_DIR}")
    set(EOSPAC_INCLUDE_DIR "${SINGULARITY_EOSPAC_INSTALL_DIR}/include")
  endif()
  find_package(EOSPAC)
  target_link_libraries(singularity-eos::libs INTERFACE EOSPAC::eospac)
endif ()

# required package ports of call
find_package(PortsofCall REQUIRED)
target_link_libraries(singularity-eos::flags INTERFACE PortsofCall::PortsofCall)

# catch2
if(SINGULARITY_BUILD_TESTS)
  add_subdirectory(${PROJECT_SOURCE_DIR}/utils/spiner/Catch2)
endif()

# JSON
if(SINGULARITY_BUILD_SESAME2SPINER)
  set(JSON_BuildTests OFF CACHE INTERNAL "")
  set(JSON_Install OFF CACHE INTERNAL "")
  add_subdirectory(${PROJECT_SOURCE_DIR}/utils/json)
endif()

if(SINGULARITY_BUILD_TESTS)
  include(CTest)
  include(${PROJECT_SOURCE_DIR}/utils/spiner/Catch2/contrib/Catch.cmake)
endif()

add_subdirectory(${PROJECT_SOURCE_DIR}/singularity-eos/eos)
if (SINGULARITY_BUILD_TESTS)
  add_subdirectory(test)
endif()
# create the libsingularity-eos.a library
if(SINGULARITY_USE_FORTRAN)
  add_library(singularity-eos
              $<TARGET_OBJECTS:eos_f>
              $<TARGET_OBJECTS:eos_objs>)
  install(DIRECTORY "${CMAKE_BINARY_DIR}/singularity-eos"
    DESTINATION "include"
    FILES_MATCHING PATTERN "*.mod")
else()
  add_library(singularity-eos
              $<TARGET_OBJECTS:eos_objs>)
endif()

if (SINGULARITY_USE_HDF5 AND SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER)
   add_subdirectory(${PROJECT_SOURCE_DIR}/stellarcollapse2spiner)
   install(TARGETS stellarcollapse2spiner DESTINATION bin)
endif()

install(TARGETS eos singularity-eos DESTINATION lib)
if(SINGULARITY_BUILD_SESAME2SPINER)
  if(SINGULARITY_USE_EOSPAC)
    if(SINGULARITY_USE_HDF5)
      add_subdirectory(${PROJECT_SOURCE_DIR}/sesame2spiner)
      install(TARGETS sesame2spiner DESTINATION bin)
   endif()
  endif()
endif()

# These are always needed.
# TODO(JMM): Should always-needed utils be moved into a separate
# directory? Or the optional ones be moved to something like
# "external"?
# TODO(JMM): Should I write a function so I'm not repeating all these
# install calls?  More cmake wizardry than I have is required for
# that, I think.
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/singularity-eos"
  DESTINATION "include"
  FILES_MATCHING PATTERN "*.hpp"
  )
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/utils/fast-math"
  DESTINATION "include"
  FILES_MATCHING PATTERN "*.hpp"
  )
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/utils/herumi-fmath"
  DESTINATION "include"
  FILES_MATCHING PATTERN "*.hpp"
  )
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/utils/root-finding-1d"
  DESTINATION "include"
  FILES_MATCHING PATTERN "*.hpp"
  )
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/utils/sp5"
  DESTINATION "include"
  FILES_MATCHING PATTERN "*.hpp"
  )
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/utils/spiner"
  DESTINATION "include"
  FILES_MATCHING PATTERN "*.hpp"
  )
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/utils/variant"
  DESTINATION "include"
  FILES_MATCHING PATTERN "*.hpp"
  )
install(CODE
  "execute_process( \
    COMMAND ${CMAKE_COMMAND} -E create_symlink \
    ${CMAKE_INSTALL_PREFIX}/include/spiner/ports-of-call \
    ${CMAKE_INSTALL_PREFIX}/include/ports-of-call \
  )"
  )
