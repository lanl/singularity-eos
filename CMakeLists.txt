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
option (SINGULARITY_BUILD_TESTS "Compile tests" OFF)
option (SINGULARITY_BUILD_EXAMPLES "Compile examples" OFF)
option (SINGULARITY_BUILD_SESAME2SPINER "Compile sesame2spiner" OFF)
option (SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER "Compile stellarcollapse2spiner" OFF)
option (SINGULARITY_BETTER_DEBUG_FLAGS
  "Better debug flags for singularity" ON)
option (SINGULARITY_HIDE_MORE_WARNINGS "hide more warnings" OFF)
option (SINGULARITY_SUBMODULE_MODE "Submodule mode" OFF)
option (SINGULARITY_BUILD_CLOSURE "Mixed cell closure" ON)
option (SINGULARITY_TEST_SESAME "Test the Sesame table readers" OFF)
option (SINGULARITY_TEST_STELLAR_COLLAPSE "Test the stellar collapse table readers" OFF)
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

project(singularity-eos VERSION 1.5.0)

#------------------------------------------------------------------------------#
# Compiler & language setup
#------------------------------------------------------------------------------#
enable_language(CXX)
include(CMakeDetermineCXXCompiler)
if(SINGULARITY_USE_FORTRAN)
  enable_language(Fortran)
  # $CMM/cmake-idiomatic-install
  # $CMM added to determine Fortran compiler
  include(CMakeDetermineFortranCompiler)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# $CMM/cmake-idiomatic-install
# $CMM prepare some generator-string snippets
set(with_cxx "$<COMPILE_LANGUAGE:CXX>")
set(build_release "$<CONFIG:Release>")
set(build_debug "$<CONFIG:Debug>")
set(on_xl "$<COMPILER_ID:XL>")

set(CMAKE_EXPORT_PACKAGE_REGISTRY "ON")

# Patches variant to be compatible with cuda
# Assumes "patch" is present on system
message(STATUS "Patching mpark::variant to support GPUs")
execute_process(COMMAND patch -N -s -V never
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/variant/include/mpark/variant.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/cuda_compatibility.patch
)
#------------------------------------------------------------------------------#

set(_tar singularity-eos)
set(_ns singularity-eos)

add_subdirectory(singularity-eos)
foreach(_header ${EOS_HEADERS})
  list(APPEND _install_headers ${_header})
  list(APPEND _headers singularity-eos/${_header})
endforeach()

foreach(_src ${EOS_SRCS})
  list(APPEND _srcs singularity-eos/${_src})
endforeach()

add_library(${_tar} ${_srcs} ${_headers})
add_library(${_ns}::${_tar} ALIAS ${_tar})

if (SINGULARITY_USE_FORTRAN)
  set_target_properties(${_tar} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  target_include_directories(${_tar} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>
  )
endif()

target_include_directories(${_tar}
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

# TODO: Is this generator expression correct/needed?
#target_include_directories(singularity-eos::flags
#                           INTERFACE
#                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
#                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utils>)

# xl fix
#target_compile_options(singularity-eos::flags INTERFACE
#                       $<$<COMPILE_LANG_AND_ID:CXX,XL>:-std=c++1y;-qxflag=disable__cplusplusOverride>)
#target_link_options(singularity-eos::flags INTERFACE
#                   $<$<COMPILE_LANG_AND_ID:CXX,XL>:-std=c++1y;-qxflag=disable__cplusplusOverride>)

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
    message(WARNING "EOSPAC not supported with cuda. You may encounter errors.")
  endif()
endif()

if (SINGULARITY_TEST_SESAME)
   target_compile_definitions(${_tar} INTERFACE SINGULARITY_TEST_SESAME)
endif()

# only use hdf5 if spiner
if (SINGULARITY_USE_HDF5)
  target_compile_definitions(${_tar} INTERFACE
                             SPINER_USE_HDF)
  if(SINGULARITY_HDF5_INSTALL_DIR)
    list(APPEND CMAKE_PREFIX_PATH "${SINGULARITY_HDF5_INSTALL_DIR}")
  endif()
  find_package(HDF5 COMPONENTS C HL)
  if(HDF5_FOUND)
    #    set_target_properties(${_tar}
    #  PROPERTIES
    #    INTERFACE_LINK_LIBRARIES "${HDF5_LIBRARIES};${HDF5_HL_LIBRARIES}"
    #  INTERFACE_COMPILE_DEFINITIONS "SINGULARITY_USE_HDF5"
    #  INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIRS}")
    target_link_libraries(${_tar} INTERFACE "${HDF5_LIBRARIES};${HDF5_HL_LIBRARIES}")
    target_compile_definitions(${_tar} INTERFACE "SINGULARITY_USE_HDF5")
    target_include_directories(${_tar} INTERFACE "${HDF5_INCLUDE_DIRS}")
    if(HDF5_IS_PARALLEL)
      if(SINGULARITY_MPI_INSTALL_DIR)
        list(APPEND CMAKE_PREFIX_PATH "${SINGULARITY_MPI_INSTALL_DIR}")
      endif()
      find_package(MPI COMPONENTS CXX)
      if(MPI_FOUND)
        target_include_directories(${_tar} INTERFACE "${MPI_CXX_INCLUDE_DIRS}")
      endif()
    endif()
  else()
    message(FATAL_ERROR "HDF5 was requested but not found. Can be disabled with -DSINGULARITY_USE_HDF5=OFF")
  endif()
endif()

if (SINGULARITY_USE_KOKKOS)
  target_compile_definitions(${_tar} INTERFACE
                             PORTABILITY_STRATEGY_KOKKOS)
  # Import Kokkos if not already available as a build target
  if (NOT TARGET Kokkos::kokkos)
    if (SINGULARITY_KOKKOS_INSTALL_DIR)
      list(APPEND CMAKE_PREFIX_PATH "${SINGULARITY_KOKKOS_INSTALL_DIR}")
      set(Kokkos_ROOT "${SINGULARITY_KOKKOS_INSTALL_DIR}/lib64/cmake/Kokkos")
    endif()
    find_package(Kokkos)
    if(NOT Kokkos_FOUND)
      if(SINGULARITY_USE_CUDA)
        set(Kokkos_ENABLE_CUDA ON CACHE BOOL "" FORCE)
        set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
        set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "" FORCE)
        set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE ON CACHE BOOL "" FORCE)
      endif()
      add_subdirectory(${PROJECT_SOURCE_DIR}/utils/kokkos)
    endif()
    get_target_property(KOKKOS_COMPILE_OPTIONS Kokkos::kokkoscore
        INTERFACE_COMPILE_OPTIONS
    )
    target_link_libraries(${_tar} PUBLIC Kokkos::kokkos)
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
  target_link_libraries(${_tar} PUBLIC
                        Kokkos::kokkoskernels)
  target_compile_definitions(${_tar} INTERFACE
                             SINGULARITY_USE_KOKKOSKERNELS)
else()
  # Eigen
  if (NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 QUIET)
    if(NOT Eigen3_FOUND)
        set(EIGEN_TEST_CXX11 OFF CACHE BOOL "" FORCE)
        set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
        set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
        set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
        add_subdirectory(${PROJECT_SOURCE_DIR}/utils/eigen)
      endif()
    target_link_libraries(${_tar} PUBLIC Eigen3::Eigen)
  endif()
endif()

if (SINGULARITY_USE_EOSPAC)
  add_subdirectory(eospac-wrapper)
  target_link_libraries(${_tar}
    PUBLIC
    eospac-wrapper)
  target_compile_definitions(${_tar}
    PUBLIC
    SINGULARITY_USE_EOSPAC)
endif ()

# fast math
if (SINGULARITY_USE_SINGLE_LOGS)
  target_compile_definitions(${_tar} INTERFACE
      SINGULARITY_USE_SINGLE_LOGS)
endif()
if (SINGULARITY_FMATH_USE_ORDER_4 AND SINGULARITY_FMATH_USE_ORDER_5)
  message(FATAL_ERROR "Order 4 and order 5 interpolation both specified. Please specify only one or zero. If no specification is made, order 7 will be used.")
endif()
if (SINGULARITY_FMATH_USE_ORDER_4)
  target_compile_definitions(${_tar} INTERFACE
    SINGULARITY_FMATH_USE_ORDER_4)
endif()
if (SINGULARITY_FMATH_USE_ORDER_5)
  target_compile_definitions(${_tar} INTERFACE
    SINGULAIRTY_FMATH_USE_ORDER_5)
endif()

# required package ports of call
# same as spiner
include(GNUInstallDirs)
add_library(PortsofCall INTERFACE IMPORTED)

# temporary
# add submodule or use subdir
list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/../installs.singularity/spiner.main.defaultconfig/lib/cmake/spiner")
find_package(spiner REQUIRED)

target_link_libraries(${_tar} PUBLIC spiner::spiner)

target_include_directories(PortsofCall INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utils>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/ports-of-call>
    )

target_link_libraries(${_tar} PUBLIC PortsofCall)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/utils/ports-of-call
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ports-of-call)


# catch2 - check for env catch first, otherwise fetch
if(SINGULARITY_BUILD_TESTS)
  #  add_subdirectory(${PROJECT_SOURCE_DIR}/utils/spiner/Catch2)
endif()

if(SINGULARITY_BUILD_TESTS)
  include(CTest)
  include(${PROJECT_SOURCE_DIR}/utils/spiner/Catch2/contrib/Catch.cmake)
endif()

if (SINGULARITY_USE_HDF5 AND SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER)
   add_subdirectory(${PROJECT_SOURCE_DIR}/stellarcollapse2spiner)
   install(TARGETS stellarcollapse2spiner DESTINATION bin)
endif()

if (SINGULARITY_BUILD_EXAMPLES)
  add_subdirectory(${PROJECT_SOURCE_DIR}/example)
endif()

if(NOT DEFINED _libdir)
  include(GNUInstallDirs)
  set(_libdir "${CMAKE_INSTALL_LIBDIR}")
endif(NOT DEFINED _libdir)

install(
  TARGETS ${_tar} 
  EXPORT ${_tar}Targets
  DESTINATION ${_libdir}
)

install(
  EXPORT ${_tar}Targets
  FILE  ${_tar}Targets.cmake
  NAMESPACE "${_ns}::"
  DESTINATION ${_libdir}/cmake/${_tar}
)

foreach(file ${_install_headers})
  get_filename_component(DIR ${file} DIRECTORY)
  install(FILES ${_tar}/${file}
  DESTINATION include/${_tar}/${DIR})
endforeach()

#----------------------------------------------------------------------------#
# Packaging
#----------------------------------------------------------------------------#

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_SOURCE_DIR}/config/${_tar}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${_tar}Config.cmake
  INSTALL_DESTINATION ${_libdir}/cmake/${_tar}
)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${_tar}Config.cmake
  DESTINATION ${_libdir}/cmake/${_tar}
)

#----------------------------------------------------------------------------#
# Local export
#----------------------------------------------------------------------------#

export(
  EXPORT    ${_tar}Targets
  FILE      ${CMAKE_CURRENT_BINARY_DIR}/cmake/${_tar}Targets.cmake
  NAMESPACE ${_ns}::
)

export(PACKAGE ${_tar})

if(SINGULARITY_USE_HDF5)
  if(SINGULARITY_USE_EOSPAC AND SINGULARITY_BUILD_SESAME2SPINER)
    add_subdirectory(sesame2spiner)
  endif()
endif()

if(SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER)
  add_subdirectory(stellarcollapse2spiner)
endif()

if(SINGULARITY_BUILD_TESTS)
  include(CTest)
  add_subdirectory(test)
endif()

# clang format
include(cmake/Format.cmake)

include(CPack)
