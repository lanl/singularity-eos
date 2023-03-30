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

cmake_minimum_required(VERSION 3.20)

project(singularity-eos VERSION 1.7.0 LANGUAGES NONE)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(SINGULARITY_GOLDFILES_VERSION "goldfiles-1.6.2")
set(SINGULARITY_GOLDFILE_HASH
  33ef74b29937cc1347b525f72662933dd2b0556550f6541f97fc1de8a01c3c2a)

#------------------------------------------------------------------------------#
# Options
#------------------------------------------------------------------------------#
option (SINGULARITY_USE_HDF5 "Pull in hdf5" ON)
option (SINGULARITY_USE_FORTRAN "Enable fortran bindings" ON)
option (SINGULARITY_USE_KOKKOS "Use Kokkos for portability" OFF)
option (SINGULARITY_USE_EOSPAC "Pull in eospac" OFF)
option (SINGULARITY_USE_CUDA "Enable cuda support" OFF)
option (SINGULARITY_USE_KOKKOSKERNELS
  "Use kokkos-kernels for linear algebra" OFF)
option (SINGULARITY_BUILD_PYTHON "Compile Python bindings" OFF)
option (SINGULARITY_BUILD_TESTS "Compile tests" OFF)
option (SINGULARITY_BUILD_EXAMPLES "Compile examples" OFF)
option (SINGULARITY_BUILD_SESAME2SPINER "Compile sesame2spiner" OFF)
option (SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER "Compile stellarcollapse2spiner" OFF)
option (SINGULARITY_BETTER_DEBUG_FLAGS
  "Better debug flags for singularity" ON)
option (SINGULARITY_HIDE_MORE_WARNINGS "hide more warnings" OFF)
option (SINGULARITY_BUILD_CLOSURE "Mixed cell closure" ON)
option (SINGULARITY_TEST_SESAME "Test the Sesame table readers" OFF)
option (SINGULARITY_TEST_STELLAR_COLLAPSE "Test the stellar collapse table readers" OFF)
option (SINGULARITY_TEST_PYTHON "Test the Python bindings" OFF)
option (SINGULARITY_USE_SINGLE_LOGS "Use single precision logs. Can harm accuracy." OFF)
option (SINGULARITY_USE_TRUE_LOG_GRIDDING "Use grids that conform to log spacing." OFF)

option (SINGULARITY_SUBMODULE_MODE "Submodule mode" OFF)
option (SINGULARITY_CUSTOM_KOKKOS "Use custom kokkos dir" OFF)

option (SINGULARITY_PATCH_MPARK_VARIANT "Apply GPU patch to mpark-variant submodule" ON)

#------------------------------------------------------------------------------#
# singularity-eos Library
#------------------------------------------------------------------------------#

include(singularity/mpark_variant)
include(singularity/Eigen3)
include(singularity/eospac)
include(singularity/hdf5)
include(singularity/kokkos)
include(singularity/spiner)
include(singularity/ports-of-call)

add_library(singularity-eos)
add_library(singularity-eos::singularity-eos ALIAS singularity-eos)

#?
target_include_directories(singularity-eos PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>)

#------------------------------------------------------------------------------#
# Compiler & language setup
#------------------------------------------------------------------------------#
enable_language(CXX)
include(CMakeDetermineCXXCompiler)
if(SINGULARITY_USE_FORTRAN)
  enable_language(Fortran)
  include(CMakeDetermineFortranCompiler)
endif()

if (SINGULARITY_BUILD_PYTHON)
  # need position independent code for Python shared library
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# checks if this is our build, or we've 
# been imported via `add_subdirectory`
# NB: this should make the `option(SINGULARITY_SUBMODULE_MODE ...)` unnecessary
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
else()
  if(NOT SINGULARITY_SUBMODULE_MODE)
    message(WARNING "Detected that `singularity-eos` is a subproject, but `SINGULARITY_SUBMODULE_MODE` is not set or \
    is set to OFF. We will force this option on.")
    set(SINGULARITY_SUBMODULE_MODE ON CACHE BOOL "" FORCE)
  endif()
endif()

#------------------------------------------------------------------------------#
# De-thaw some options
#------------------------------------------------------------------------------#

# using this space for "immediate" logic that can be avoid messages getting
# lost inside the brrrrrrrr of cmake.
if(SINGULARITY_USE_KOKKOSKERNELS AND NOT SINGULARITY_USE_KOKKOS)
  message(FATAL_ERROR "Kokkos is required for KokkosKernels.")
endif()

if(SINGULARITY_USE_CUDA)
  if(NOT SINGULARITY_USE_KOKKOS)
    message(FATAL_ERROR "CUDA without Kokkos is not currently supported.\n\
    If available, reconfigure with '-DSINGULARITY_USE_KOKKOS=ON'. Otherwise, disable CUDA offloading '-DSINGULARITY_USE_CUDA=OFF'")
  endif()
  if(SINGULARITY_BUILD_CLOSURE AND NOT SINGULARITY_USE_KOKKOSKERNELS)
    message(FATAL_ERROR "Mixed cell closures require linear algebra. With CUDA offload, this must must utilize KokkosKernels.\n\
    If available, reconfigure with '-DSINGULARITY_USE_KOKKOSKERNELS=ON'.\n\
    Otherwise, disable CUDA offloading '-DSINGULARITY_USE_CUDA=OFF' or disable mixed cell closures with '-DSINGULARITY_BUILD_CLOSURE=OFF'.")
  endif()
endif()

if (SINGULARITY_USE_EOSPAC)
  if (NOT SINGULARITY_USE_HDF5)
    message(WARNING "EOSPAC enabled without hdf5.")
  endif()
  if (SINGULARITY_USE_CUDA)
    message(WARNING "EOSPAC not supported with cuda. You may encounter errors.")
  endif()
endif()

#------------------------------------------------------------------------------#
# Process some options
#------------------------------------------------------------------------------#

if (SINGULARITY_SUBMODULE_MODE)
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  set(CMAKE_FIND_USE_PACKAGE_REGISTRY OFF CACHE BOOL "" FORCE)
  set(CMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY OFF CACHE BOOL "" FORCE)
  set(SINGULARITY_BETTER_DEBUG_FLAGS OFF CACHE BOOL "" FORCE)
  set(SINGULARITY_HIDE_MORE_WARNINGS ON CACHE BOOL "" FORCE)
endif()

# side-projects
if(SINGULARITY_USE_HDF5)
  if(SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER)
    add_subdirectory(stellarcollapse2spiner)
  endif()
  if(SINGULARITY_USE_EOSPAC AND SINGULARITY_BUILD_SESAME2SPINER)
    add_subdirectory(sesame2spiner)
  endif()
  if(SINGULARITY_BUILD_EXAMPLES)
    add_subdirectory(example)
  endif()
endif()

if (SINGULARITY_BUILD_PYTHON)
  add_subdirectory(python)
endif()

if(SINGULARITY_USE_SINGLE_LOGS)
  target_compile_definitions(singularity-eos PUBLIC SINGULARITY_USE_SINGLE_LOGS)
endif()

#------------------------------------------------------------------------------#
# Handle dependencies
#------------------------------------------------------------------------------#

singularity_enable_mpark_variant(singularity-eos)

singularity_enable_ports_of_call(singularity-eos)

#----------------------------------------------------------------------------#
# Optional dependencies
#----------------------------------------------------------------------------#

# need this for consistant Catch2 up/downstream.
# depending on how this develops, this could be wrapped
# into the dependcy logic (below)

if(SINGULARITY_BUILD_TESTS)
  if (NOT TARGET Catch2::Catch2)
    find_package(Catch2 QUIET)
  endif()

  if(NOT Catch2_FOUND)
    message(STATUS "Fetching Catch2 as needed")
    # idiomatic FetchContent
    include(FetchContent)
    FetchContent_Declare(
      Catch2
      GIT_REPOSITORY https://github.com/catchorg/Catch2.git
      # MAUNEYC: On my side, this version works. The ones
      #           pinned onto upstream codes do not like to play
      #           with clang and/or newer compilers.
      GIT_TAG        v2.13.7)
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/contrib)
  endif()
endif()


if(SINGULARITY_USE_KOKKOS)
  singularity_enable_kokkos(singularity-eos)
endif()

if(SINGULARITY_USE_HDF5)
  singularity_enable_hdf5(singularity-eos)
endif()

if(SINGULARITY_USE_SPINER)
  singularity_enable_spiner(singularity-eos)
endif()

if(SINGULARITY_USE_EOSPAC)
  singularity_enable_eospac(singularity-eos)
endif()

#TODO: this feels out-of-place here, though i haven't crawled into the source yet to see
# what the point is.
if (SINGULARITY_TEST_SESAME)
  target_compile_definitions(singularity-eos PRIVATE SINGULARITY_TEST_SESAME)
endif()
if (SINGULARITY_BUILD_CLOSURE)
  target_compile_definitions(singularity-eos PRIVATE SINGULARITY_BUILD_CLOSURE)
endif()
###########################OLD
#------------------------------------------------------------------------------#
# singularity-eos library
#------------------------------------------------------------------------------#

# this subdirectory populates `EOS_HEADERS/EOS_SRCS`
# NOTE: these include path prefixes of subdirectories on files (e.g. eos/eos.hpp)
# see singularity-eos/CMakeLists.txt

add_subdirectory(singularity-eos)

foreach(_header ${EOS_HEADERS})
  list(APPEND _install_headers ${_header})
  list(APPEND _headers singularity-eos/${_header})
endforeach()

foreach(_mod ${EOS_MODS})
  list(APPEND _install_mods ${CMAKE_CURRENT_BINARY_DIR}/${_mod})
endforeach()

foreach(_src ${EOS_SRCS})
  list(APPEND _srcs singularity-eos/${_src})
endforeach()

target_sources(singularity-eos PRIVATE ${_srcs} ${_headers})

# make sure .mods are placed in build path, and installed
# along with includes
if (SINGULARITY_USE_FORTRAN)
  set_target_properties(singularity-eos PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  target_include_directories(singularity-eos INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include/singularity-eos/eos>
  )
endif() # SINGULARITY_USE_FORTRAN


target_include_directories(singularity-eos
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)


# plug in collected includes/libs/definitions

target_include_directories(singularity-eos
  PUBLIC
    ${SINGULARITY_PUBLIC_INCS}
  PRIVATE
    ${SINGULARITY_PRIVATE_INCS}
)

target_link_libraries(singularity-eos
  PUBLIC
    ${SINGULARITY_PUBLIC_LIBS}
  PRIVATE
    ${SINGULARITY_PRIVATE_LIBS}
)

target_compile_definitions(singularity-eos
  PUBLIC
    ${SINGULARITY_PUBLIC_DEFINES}
  PRIVATE
    ${SINGULARITY_PRIVATE_DEFINES}
)

#------------------------------------------------------------------------------#
# build flags
#------------------------------------------------------------------------------#

set(with_cxx "$<COMPILE_LANGUAGE:CXX>")
set(with_xlcxx  "$<COMPILE_LANG_AND_ID:CXX,XL>")
set(build_debug "$<CONFIG:Debug>")
set(build_release "$<CONFIG:Release>")
set(with_cuda "$<BOOL:${SINGULARITY_USE_CUDA}>")
set(with_warnings "$<NOT:$<BOOL:${SINGULARITY_HIDE_MORE_WARNINGS}>>")
set(with_better_flags "$<BOOL:${SINGULARITY_BETTER_DEBUG_FLAGS}>")

set(xlfix "$<${with_xlcxx}:-std=c++1y;-qxflag=disable__cplusplusOverride>")

target_compile_options(singularity-eos
  PRIVATE
    ${xlfix}
  PUBLIC
    $<${with_cuda}:
      $<${with_cxx}:
        --expt-relaxed-constexpr
        $<${with_warnings}:
          -Xcudafe;--diag_suppress=esa_on_defaulted_function_ignored
        > #with_warnings
        $<${build_debug}:
          $<${with_better_flags}:
            -G;-lineinfo
          > #better_flags
        > #debug
      > #cxx
      $<${build_release}:
        -use_fast_math
      > #release
    > #cuda
)

target_link_options(singularity-eos
  PRIVATE
    ${xlfix}
)

include(GNUInstallDirs)

#----------------------------------------------------------------------------#
# Export/install
#----------------------------------------------------------------------------#

# placed into a module file.
# NB: exporting is done too, so 
#     even submodule mode should
#     include this file.
include(cmake/install.cmake)

#----------------------------------------------------------------------------#
# Testing
#----------------------------------------------------------------------------#

#NOTE: `include(CTest)` comes "pre-set" with options and modules
# to enable tests. It creates a "BUILD_TESTING" option, and runs
# `enable_testing()`. See https://cmake.org/cmake/help/latest/module/CTest.html
if(SINGULARITY_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

#----------------------------------------------------------------------------#
# Misc.
#----------------------------------------------------------------------------#

# clang format
include(cmake/Format.cmake)


