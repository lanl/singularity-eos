# ------------------------------------------------------------------------------#
# © 2021-2023. Triad National Security, LLC. All rights reserved.  This program
# was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
# National Laboratory (LANL), which is operated by Triad National Security, LLC
# for the U.S.  Department of Energy/National Nuclear Security Administration.
# All rights in the program are reserved by Triad National Security, LLC, and
# the U.S. Department of Energy/National Nuclear Security Administration. The
# Government is granted for itself and others acting on its behalf a
# nonexclusive, paid-up, irrevocable worldwide license in this material to
# reproduce, prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.
# ------------------------------------------------------------------------------#

# We want to stay near to latest, as more compilers and third-party libraries
# are supported. However, uptake of new CMake releases can be slow. 20230403:
# version 3.19 is available on all testing platforms
cmake_minimum_required(VERSION 3.19)

# Disable "in-source" builds
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(
    FATAL_ERROR
      "ERROR: Will not configure \"in-source\" builds. Create a seperate directory for building"
  )
endif()

# declare the project name
project(
  singularity-eos
  VERSION 1.7.0
  LANGUAGES NONE)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(SINGULARITY_GOLDFILES_VERSION "goldfiles-1.6.2")
set(SINGULARITY_GOLDFILE_HASH
    33ef74b29937cc1347b525f72662933dd2b0556550f6541f97fc1de8a01c3c2a)

# ------------------------------------------------------------------------------#
# Options
# ------------------------------------------------------------------------------#
include(CMakeDependentOption)

# dependcy options
option(SINGULARITY_USE_SPINER "Use spiner EOS" OFF)
cmake_dependent_option(
  SINGULARITY_USE_SPINER_WITH_HDF5 "Use HDF5 component of spiner" ON
  "SINGULARITY_USE_SPINER" OFF)

option(SINGULARITY_USE_FORTRAN "Enable fortran bindings" ON)
option(SINGULARITY_USE_KOKKOS "Use Kokkos for portability" OFF)
option(SINGULARITY_USE_EOSPAC "Enable eospac backend" OFF)

option(SINGULARITY_BUILD_CLOSURE "Mixed Cell Closure" ON)

cmake_dependent_option(SINGULARITY_USE_CUDA "Use CUDA backend of Kokkos" ON
                       "SINGULARITY_USE_KOKKOS" OFF)
cmake_dependent_option(
  SINGULARITY_USE_KOKKOSKERNELS "Use KokkosKernels for LA routines" ON
  "SINGULARITY_USE_KOKKOS;SINGULARITY_BUILD_CLOSURE" OFF)

# extra build options
option(SINGULARITY_BUILD_PYTHON "Compile Python bindings" OFF)
option(SINGULARITY_BUILD_EXAMPLES "Compile examples" OFF)
cmake_dependent_option(
  SINGULARITY_BUILD_SESAME2SPINER
  "Compile sesame2spiner"
  ON
  "SINGULARITY_USE_SPINER;SINGULARITY_USE_SPINER_WITH_HDF5;SINGULARITY_USE_EOSPAC"
  OFF)
cmake_dependent_option(
  SINGUALRITY_BUILD_STELLARCOLLAPSE2SPINER "Compile stellarcollapse2spiner" ON
  "SINGULARITY_USE_SPINER;SINGULARITY_USE_SPINER_WITH_HDF5" OFF)

# testing options
option(SINGULARITY_BUILD_TESTS "Compile tests" OFF)

cmake_dependent_option(
  SINGULARITY_TEST_SESAME "Test Sesame table readers" ON
  "SINGULARITY_BUILD_TESTS;SINGULARITY_BUILD_SESAME2SPINER" OFF)
cmake_dependent_option(
  SINGULARITY_TEST_STELLAR_COLLAPSE "Test stellar collapse table readers" ON
  "SINGULARITY_BUILD_TESTS;SINGUALRITY_BUILD_STELLARCOLLAPSE2SPINER" OFF)
cmake_dependent_option(SINGULARITY_TEST_PYTHON "Test the Python bindings" ON
                       "SINGULARITY_BUILD_TESTS;SINGULARITY_BUILD_PYTHON" OFF)

# modify flags options
option(SINGULARITY_BETTER_DEBUG_FLAGS "Better debug flags for singularity" ON)
option(SINGULARITY_HIDE_MORE_WARNINGS "hide more warnings" OFF)

# toggle code options
option(SINGULARITY_USE_SINGLE_LOGS
       "Use single precision logs. Can harm accuracy." OFF)
option(SINGULARITY_USE_TRUE_LOG_GRIDDING
       "Use grids that conform to log spacing." OFF)

# misc options
option(SINGULARITY_FORCE_SUBMODULE_MODE "Submodule mode" OFF)
# TODO This is an edge-case, but still used (e.g. github CI), but need to work
# out a way to do SUBMODULE_MODE (use submodules for upstream deps) but also do
# a proper install (FORCE_ disables the install/export stages). I may come back
# to this, but there isn't a lot of clamoring for this right now so could return
# to it in a later PR. something like
# cmake_dependent_option(SINGULARITY_SUBMODULE_MODE_INSTALL ...)

# TODO This should depend on using offloading (`SINGULARITY_USE_CUDA` &t) Also,
# patching in the source tree is messy, this should be done in the build dir
option(SINGULARITY_PATCH_MPARK_VARIANT
       "Apply GPU patch to mpark-variant submodule" ON)

# ------------------------------------------------------------------------------#
# singularity-eos Library
# ------------------------------------------------------------------------------#

include(singularity-eos/mpark_variant)
include(singularity-eos/Eigen3)
include(singularity-eos/eospac)
include(singularity-eos/hdf5)
include(singularity-eos/kokkos)
include(singularity-eos/spiner)
include(singularity-eos/ports-of-call)

add_library(singularity-eos)
add_library(singularity-eos::singularity-eos ALIAS singularity-eos)

# ?
target_include_directories(singularity-eos
                           PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>)

# ------------------------------------------------------------------------------#
# Compiler & language setup
# ------------------------------------------------------------------------------#

enable_language(CXX)
include(CMakeDetermineCXXCompiler)
enable_language(C)
include(CMakeDetermineCCompiler)
if(SINGULARITY_USE_FORTRAN)
  enable_language(Fortran)
  include(CMakeDetermineFortranCompiler)
endif()

include(GNUInstallDirs)

if(SINGULARITY_BUILD_PYTHON)
  # need position independent code for Python shared library
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# checks if this is our build, or we've been imported via `add_subdirectory` NB:
# this should make the `option(SINGULARITY_SUBMODULE_MODE ...)` unnecessary
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
else()
  message(
    STATUS
      "Detected that `singularity-eos` is a subproject, will configure build in submodule mode"
  )
  set(SINGULARITY_SUBMODULE_MODE ON)
endif()

if(SINGULARITY_FORCE_SUBMODULE_MODE)
  message(STATUS "Building as though project was a submodule.")
  set(SINGULARITY_SUBMODULE_MODE ON)
endif()

# Use to determine if Eigen is used or not
set(SINGULARITY_USE_EIGEN
    OFF
    CACHE BOOL "" FORCE)
if(SINGULARITY_BUILD_CLOSURE AND NOT SINGULARITY_USE_KOKKOSKERNELS)
  if(SINGULARITY_USE_CUDA)
    message(
      FATAL_ERROR
        "\"SINGULARITY_BUILD_CLOSURE=ON\" and \"SINGULARITY_USE_CUDA=ON\" requires \"SINGULARITY_USE_KOKKOSKERNELS=ON\""
    )
  endif()
  set(SINGULARITY_USE_EIGEN
      ON
      CACHE BOOL "" FORCE)
endif()

# ------------------------------------------------------------------------------#
# De-thaw some options
# ------------------------------------------------------------------------------#

# TODO: MAUNEYC place this in spiner resolution if (NOT SINGULARITY_USE_HDF5)
# message(WARNING "EOSPAC enabled without hdf5.") endif()

# ------------------------------------------------------------------------------#
# Process some options
# ------------------------------------------------------------------------------#

# if building in-tree as a subproject, disable standalone options
if(SINGULARITY_SUBMODULE_MODE)
  set(BUILD_TESTING
      OFF
      CACHE BOOL "" FORCE)
  set(CMAKE_FIND_USE_PACKAGE_REGISTRY
      OFF
      CACHE BOOL "" FORCE)
  set(CMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY
      OFF
      CACHE BOOL "" FORCE)
  set(SINGULARITY_BETTER_DEBUG_FLAGS
      OFF
      CACHE BOOL "" FORCE)
  set(SINGULARITY_HIDE_MORE_WARNINGS
      ON
      CACHE BOOL "" FORCE)
endif()

# side-projects TODO: hdf5 again if(SINGULARITY_BUILD_EXAMPLES)
# add_subdirectory(example) endif()

# add subdirs
if(SINGULARITY_BUILD_PYTHON)
  add_subdirectory(python)
endif()

if(SINGULARITY_BUILD_SESAME2SPINER)
  add_subdirectory(sesame2spiner)
endif()

# defines
if(SINGULARITY_USE_SINGLE_LOGS)
  target_compile_definitions(singularity-eos PUBLIC SINGULARITY_USE_SINGLE_LOGS)
endif()

if(SINGULARITY_TEST_SESAME)
  target_compile_definitions(singularity-eos PRIVATE SINGULARITY_TEST_SESAME)
endif()
if(SINGULARITY_BUILD_CLOSURE)
  target_compile_definitions(singularity-eos PRIVATE SINGULARITY_BUILD_CLOSURE)
endif()

# ------------------------------------------------------------------------------#
# Handle dependencies
# ------------------------------------------------------------------------------#

# There are two modes for building: 1.) Submodule mode this mode is for projects
# where singularity-eos, and it's dependencies, are present in the source tree.
# This mode is only appropriate for projects that have been designed around
# using git submodules for dependency management. Submodule mode disables all
# export/install steps. 2.) Standalone mode this mode will build and install
# singularity-eos as as a complete software package for the intended platform.
# In standalone mode, all dependencies are expected to be found with
# `find_package`, and an error will be produced if the packages required for
# building are not located outside of the source directory, except for explicit
# cases.

if(SINGULARITY_USE_SPINER_WITH_HDF5)
  singularity_enable_hdf5(singularity-eos)
endif()

if(SINGULARITY_USE_EOSPAC)
  # NB This will add the `eospac-wrapper` directory.
  singularity_enable_eospac(singularity-eos)
endif()

if(SINGULARITY_SUBMODULE_MODE)
  # add all submodules
  message(STATUS "singularity-eos configuring in submodule mode.")
  singularity_import_mpark_variant()
  singularity_import_ports_of_call()
  if(SINGULARITY_USE_SPINER)
    singularity_import_spiner()
  endif()
  if(SINGULARITY_USE_KOKKOS)
    singularity_import_kokkos()
    if(SINGULARITY_USE_KOKKOSKERNELS)
      singularity_import_kokkoskernels()
    endif()
  endif()

  if(SINGULARITY_USE_EIGEN)
    singularity_import_eigen()
  endif()
else()
  # use system packages
  singularity_find_mpark_variant()
  singularity_find_ports_of_call()
  if(SINGULARITY_USE_SPINER)
    singularity_find_spiner()
  endif()

  if(SINGULARITY_USE_KOKKOS)
    singularity_find_kokkos()
    if(SINGULARITY_USE_KOKKOSKERNELS)
      singularity_find_kokkoskernels()
    endif()
  endif()
  if(SINGULARITY_USE_EIGEN)
    singularity_find_eigen()
  endif()

endif()

singularity_enable_mpark_variant(singularity-eos)
singularity_enable_ports_of_call(singularity-eos)

if(SINGULARITY_USE_SPINER)
  singularity_enable_spiner(singularity-eos)
  # if(SINGULARITY_USE_SPINER_WITH_HDF5)
  # singularity_enable_hdf5(singularity-eos) endif()
endif()

if(SINGULARITY_USE_KOKKOS)
  singularity_enable_kokkos(singularity-eos)
  if(SINGULARITY_USE_KOKKOSKERNELS)
    singularity_enable_kokkoskernels(singularity-eos)
  endif()
endif()
if(SINGULARITY_USE_EIGEN)
  singularity_enable_eigen(singularity-eos)
endif()

# ----------------------------------------------------------------------------#
# Optional dependencies
# ----------------------------------------------------------------------------#

# need this for consistant Catch2 up/downstream.

if(SINGULARITY_BUILD_TESTS)
  if(NOT TARGET Catch2::Catch2)
    find_package(Catch2 QUIET)
  endif()

  if(NOT Catch2_FOUND)
    message(STATUS "Fetching Catch2 as needed")
    # idiomatic FetchContent
    include(FetchContent)
    FetchContent_Declare(
      Catch2
      GIT_REPOSITORY https://github.com/catchorg/Catch2.git
      # MAUNEYC: On my side, this version works. The ones pinned onto upstream
      # codes do not like to play with clang and/or newer compilers.
      GIT_TAG v2.13.7)
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/contrib)
  endif()
endif()

# ##########################OLD
# ------------------------------------------------------------------------------#
# singularity-eos library
# ------------------------------------------------------------------------------#

# this subdirectory populates `EOS_HEADERS/EOS_SRCS` NOTE: these include path
# prefixes of subdirectories on files (e.g. eos/eos.hpp) see
# singularity-eos/CMakeLists.txt

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

# make sure .mods are placed in build path, and installed along with includes
if(SINGULARITY_USE_FORTRAN)
  set_target_properties(singularity-eos PROPERTIES Fortran_MODULE_DIRECTORY
                                                   ${CMAKE_CURRENT_BINARY_DIR})
  target_include_directories(
    singularity-eos INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                              $<INSTALL_INTERFACE:include/singularity-eos/eos>)
endif() # SINGULARITY_USE_FORTRAN

target_include_directories(
  singularity-eos PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                         $<INSTALL_INTERFACE:include>)

# plug in collected includes/libs/definitions

target_include_directories(
  singularity-eos
  PUBLIC ${SINGULARITY_PUBLIC_INCS}
  PRIVATE ${SINGULARITY_PRIVATE_INCS})

target_link_libraries(
  singularity-eos
  PUBLIC ${SINGULARITY_PUBLIC_LIBS}
  PRIVATE ${SINGULARITY_PRIVATE_LIBS})

target_compile_definitions(
  singularity-eos
  PUBLIC ${SINGULARITY_PUBLIC_DEFINES}
  PRIVATE ${SINGULARITY_PRIVATE_DEFINES})

# ------------------------------------------------------------------------------#
# build flags
# ------------------------------------------------------------------------------#

set(with_cxx "$<COMPILE_LANGUAGE:CXX>")
set(with_xlcxx "$<COMPILE_LANG_AND_ID:CXX,XL>")
set(build_debug "$<CONFIG:Debug>")
set(build_release "$<CONFIG:Release>")
set(with_cuda "$<BOOL:${SINGULARITY_USE_CUDA}>")
set(with_warnings "$<NOT:$<BOOL:${SINGULARITY_HIDE_MORE_WARNINGS}>>")
set(with_better_flags "$<BOOL:${SINGULARITY_BETTER_DEBUG_FLAGS}>")

set(xlfix "$<${with_xlcxx}:-std=c++1y;-qxflag=disable__cplusplusOverride>")

target_compile_options(
  singularity-eos
  PRIVATE ${xlfix}
  PUBLIC $<${with_cuda}:
         $<${with_cxx}:
         --expt-relaxed-constexpr
         $<${with_warnings}:
         -Xcudafe;--diag_suppress=esa_on_defaulted_function_ignored
         > # with_warnings
         $<${build_debug}:
         $<${with_better_flags}:
         -G;-lineinfo
         > # better_flags
         > # debug
         > # cxx
         $<${build_release}:
         -use_fast_math
         > # release
         > # cuda
  PUBLIC -fmax-errors=3)

target_link_options(singularity-eos PRIVATE ${xlfix})

# ----------------------------------------------------------------------------#
# Export/install
# ----------------------------------------------------------------------------#

# if in submod mode, assume `add_subdirectory` has been done and do not generate
# export files or install anything
if(NOT SINGULARITY_SUBMODULE_MODE)
  include(cmake/install.cmake)
endif()

# ----------------------------------------------------------------------------#
# Testing
# ----------------------------------------------------------------------------#

if(SINGULARITY_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

# ----------------------------------------------------------------------------#
# Misc.
# ----------------------------------------------------------------------------#

# clang format
include(cmake/Format.cmake)
