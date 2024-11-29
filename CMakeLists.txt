# ------------------------------------------------------------------------------#
# © 2021-2024. Triad National Security, LLC. All rights reserved.  This program
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
  VERSION 1.9.1
  LANGUAGES NONE)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set(SINGULARITY_GOLDFILES_VERSION "goldfiles-1.8.0")
set(SINGULARITY_GOLDFILE_HASH
    249772f3314c4b6b9386aa08895280698ae905ef188f4383b819f28c1484603b)
set(DOWNLOAD_EXTRACT_TIMESTAMP ON)

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
# Optionally build these if you want to build/test the fortran
# infrastructure without invoking a fortran compiler If fortran is
# off, you can set this. If fortran is on, it's forced to ON and is
# not set-able.
cmake_dependent_option(SINGULARITY_BUILD_FORTRAN_BACKEND
  "Build the C++ code to which the fortran bindings bind"
  OFF "NOT SINGULARITY_USE_FORTRAN" ON)

option(SINGULARITY_USE_KOKKOS "Use Kokkos for portability" OFF)
option(SINGULARITY_USE_EOSPAC "Enable eospac backend" OFF)
option(SINGULARITY_EOSPAC_ENABLE_SHMEM
  "Support shared memory with EOSPAC backend. Requires EOSPAC version 6.5.7."
  OFF)

# TODO This should be dependent option (or fortran option)
option(SINGULARITY_BUILD_CLOSURE "Mixed Cell Closure" ON)

cmake_dependent_option(SINGULARITY_USE_CUDA "Use CUDA backend of Kokkos" OFF
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
  SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER "Compile stellarcollapse2spiner" ON
  "SINGULARITY_USE_SPINER;SINGULARITY_USE_SPINER_WITH_HDF5" OFF)
cmake_dependent_option(
  SINGULARITY_USE_HELMHOLTZ "Include Helmholtz equation of state" OFF
  "SINGULARITY_USE_SPINER;SINGULARITY_USE_SPINER_WITH_HDF5" OFF)

# testing options
option(SINGULARITY_BUILD_TESTS "Compile tests" OFF)

cmake_dependent_option(
  SINGULARITY_TEST_SESAME "Test Sesame table readers" ON
  "SINGULARITY_BUILD_TESTS;SINGULARITY_BUILD_SESAME2SPINER" OFF)
cmake_dependent_option(
  SINGULARITY_TEST_STELLAR_COLLAPSE "Test stellar collapse table readers" ON
  "SINGULARITY_BUILD_TESTS;SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER" OFF)
cmake_dependent_option(SINGULARITY_TEST_PYTHON "Test the Python bindings" ON
                       "SINGULARITY_BUILD_TESTS;SINGULARITY_BUILD_PYTHON" OFF)
cmake_dependent_option(
  SINGULARITY_TEST_HELMHOLTZ "Test the Helmholtz equation of state" ON
  "SINGULARITY_BUILD_TESTS;SINGULARITY_USE_HELMHOLTZ" OFF)

# modify flags options
option(SINGULARITY_BETTER_DEBUG_FLAGS "Better debug flags for singularity" ON)
option(SINGULARITY_HIDE_MORE_WARNINGS "hide more warnings" OFF)
option(SINGULARITY_STRICT_WARNINGS "Make warnings strict" OFF)

# toggle code options
option(SINGULARITY_USE_TRUE_LOG_GRIDDING
  "Use grids that conform to log spacing." OFF)
cmake_dependent_option(SINGULARITY_USE_SINGLE_LOGS
  "Use single precision logs. Only available for true log gridding. Can harm accuracy."
  OFF "SINGULARITY_USE_TRUE_LOG_GRIDDING" OFF)
option(SINGULARITY_NQT_ORDER_1
  "In NQT logs, use first order. Faster but less accurate."
  OFF)
option(SINGULARITY_NQT_PORTABLE
  "In NQT logs, use portable, rather than bithacked implementation. Slower, but more likely to function on exotic architectures."
  OFF)

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

# Plugins
set(SINGULARITY_PLUGINS "" CACHE STRING "List of paths to plugin directories")
set(SINGULARITY_VARIANT "singularity-eos/eos/default_variant.hpp" CACHE STRING
  "The include path for the file containing the definition of singularity::EOS.")

# Detect ARM architecture
set(SINGULARITY_ON_ARM OFF CACHE BOOL "We are running on an ARM system")
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm|aarch64)")
    if (NOT SINGULARITY_USE_CUDA)
      message(STATUS
        "ARM architecture detected: ${CMAKE_SYSTEM_PROCESSOR}")
      set(SINGULARITY_ON_ARM ON CACHE BOOL
        "We are running on an ARM system")
    endif()
endif()

if (SINGULAIRTY_ON_ARM)
  if (NOT SINGULARITY_USE_TRUE_LOG_GRIDDING)
    message(WARNING
      "Fast logs not necessarily better on ARM CPU systems. "
      "You may wish to build with "
      "-DSINGULARITY_USE_TRUE_LOG_GRIDDING=ON.")
  endif()
endif()

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

add_library(singularity-eos INTERFACE)
add_library(singularity-eos::singularity-eos ALIAS singularity-eos)

# interface target to collect runtime libraries
add_library(singularity-eos_Common INTERFACE)
add_library(singularity-eos::singularity-eos_Common ALIAS singularity-eos_Common)

add_library(singularity-eos_Interface INTERFACE)
add_library(singularity-eos::singularity-eos_Interface ALIAS singularity-eos_Interface)
target_link_libraries(singularity-eos_Interface INTERFACE singularity-eos_Common)
target_link_libraries(singularity-eos INTERFACE singularity-eos_Interface)

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

# Big endianness
include(TestBigEndian)
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)
if (BIG_ENDIAN)
  message(WARNING "Big endian detected! "
    "Integer aliasing as currently implemented will not function. "
    "Please set -DSINGULARITY_NQT_PORTABLE=ON.")
endif()

include(GNUInstallDirs)

if(SINGULARITY_BUILD_PYTHON)
  # need position independent code for Python shared library
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# require at least C++14, but allow newer versions to become a client requirement if
# explicitly set at build time (needed for building with newer Kokkos)
if(CMAKE_CXX_STANDARD)
  target_compile_features(singularity-eos_Interface INTERFACE cxx_std_${CMAKE_CXX_STANDARD})
else()
  target_compile_features(singularity-eos_Interface INTERFACE cxx_std_17)
endif()

# checks if this is our build, or we've been imported via `add_subdirectory` NB:
# this should make the `option(SINGULARITY_SUBMODULE_MODE ...)` unnecessary
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
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

if(SINGULARITY_BUILD_EXAMPLES)
  add_subdirectory(example)
endif()

# add subdirs
if(SINGULARITY_BUILD_PYTHON)
  add_subdirectory(python)
endif()

if(SINGULARITY_BUILD_SESAME2SPINER)
  add_subdirectory(sesame2spiner)
endif()

# Define the full version as a string macro
target_compile_definitions(singularity-eos_Interface INTERFACE
    SINGULARITY_VERSION=\"${PROJECT_VERSION}\"
)

# Optionally, define major, minor, and patch versions separately
target_compile_definitions(singularity-eos_Interface INTERFACE
    SINGULARITY_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
    SINGULARITY_VERSION_MINOR=${PROJECT_VERSION_MINOR}
    SINGULARITY_VERSION_PATCH=${PROJECT_VERSION_PATCH}
)

# defines
if (SINGULARITY_USE_TRUE_LOG_GRIDDING)
  target_compile_definitions(singularity-eos_Interface
    INTERFACE SINGULARITY_USE_TRUE_LOG_GRIDDING)
endif()
if(SINGULARITY_USE_SINGLE_LOGS)
  target_compile_definitions(singularity-eos_Interface
    INTERFACE SINGULARITY_USE_SINGLE_LOGS)
endif()
if(SINGULARITY_NQT_ORDER_1)
  target_compile_definitions(singularity-eos_Interface
    INTERFACE SINGULARITY_NQT_ORDER_1)
endif()
if(SINGULARITY_NQT_PORTABLE)
  target_compile_definitions(singularity-eos_Interface
    INTERFACE SINGULARITY_NQT_PORTABLE)
endif()
if(SINGULARITY_TEST_SESAME)
  target_compile_definitions(singularity-eos_Interface INTERFACE SINGULARITY_TEST_SESAME)
endif()
if(SINGULARITY_BUILD_CLOSURE)
  target_compile_definitions(singularity-eos_Interface INTERFACE SINGULARITY_BUILD_CLOSURE)
endif()
if(SINGULARITY_USE_HELMHOLTZ)
  target_compile_definitions(singularity-eos_Interface INTERFACE SINGULARITY_USE_HELMHOLTZ)
endif()
if(SINGULARITY_USE_SPINER_WITH_HDF5)
  target_compile_definitions(singularity-eos_Interface INTERFACE SINGULARITY_USE_SPINER_WITH_HDF5)
endif()
if (SINGULARITY_USE_EOSPAC AND SINGULARITY_EOSPAC_ENABLE_SHMEM)
  target_compile_definitions(singularity-eos_Interface INTERFACE SINGULARITY_EOSPAC_ENABLE_SHARED_MEMORY)
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
  singularity_enable_hdf5(singularity-eos_Common)
endif()

if(SINGULARITY_USE_EOSPAC)
  # NB This will add the `eospac-wrapper` directory.
  singularity_enable_eospac(singularity-eos_Common)
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

singularity_enable_mpark_variant(singularity-eos_Interface)
singularity_enable_ports_of_call(singularity-eos_Interface)

if(SINGULARITY_USE_SPINER)
  singularity_enable_spiner(singularity-eos_Interface)
  # if(SINGULARITY_USE_SPINER_WITH_HDF5)
  # singularity_enable_hdf5(singularity-eos_Interface) endif()
endif()

# Both the interface (headers) and the library (compiled) need to link to Kokkos
# see get_sg_eos.cpp
if(SINGULARITY_USE_KOKKOS)
  singularity_enable_kokkos(singularity-eos_Common)
  if(SINGULARITY_USE_KOKKOSKERNELS)
    singularity_enable_kokkoskernels(singularity-eos_Common)
  endif()
endif()
if(SINGULARITY_USE_EIGEN)
  singularity_enable_eigen(singularity-eos_Interface)
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
      # or later is fine too
      GIT_TAG v3.0.1)
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/contrib)
  endif()
endif()

# ------------------------------------------------------------------------------#
# Plugin infrastructure
# ------------------------------------------------------------------------------#
include(cmake/plugins.cmake)
set(PLUGIN_BUILD_ROOT "${CMAKE_CURRENT_BINARY_DIR}/plugins")

# ------------------------------------------------------------------------------#
# singularity-eos library
# ------------------------------------------------------------------------------#

# Here we populate `EOS_HEADERS/EOS_SRCS` NOTE: these include path
# prefixes of subdirectories on files (e.g. eos/eos.hpp) see
# singularity-eos/CMakeLists.txt
add_subdirectory(singularity-eos)

foreach(_plugin ${SINGULARITY_PLUGINS})
  message(STATUS "Adding plugin ${_plugin}...")
  get_filename_component(BINDIR_LOC ${_plugin} NAME)
  add_subdirectory(${_plugin} ${PLUGIN_BUILD_ROOT}/${BINDIR_LOC})
endforeach()

# TODO(JMM): Kind of nice to have?
get_property(eos_headers GLOBAL PROPERTY EOS_HEADERS)
target_sources(singularity-eos_Interface PRIVATE ${eos_headers})
message(VERBOSE "EOS Headers:\n\t${eos_headers}")

if(SINGULARITY_BUILD_CLOSURE)
  get_property(eos_srcs GLOBAL PROPERTY EOS_SRCS)

  if(eos_srcs)
    add_library(singularity-eos_Library)
    set_target_properties(singularity-eos_Library PROPERTIES OUTPUT_NAME singularity-eos)
    add_library(singularity-eos::singularity-eos_Library ALIAS singularity-eos_Library)
    # Public scope ensures explicit Kokkos dependency for library and anything
    # that links to it
    target_link_libraries(singularity-eos_Library PUBLIC singularity-eos_Common)
    target_link_libraries(singularity-eos INTERFACE singularity-eos_Library)
    target_sources(singularity-eos_Library PRIVATE ${eos_srcs})
    message(VERBOSE "EOS Library Sources:\n\t${eos_srcs}")
  endif()
endif()

if(TARGET singularity-eos_Library AND SINGULARITY_USE_FORTRAN)
  # Turn on preprocessor for fortran files
  set_target_properties(singularity-eos_Library PROPERTIES Fortran_PREPROCESS ON)
  # make sure .mods are placed in build path, and installed along with includes
  set_target_properties(singularity-eos_Library PROPERTIES Fortran_MODULE_DIRECTORY
                                                   ${CMAKE_CURRENT_BINARY_DIR}/fortran)
  target_include_directories(
    singularity-eos_Library INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran>
                              $<INSTALL_INTERFACE:include/singularity-eos/eos>)
endif() # SINGULARITY_USE_FORTRAN

target_include_directories(singularity-eos_Interface INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/generated>
  $<INSTALL_INTERFACE:include>
)

get_property(plugin_include_paths GLOBAL PROPERTY PLUGIN_INCLUDE_PATHS)
foreach(path ${plugin_include_paths})
  target_include_directories(singularity-eos_Interface INTERFACE $<BUILD_INTERFACE:${path}>)
endforeach()

# plug in collected includes/libs/definitions

target_include_directories(
  singularity-eos_Interface
  PUBLIC ${SINGULARITY_PUBLIC_INCS}
  PRIVATE ${SINGULARITY_PRIVATE_INCS})

target_link_libraries(
  singularity-eos_Interface
  PUBLIC ${SINGULARITY_PUBLIC_LIBS}
  PRIVATE ${SINGULARITY_PRIVATE_LIBS})

target_compile_definitions(
  singularity-eos_Interface
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
  singularity-eos_Interface
  INTERFACE $<${with_cuda}:
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
)
if (SINGULARITY_STRICT_WARNINGS)
  target_compile_options(singularity-eos_Interface INTERFACE
    -Wall -Werror)
endif()

if(TARGET singularity-eos_Library)
  target_compile_options(singularity-eos_Library PRIVATE ${xlfix})
  target_link_options(singularity-eos_Library PRIVATE ${xlfix})

  # target_link_libraries(singularity-eos_Library PRIVATE singularity-eos_Interface)
  #
  # Can not use PRIVATE, since it would add singularity-eos_Interface as LINK_ONLY
  # in the singularity-eos_Library export
  #
  # CMake 3.26 gives us $<BUILD_LOCAL_INTERFACE:...>, for now we need to do this workaround
  if(CMAKE_VERSION VERSION_LESS "3.26.0")
    target_link_options(singularity-eos_Library PRIVATE $<TARGET_PROPERTY:singularity-eos_Interface,INTERFACE_LINK_OPTIONS>)
    target_include_directories(singularity-eos_Library PRIVATE $<TARGET_PROPERTY:singularity-eos_Interface,INTERFACE_INCLUDE_DIRECTORIES>)
    target_compile_options(singularity-eos_Library PRIVATE $<TARGET_PROPERTY:singularity-eos_Interface,INTERFACE_COMPILE_OPTIONS>)
    target_compile_definitions(singularity-eos_Library PRIVATE $<TARGET_PROPERTY:singularity-eos_Interface,INTERFACE_COMPILE_DEFINITIONS>)
    target_sources(singularity-eos_Library PRIVATE $<TARGET_PROPERTY:singularity-eos_Interface,INTERFACE_SOURCES>)
  else()
    target_link_libraries(singularity-eos_Library PRIVATE $<BUILD_LOCAL_INTERFACE:singularity-eos_Interface>)
  endif()
endif()

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
  include(CTest)
  enable_testing()
  add_subdirectory(test)
endif()

# ----------------------------------------------------------------------------#
# Misc.
# ----------------------------------------------------------------------------#

# clang format
include(cmake/Format.cmake)
