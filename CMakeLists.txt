#------------------------------------------------------------------------------#
# © 2021-2022. Triad National Security, LLC. All rights reserved.  This
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

cmake_minimum_required(VERSION 3.23)

project(singularity-eos
  VERSION 1.7.1 
  LANGUAGES NONE)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(SINGULARITY_GOLDFILES_VERSION "goldfiles-1.6.2")
set(SINGULARITY_GOLDFILE_HASH
  33ef74b29937cc1347b525f72662933dd2b0556550f6541f97fc1de8a01c3c2a)

#------------------------------------------------------------------------------#
# Use better messages
#------------------------------------------------------------------------------#
include(cmake/singularity_utils.cmake)

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
option (SINGULARITY_SUBMODULE_MODE "Submodule mode" OFF)
option (SINGULARITY_BUILD_CLOSURE "Mixed cell closure" ON)
option (SINGULARITY_TEST_SESAME "Test the Sesame table readers" OFF)
option (SINGULARITY_TEST_STELLAR_COLLAPSE "Test the stellar collapse table readers" OFF)
option (SINGULARITY_TEST_PYTHON "Test the Python bindings" OFF)
option (SINGULARITY_USE_SINGLE_LOGS "Use single precision logs. Can harm accuracy." OFF)
option (SINGULARITY_USE_TRUE_LOG_GRIDDING "Use grids that conform to log spacing." OFF)
option (SINGULARITY_PATCH_MPARK_VARIANT "Apply GPU patch to mpark-variant submodule" ON)

# keep config data stored in easily accessable vars,
# only attaching them to targets as late as possible.
# this allows for easier handling and manipulating, as
# well as later info printing and debugging
set(SINGULARITY_PUBLIC_DEFINES "")    # list to `target_compile_definitions()`
set(SINGULARITY_PRIVATE_DEFINES "")   # """
set(SINGULARITY_PUBLIC_INCS "")       # list to `target_include_directories()`
set(SINGULARITY_PRIVATE_INCS "")      # """ 
set(SINGULARITY_PUBLIC_LIBS "")       # list to `target_link_libraries()`
set(SINGULARITY_PRIVATE_LIBS "")      # """
set(SINGULARITY_EXPORT_TARGETS "")    # probably not useful
set(SINGULARITY_COMPILEFLAGS "")      # TODO: actually use these
set(SINGULARITY_DEP_PKGS "")          # a list of packages, used to generate import rules for installed targets. (may want to rename)

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

# Don't allow in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR
    "You cannot build in a source directory (or any directory with a CMakeLists.txt file). "
    "Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
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


if (SINGULARITY_BUILD_PYTHON)
  list(APPEND _subdirs2add python)
endif()

if(SINGULARITY_USE_SINGLE_LOGS)
  list(APPEND TAR_DEFINES SINGULARITY_USE_SINGLE_LOGS)
endif()

#------------------------------------------------------------------------------#
# Handle dependencies
#------------------------------------------------------------------------------#

# need this for consistant Catch2 up/downstream.
# depending on how this develops, this could be wrapped
# into the dependcy logic (below)
include(FetchContent)

# If we are on version 3.24+, then set FetchContent to always try 
# to `find_package` before trying a download method
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE ALWAYS)
  message(STATUS "FetchContent routines will try `find_package` first")
else()
  message(DEPRECATION 
    "Detected cmake version (${CMAKE_VERSION}) older then 3.24. `spiner`"
    "will begin requiring this version soon.  The current depedency "
    "resolution is a soft-copy of the pattern of FetchContent introduced" 
    "in cmake 3.24."
  )
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

include(content)
# PoC
singularityeos_content_declare(
  ports-of-call
  GIT_REPO    https://github.com/lanl/ports-of-call
  GIT_TAG     5fea4f6a0ad75ca2326e1c32507e24a3e09e2f8b
  NAMESPACE singularityeosDeps
)

# spiner
if(SINGULARITY_USE_HDF5)
  set(_singularity_spiner_opt "SPINER_USE_HDF")
endif()
singularityeos_content_declare(
  spiner
  EXPECTED_TARGETS spiner::spiner
  GIT_REPO    https://github.com/lanl/spiner
  GIT_TAG     9f3ccf017083af9111f4ea7e2dfb887d2f3579cb
  NAMESPACE singularityeosDeps
  ENABLE_OPTS ${_singularity_spiner_opt}
  DISABLE_OPTS "BUILD_TESTING"
)
# variant
singularityeos_content_declare(
  mpark_variant
  GIT_REPO https://github.com/mpark/variant
  GIT_TAG 4988879a9f5a95d72308eca2b1779db6ed9b135d
  PATCH ${CMAKE_CURRENT_SOURCE_DIR}/utils/cuda_compatibility.patch
  EXPECTED_TARGETS mpark_variant
  NAMESPACE singularityeosDeps
)

if(SINGULARITY_USE_KOKKOS)
  # kokkos offloading
  # put together a list of kokkos configure options
  list(APPEND _singularity_kokkos_opts "Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION")
  if (singularity_kokkos_opts_TEST_USE_KOKKOS_CUDA)
    list(APPEND _singularity_kokkos_opts "Kokkos_ENABLE_CUDA")
    list(APPEND _singularity_kokkos_opts "Kokkos_ENABLE_CUDA_LAMBDA")
    list(APPEND _singularity_kokkos_opts "Kokkos_ENABLE_CUDA_CONSTEXPR")
  endif()
  if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    list(APPEND _singularity_kokkos_opts "Kokkos_ENABLE_DEBUG")
    list(APPEND _singularity_kokkos_opts "Kokkos_ENABLE_DEBUG_BOUNDS_CHECK")
    list(APPEND _singularity_kokkos_opts "Kokkos_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK")
  endif()
 
  singularityeos_content_declare(
    Kokkos
    GIT_REPO https://github.com/kokkos/kokkos
    GIT_TAG        3.7.01
    EXPECTED_TARGETS Kokkos::kokkos
    NAMESPACE singularityeosDeps
    ENABLE_OPTS ${_singularity_kokkos_opts}
   )
  list(APPEND SINGULARITY_PUBLIC_DEFINES PORTABILITY_STRATEGY_KOKKOS)
  ####
  # TODO: Once it's oppurtune, we will
  #       not support Kokkos < 3.3
  if((${Kokkos_VERSION_MAJOR} MATCHES "^[0-9]+$") AND (${Kokkos_VERSION_MINOR} MATCHES "^[0-9]+$"))
    if(${Kokkos_VERSION_MAJOR}.${Kokkos_VERSION_MINOR} VERSION_LESS "3.3")
      message(WARNING "`Kokkos` version [${Kokkos_VERSION_MAJOR}.${Kokkos_VERSION_MINOR}] is DEPRECATED, and `singularity-eos` will not support versions < '3.3' in the very near future.")
      if(SINGULARITY_USE_CUDA)
        get_filename_component(_compiler_exe ${CMAKE_CXX_COMPILER} NAME)
        if(NOT _compiler_exe STREQUAL "nvcc_wrapper")
          message(FATAL_ERROR "To use CUDA offloading with `Kokkos` version <= 3.2, you must use `-DCMAKE_CXX_COMPILER=/path/to/kokkos/bin/nvcc_wrapper` when invoking `cmake`. You may also need to set a Kokkos target architecture via: https://github.com/kokkos/kokkos/wiki/Compiling.")
        endif()
     endif() #SINGULARITY_USE_CUDA
    endif() #Kokkos_VERSION_MAJOR
  else()
    message(WARNING "`Kokkos` version could not be inferred. Note that versions less than `3.3` are DEPRECATED and may not work properly.")
  endif() #Kokkos_VERSION_MAJOR
  ####
endif() #SINGULARITY_USE_KOKKOS

# TODO(JMM): If linear algebra needs grow beyond mixed cell closures,
# we may need to rename this variable.
if (SINGULARITY_BUILD_CLOSURE)
  if(SINGULARITY_USE_KOKKOSKERNELS)
    # kokkoskernels solvers
    # Disable TPLs
    list(APPEND _kk_tpls "BLAS" "MKL" "LAPACK" "CUBLAS" "CUPARSE" "MAGMA")
    foreach(_tpl ${_kk_tpls})
      list(APPEND _singularity_kokkoskernels_opt_off "KokkosKernels_ENABLE_TPL_${_tpl}")
    endforeach()
    # Disable ETIs
    list(APPEND _kk_etis "COMPLEX_DOUBLE" "COMPLEX_FLOAT" "DOUBLE" "EXECSPACE_OPENMP" "EXECSPACE_SERIAL")
    list(APPEND _kk_etis "EXECSPACE_THREADS" "EXECSPACE_CUDA" "FLOAT" "LAYOUTLEFT" "LAYOUTRIGHT")
    list(APPEND _kk_etis "MEMSPACE_HOSTSPACE" "OFFSET_INT" "OFFSET_SIZE_T" "ORDINAL_INT")
    foreach(_eti ${_kk_etis})
      list(APPEND _singularity_kokkoskernels_opt_off "KokkosKernels_INST_${_eti}")
    endforeach()

    singularityeos_content_declare(
      KokkosKernels
      GIT_REPO https://github.com/kokkos/kokkos-kernels
      GIT_TAG        3.7.01
      EXPECTED_TARGETS Kokkos::kokkoskernels
      NAMESPACE singularityeosDeps
      DISABLE_OPTS ${_singularity_kokkoskernels_opt_off}
    )
    list(APPEND SINGULARITY_PUBLIC_DEFINES SINGULARITY_USE_KOKKOSKERNELS)
  else()
    list(APPEND _singularity_eigen_opt_off "EIGEN_TEST_CXX11" "EIGEN_BUILD_PKGCONFIG" "EIGEN_BUILD_DOC")
    singularityeos_content_declare(
      Eigen3
      GIT_REPO https://gitlab.com/libeigen/eigen
      GIT_TAG  e7248b26a1ed53fa030c5c459f7ea095dfd276ac
      EXPECTED_TARGETS Eigen3::Eigen
      NAMESPACE singularityeosDeps
      DISABLE_OPTS ${_singularity_eigen_opt_off}
    )
  endif() #SINGULARITY_USE_KOKKOS_KERNELS
endif() #SINGULARITY_BUILD_CLOSURE


# only use hdf5 if spiner

if (SINGULARITY_USE_HDF5)

  list(APPEND SINGULARITY_PUBLIC_DEFINES SPINER_USE_HDF)
  list(APPEND SINGULARITY_PUBLIC_DEFINES SINGULARITY_USE_HDF)
endif() # SINGULARITY_USE_HDF5

# side-projects
if(SINGULARITY_USE_HDF5)
  if(SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER)
    list(APPEND _subdirs2add stellarcollapse2spiner)
  endif()
  if(SINGULARITY_USE_EOSPAC AND SINGULARITY_BUILD_SESAME2SPINER)
    list(APPEND _subdirs2add sesame2spiner)
  endif()
  if(SINGULARITY_BUILD_EXAMPLES)
    list(APPEND _subdirs2add example)
  endif()
endif()

#TODO: Check: Does EOSPAC only pass into side-targets?
if (SINGULARITY_USE_EOSPAC)
  #TODO: should be strictly module mode
  find_package(EOSPAC REQUIRED)
  list(APPEND SINGULARITY_PUBLIC_DEFINES SINGULARITY_USE_EOSPAC)
  list(APPEND SINGULARITY_PUBLIC_LIBS EOSPAC::eospac)
  list(APPEND SINGULARITY_PUBLIC_LIBS eospac-wrapper)
  list(APPEND SINGULARITY_DEPS_PKGS EOSPAC)
  list(APPEND _subdirs2add eospac-wrapper)
endif ()


#TODO: this feels out-of-place here, though i haven't crawled into the source yet to see
# what the point is.
if (SINGULARITY_TEST_SESAME)
  list(APPEND SINGULARITY_PRIVATE_DEFINES SINGULARITY_TEST_SESAME)
endif()
if (SINGULARITY_BUILD_CLOSURE)
  list(APPEND SINGULARITY_PRIVATE_DEFINES SINGULARITY_BUILD_CLOSURE)
endif()

singularityeos_content_populate(
  NAMESPACE singularityeosDeps
)

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

# define the 
add_library(singularity-eos ${_srcs} ${_headers})
add_library(singularity-eos::singularity-eos ALIAS singularity-eos)

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
    ${singularityeosDeps_POPULATED_TARGETS}
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

foreach(_sdir ${_subdirs2add})
  singularity_msg(STATUS "Adding ${_sdir}...")
  add_subdirectory(${_sdir})
endforeach()


# TODO either expand this or do feature_summary
#singularity_msg(STATUS "Configuration")
#message(STATUS "[== INFO: status so far ===============================")
#message(STATUS "[== dependency list: ${SINGULARITY_DEP_PKGS}")
#message(STATUS "[== library targets: ${SINGULARITY_PUBLIC_LIBS} ${SINGULARITY_PRIVATE_LIBS}")
#message(STATUS "[== compiler definitions: ${SINGULARITY_PUBLIC_DEFINES} ${SINGULARITY_PRIVATE_DEFINES}")
#message(STATUS "[======================================================")


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


