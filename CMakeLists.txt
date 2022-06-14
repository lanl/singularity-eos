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

cmake_minimum_required(VERSION 3.14)

project(singularity-eos VERSION 1.5.0)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

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
option (SINGULARITY_USE_TRUE_LOG_GRIDDING "Use grids that conform to log spacing." OFF)


# Patches variant to be compatible with cuda
# Assumes "patch" is present on system
message(STATUS "Patching mpark::variant to support GPUs")
execute_process(COMMAND patch -N -s -V never
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/variant/include/mpark/variant.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/cuda_compatibility.patch
)

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

enable_language(CXX)
include(CMakeDetermineCXXCompiler)
if(SINGULARITY_USE_FORTRAN)
  enable_language(Fortran)
  include(CMakeDetermineFortranCompiler)
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
  set(CMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY OFF CACHE BOOL "" FORCE)
  set(SINGULARITY_BETTER_DEBUG_FLAGS OFF CACHE BOOL "" FORCE)
  set(SINGULARITY_HIDE_MORE_WARNINGS ON CACHE BOOL "" FORCE)
endif()

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

if(SINGULARITY_USE_SINGLE_LOGS)
  list(APPEND TAR_DEFINES SINGULARITY_USE_SINGLE_LOGS)
endif()

#------------------------------------------------------------------------------#
# Handle dependencies
#------------------------------------------------------------------------------#

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

include(cmake/dependency.cmake)

if(SINGULARITY_USE_KOKKOS)
  # kokkos offloading
  # type: soft submodule 
  singularity_import_dependency(
    PKG Kokkos
    TARGET Kokkos::kokkos
    SUBDIR ${PROJECT_SOURCE_DIR}/utils/kokkos
    )
  set(PORTABILITY_STRATEGY "Kokkos")
  set(SPINER_USE_KOKKOS ON)
  list(APPEND SINGULARITY_PUBLIC_DEFINES PORTABILITY_STRATEGY_KOKKOS)
  list(APPEND SINGULARITY_PUBLIC_DEFINES SPINER_USE_KOKKOS)
  ####
  # TODO: Once it's oppurtune, we will
  #       not support Kokkos < 3.3
  if(Kokkos_VERSION_MAJOR VERSION_LESS "3.3")
    message(WARNING "`Kokkos` version [${Kokkos_VERSION_MAJOR}] is DEPRECATED, and `singularity-eos` will not support versions < '3.3' in the very near future.")
    message(WARNING "You may need to set `-DCMAKE_CXX_COMPILER=/path/to/nvcc_wrapper` by hand. You may also need to set a Kokkos target architecture via: https://github.com/kokkos/kokkos/wiki/Compiling")
    if(SINGULARITY_USE_CUDA)
      if(NOT "$ENV{CXX}x" STREQUAL "x" AND
          "$ENV{NVCC_WRAPPER_DEFAULT_COMPILER}x" STREQUAL "x")
        set(ENV{NVCC_WRAPPER_DEFAULT_COMPILER} "$ENV{CXX}")
      endif()
      message(STATUS "Looking for nvcc_wrapper in ${Kokkos_DIR} and ${Kokkos_SOURCE_DIR}")
      find_program(SINGULARITY_NVCC_WRAPPER
        NAMES           nvcc_wrapper
        HINTS           "${Kokkos_DIR}" "${Kokkos_SOURCE_DIR}"
        PATHS           "${Kokkos_DIR}" "${Kokkos_SOURCE_DIR}"
        PATH_SUFFIXES   bin)
      # fatal if we can't nvcc_wrapper
      if(NOT SINGULARITY_NVCC_WRAPPER)
        message(FATAL_ERROR "Could not find nvcc_wrapper.")
      endif()
      set(CMAKE_CXX_COMPILER ${SINGULARITY_NVCC_WRAPPER} CACHE STRING "")
    endif() #SINGULARITY_USE_CUDA
  endif() #Kokkos_VERSION_MAJOR
  ####
endif() #SINGULARITY_USE_KOKKOS

# TODO(JMM): If linear algebra needs grow beyond mixed cell closures,
# we may need to rename this variable.
if (SINGULARITY_BUILD_CLOSURE)
  if(SINGULARITY_USE_KOKKOSKERNELS)
    # type: external
    singularity_import_dependency(
      PKG KokkosKernels
      TARGET Kokkos::kokkoskernels
      )
    list(APPEND SINGULARITY_PUBLIC_DEFINES SINGULARITY_USE_KOKKOSKERNELS)
  else()
    # type: soft submodule
    singularity_import_dependency(
      PKG Eigen3
      TARGET Eigen3::Eigen
      SUBDIR ${PROJECT_SOURCE_DIR}/utils/eigen
      )
  endif() #SINGULARITY_USE_KOKKOS_KERNELS
endif() #SINGULARITY_BUILD_CLOSURE

# variant for prior C++17
# type: hard submodule
singularity_import_dependency(
  SUBMODULE_ONLY
  PKG mpark_variant
  TARGET mpark_variant
  SUBDIR ${PROJECT_SOURCE_DIR}/utils/variant
)

# ports-of-call
# type: ghost submodule
singularity_import_dependency(
  SUBMODULE_ONLY
  PKG ports-of-call
  TARGET ports-of-call::ports-of-call
  SUBDIR ${PROJECT_SOURCE_DIR}/utils/ports-of-call
)

# spiner
# type: soft? submodule
singularity_import_dependency(
  PKG spiner
  TARGET spiner::spiner
  SUBDIR ${PROJECT_SOURCE_DIR}/utils/spiner
)

# only use hdf5 if spiner
if (SINGULARITY_USE_HDF5)
  # hdf5, for data interopability
  # type: external
  singularity_import_dependency(
    PKG HDF5
    COMPONENTS C HL
    #NB: 'newer' CMake versions, HDF5s will
    #    do the `idiomatic` cmake setup,
    #    and `find_package(HDF5)` will define
    #    these targets.
    #    however, some downstream codes
    #    are still struggling through the
    #    dark times, so we use the 
    #    'standard' hack for now
    #    TARGETS HDF5::HDF5 hdf5::hdf5_hl
  )

  # make sure we use MPI,
  # as this requirement is not
  # transitive in hdf5
  if(HDF5_IS_PARALLEL)
    # Message passing interface
    # type: external
    singularity_import_dependency(
      PKG MPI
      COMPONENTS CXX
      TARGETS MPI::MPI_CXX
    )
  endif()
  ####
  #TODO: when downstream is rescued,
  #      we can discard these
  list(APPEND SINGULARITY_PUBLIC_LIBS ${HDF5_LIBRARIES})
  list(APPEND SINGULARITY_PUBLIC_LIBS ${HDF5_HL_LIBRARIES})
  list(APPEND SINGULARITY_PUBLIC_INCS ${HDF5_INCLUDE_DIRS})
  ####
  list(APPEND SINGULARITY_PUBLIC_DEFINES SPINER_USE_HDF)
  list(APPEND SINGULARITY_PUBLIC_DEFINES SINGULARITY_USE_HDF)
endif() # SINGULARITY_USE_HDF5

#NB: An edge-case (but an easy one), that doesn't slot
#   easily into the `import_dependency()` setup.
#   TODO: just do it correctly
#TODO: Check: Does EOSPAC only pass into side-targets?
if (SINGULARITY_USE_EOSPAC)
  #TODO: should be strictly module mode
  find_package(EOSPAC REQUIRED)
  list(APPEND SINGULARITY_PUBLIC_DEFINES SINGULARITY_USE_EOSPAC)
  list(APPEND SINGULARITY_PUBLIC_LIBS EOSPAC::eospac)
  list(APPEND SINGULARITY_PUBLIC_LIBS eospac-wrapper)
  list(APPEND SINGULARITY_DEPS_PKGS EOSPAC)

  add_subdirectory(eospac-wrapper)
endif ()


#TODO: this feels out-of-place here, though i haven't crawled into the source yet to see
# what the point is.
if (SINGULARITY_TEST_SESAME)
  list(APPEND SINGULARITY_PRIVATE_DEFINES SINGULARITY_TEST_SESAME)
endif()

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
    $<INSTALL_INTERFACE:include>
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

foreach(_sdir ${subdirs2add})
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


