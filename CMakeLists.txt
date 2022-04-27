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

project(singularity-eos VERSION 1.5.0)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

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
option (SINGULARITY_FMATH_USE_ORDER_4 "4th order interpolant for fast logs. Default is 7th order." OFF)
option (SINGULARITY_FMATH_USE_ORDER_5 "5th order interpolant for fast logs. Default is 7th order." OFF)

if (SINGULARITY_SUBMODULE_MODE)
  set(SINGULARITY_BETTER_DEBUG_FLAGS OFF CACHE BOOL "" FORCE)
  set(SINGULARITY_HIDE_MORE_WARNINGS ON CACHE BOOL "" FORCE)
endif()

# Patches variant to be compatible with cuda
# Assumes "patch" is present on system
message(STATUS "Patching mpark::variant to support GPUs")
execute_process(COMMAND patch -N -s -V never
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/variant/include/mpark/variant.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/cuda_compatibility.patch
)

#------------------------------------------------------------------------------#
# Compiler & language setup
#------------------------------------------------------------------------------#

# checks if this is our build, or we've 
# been imported via `add_subdirectory`
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

enable_language(CXX)
include(CMakeDetermineCXXCompiler)
if(SINGULARITY_USE_FORTRAN)
  enable_language(Fortran)
  include(CMakeDetermineFortranCompiler)
endif()

# NOTE: this may not work with kokkos, we'll see.
#if(SINGULARITY_USE_CUDA)
#  enable_language(CUDA)
#  include(CMakeDetermineCUDACompiler)
#endif()

#------------------------------------------------------------------------------#
# De-thaw some options
#------------------------------------------------------------------------------#

# using this space for "immediate" logic that can be avoid messages getting
# lost inside the brrrrrrrr of cmake.
if(SINGULARITY_USE_KOKKOSKERNELS AND NOT SINGULARITY_USE_KOKKOS)
  message(FATAL_ERROR "Kokkos is required for KokkosKernels.")
endif()

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

#------------------------------------------------------------------------------#
# Process some options
#------------------------------------------------------------------------------#

#TODO: process. i want to try to 'reduce' some of the boolean logic of before
# to single T/F vars, e.g.
#   before:
#     if(SINGULARITY_BUILD_SESAME2SPINER)
#       if(SINGULARITY_USE_EOSPAC)
#         if(SINGULARITY_USE_HDF5)
#           add_subdirectory(${PROJECT_SOURCE_DIR}/sesame2spiner)
#   after:
#   # in this section
#     if(...) ...
#       list(APPEND _subdirs2add sesame2spiner)
#     ...
#   # later...
#     foreach(_sdir ${_subdirs2add})
#       add_subdirectory(${_sdir})
#     endforeach()
# this way, we can decouple the input processing from the configuration 
# execution. 

#------------------------------------------------------------------------------#
# Handle dependencies
#------------------------------------------------------------------------------#

include(cmake/dependency.cmake)

# keep config data stored in easily accessable vars,
# only attaching them to targets as late as possible.
# this allows for easier handling and manipulating, as
# well as later info printing and debugging
set(TAR_PUBLIC_DEFINES "")  # list to `target_compile_definitions()`
set(TAR_PRIVATE_DEFINES "") # """
set(TAR_PUBLIC_LIB "")      # list to `target_link_libraries()`
set(TAR_PRIVATE_LIB "")     # """
set(TAR_COMPILEFLAGS "")    # TODO: actually use these
set(TAR_LIB_DEPS "")        # a list of packages, used to generate import rules for installed targets. (may want to rename)

if(SINGULARITY_USE_KOKKOS)
  list(APPEND TAR_PUBLIC_DEFINES PORTABILITY_STRATEGY_KOKKOS)
  append_target_dependency(TAR_PUBLIC_LIB TAR_LIB_DEPS
    PKG Kokkos
    TARGETS Kokkos::kokkos
    )
endif() #SINGULARITY_USE_KOKKOS

if(SINGULARITY_USE_KOKKOSKERNELS)
  list(APPEND TAR_PUBLIC_DEFINES SINGULARITY_USE_KOKKOSKERNELS)
  append_target_dependency(TAR_PUBLIC_LIB TAR_LIB_DEPS
    PKG KokkosKernels
    TARGETS Kokkos::kokkoskernels
  )
else()
  append_target_dependency(TAR_PUBLIC_LIB TAR_LIB_DEPS
    PKG Eigen3
    TARGETS Eigen3::Eigen
    SUBDIR ${PROJECT_SOURCE_DIR}/utils/eigen
  )
endif() #SINGULARITY_USE_KOKKOSKERNELS

append_target_dependency(TAR_PUBLIC_LIB TAR_LIB_DEPS
  PKG mpark_variant
  TARGETS mpark_variant
  SUBDIR ${PROJECT_SOURCE_DIR}/utils/variant
)

append_target_dependency(TAR_PUBLIC_LIB TAR_LIB_DEPS
  PKG spiner
  TARGETS spiner::spiner
)

append_target_dependency(TAR_PUBLIC_LIB TAR_LIB_DEPS
  PKG ports-of-call
  TARGETS ports-of-call::ports-of-call
)

# only use hdf5 if spiner
if (SINGULARITY_USE_HDF5)
  list(APPEND TAR_PUBLIC_DEFINES SPINER_USE_HDF)
  list(APPEND TAR_PUBLIC_DEFINES SINGULARITY_USE_HDF)
  append_target_dependency(TAR_PUBLIC_LIB TAR_LIB_DEPS
    PKG HDF5
    TARGETS HDF5::HDF5 hdf5::hdf5_hl
    COMPONENTS C HL
  )
  if(HDF5_IS_PARALLEL)
    append_target_dependency(TAR_PUBLIC_LIB TAR_LIB_DEPS
      PKG MPI
      TARGETS MPI::MPI_CXX
      COMPONENTS CXX
    )
  endif()
endif()

#TODO: this feels out-of-place here, though i haven't crawled into the source yet to see
# what the point is.
if (SINGULARITY_TEST_SESAME)
  list(APPEND TAR_PUBLIC_DEFINES SINGULARITY_TEST_SESAME)
endif()

message(STATUS "[== INFO: status so far ===============================")
message(STATUS "[== dependency list: ${TAR_LIB_DEPS}")
message(STATUS "[== library targets: ${TAR_PUBLIC_LIB} ${TAR_PRIVATE_LIB}")
message(STATUS "[== compiler definitions: ${TAR_PUBLIC_DEFINES} ${TAR_PRIVATE_DEFINES}")
message(STATUS "[======================================================")
#==============================================================================#
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

target_link_libraries(singularity-eos
  PUBLIC
    ${TAR_PUBLIC_LIB}
  PRIVATE
    ${TAR_PRIVATE_LIB}
)

target_compile_definitions(singularity-eos
  PUBLIC
    ${TAR_PUBLIC_DEFINES}
  PRIVATE
    ${TAR_PRIVATE_DEFINES}
)

#------------------------------------------------------------------------------#
# "Unincorperated real-estate"
# this section includes cmake that hasn't been reworked, and is largely 
# unmodified at present.
# TODO: for the sake of harmony, these bits could be given a polish. however,
# not much would functionally change.
#------------------------------------------------------------------------------#

# xl fix
#TODO: should these be private? I think so, but not sure
target_compile_options(singularity-eos
  PRIVATE
    $<$<COMPILE_LANG_AND_ID:CXX,XL>:-std=c++1y;-qxflag=disable__cplusplusOverride>)
target_link_options(singularity-eos
  PRIVATE
  $<$<COMPILE_LANG_AND_ID:CXX,XL>:-std=c++1y;-qxflag=disable__cplusplusOverride>)

# cuda specific flags
if(SINGULARITY_USE_CUDA)
  if(SINGULARITY_HIDE_MORE_WARNINGS)
    target_compile_options(singularity-eos
       PUBLIC # Generator expression shamelessly copied from EAP
         "$<$<COMPILE_LANGUAGE:CXX>:--expt-relaxed-constexpr;>"
    )
  else()
    target_compile_options(singularity-eos
       PUBLIC # Generator expression shamelessly copied from EAP
         "$<$<COMPILE_LANGUAGE:CXX>:--expt-relaxed-constexpr;-Xcudafe;--diag_suppress=esa_on_defaulted_function_ignored;>"
    )
  endif()
  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(singularity-eos
      PUBLIC
        -use_fast_math
    )
  elseif(CMAKE_BUILD_TYPE STREQUAL "Debug" AND SINGULARITY_BETTER_DEBUG_FLAGS)
    target_compile_options(singularity-eos
      PUBLIC
        "$<$<COMPILE_LANGUAGE:CXX>:-G;-lineinfo;>"
    )
  endif()
endif()

#TODO: Check: Does EOSPAC only pass into side-targets?
if (SINGULARITY_USE_EOSPAC)
  #TODO: should be strictly module mode
  find_package(EOSPAC REQUIRED)
  list(APPEND TAR_LIB_DEPS EOSPAC)

  add_subdirectory(eospac-wrapper)
  target_link_libraries(singularity-eos
    PUBLIC
    eospac-wrapper)
  target_compile_definitions(singularity-eos
    PUBLIC
    SINGULARITY_USE_EOSPAC)
endif ()

# fast math
if (SINGULARITY_USE_SINGLE_LOGS)
  target_compile_definitions(singularity-eos PRIVATE
      SINGULARITY_USE_SINGLE_LOGS)
endif()
if (SINGULARITY_FMATH_USE_ORDER_4 AND SINGULARITY_FMATH_USE_ORDER_5)
  message(FATAL_ERROR "Order 4 and order 5 interpolation both specified. Please specify only one or zero. If no specification is made, order 7 will be used.")
endif()
if (SINGULARITY_FMATH_USE_ORDER_4)
  target_compile_definitions(singularity-eos PRIVATE
    SINGULARITY_FMATH_USE_ORDER_4)
endif()
if (SINGULARITY_FMATH_USE_ORDER_5)
  target_compile_definitions(singularity-eos PRIVATE
    SINGULAIRTY_FMATH_USE_ORDER_5)
endif()

if (SINGULARITY_USE_HDF5)
  if(SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER)
    add_subdirectory(stellarcollapse2spiner)
  endif()
  if(SINGULARITY_USE_EOSPAC AND SINGULARITY_BUILD_SESAME2SPINER)
    add_subdirectory(sesame2spiner)
  endif()
endif()

if (SINGULARITY_BUILD_EXAMPLES)
  add_subdirectory(example)
endif()


#----------------------------------------------------------------------------#
# Generate config file
#----------------------------------------------------------------------------#

# get all available `SINGULARITY_` cmake variables set during configuration
get_cmake_property(_variableNames VARIABLES)
string (REGEX MATCHALL "(^|;)SINGULARITY_[A-Za-z0-9_]*"
  _matchedVars "${_variableNames}")

# use config template to generate the configuration of the build
# not sure why CMake doesn't do this automatically, but w/e
foreach(_variableName ${_matchedVars})
  set(SINGULARITY_CONFIG_CODE
    "${SINGULARITY_CONFIG_CODE}\nset(${_variableName} \"${${_variableName}}\")")
endforeach()

# for all the upstream packages collected in the configuration, replicate
# the `find_package()` in the install configuration. note that
# if a package was consumed with `add_subdirectory` in configuration,
# the downstream `find_package()` should find it's way to the consumed
# package install path
foreach(_depName ${TAR_LIB_DEPS})
  set(_components "")
  if(TAR_LIB_DEPS_${_depName})
    set(_components "COMPONENTS ${TAR_LIB_DEPS_${_depName}}")
  endif()
  set(SINGULARITY_CONFIG_DEPENDENCIES
    "${SINGULARITY_CONFIG_DEPENDENCIES}\nfind_package(${_depName} ${_components} REQUIRED)")
endforeach()

# generate the configuration file
# NOTE: for some lunatic reason, `CMAKE_FILES_DIRECTORY` is set with a leading `/`
configure_file(${PROJECT_SOURCE_DIR}/config/singularity-eosConfig.cmake.in
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/singularity-eosConfig.cmake @ONLY)

#----------------------------------------------------------------------------#
# Install library
#----------------------------------------------------------------------------#

include(GNUInstallDirs)

# declare library install, and associate it with the export targets.
# NOTE: the `DESTINATION` here is the implicit default, tho I think we 
# *have* to make it explicit since we (may) also install a fortran module.
install(
  TARGETS singularity-eos 
  EXPORT singularity-eosTargets
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# export install; generates CMake code that instructs other projects how to 
# import targets from this source tree.
install(
  EXPORT singularity-eosTargets
  FILE  singularity-eosTargets.cmake
  NAMESPACE "singularity-eos::"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/singularity-eos
)

#----------------------------------------------------------------------------#
# Install headers/cmake modules + Packaging
#----------------------------------------------------------------------------#

# install singularity-eos headers
foreach(file ${_install_headers})
  get_filename_component(DIR ${file} DIRECTORY)
  install(
    FILES singularity-eos/${file}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/singularity-eos/${DIR}
  )
endforeach()

# install the generated config file
install(
  FILES ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/singularity-eosConfig.cmake
  DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}/cmake/singularity-eos
)

# install cmake modules
# NOTE: I don't appear to need `FindPortsofCall.cmake`, probably due to recent updates to that code
install(
  FILES 
    ${PROJECT_SOURCE_DIR}/cmake/FindEOSPAC.cmake
    ${PROJECT_SOURCE_DIR}/cmake/FindHDF5.cmake
  DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}/cmake/singularity-eos
)

#----------------------------------------------------------------------------#
# Local export
#----------------------------------------------------------------------------#

# same as install step, but just places the file in the build tree.
# useful for downstream projects that use the source directly
export(
  EXPORT    singularity-eosTargets
  FILE      ${CMAKE_CURRENT_BINARY_DIR}/cmake/singularity-eosTargets.cmake
  NAMESPACE singularity-eos::
)

# apparently, this doesn't do anything anymore,
# see https://cmake.org/cmake/help/latest/policy/CMP0090.html
# to re-enable, uncomment 
## set(CMAKE_EXPORT_PACKAGE_REGISTRY "ON")
export(PACKAGE singularity-eos)

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
