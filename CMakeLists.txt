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
cmake_minimum_required(VERSION 3.19)

project(singularity-eos VERSION 1.5.0)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Patches variant to be compatible with cuda
# Assumes "patch" is present on system
message(STATUS "Patching mpark::variant to support GPUs")
execute_process(COMMAND patch -N -s -V never
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/variant/include/mpark/variant.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/cuda_compatibility.patch
)

include(compilers)
include(cuda)
include(hdf5)
include(kokkos)
include(kokkoskernels)
include(eigen3)
include(eospac)
include(portsofcall)
include(variant)
include(mathselect)

# 20210930.mauneyc
# This allows the exported "package" (singularity-eos)
# to be automatically found by other projects.
# NOTE: This is disabled by default, may want to
#       control at command line
set(CMAKE_EXPORT_PACKAGE_REGISTRY "ON")
include(library)

# all-in-on library configuration
add_library_target(${PROJECT_NAME} singularity-eos
  EXPORT_TARGET
    ${PROJECT_NAME}Targets
  VERSION
    ${PROJECT_VERSION}
  NAMESPACE
    ${PROJECT_NAME}
  DEFINE_PUBLIC
    ${TPL_DEFINES}
  INCLUDE_PUBLIC
    ${TPL_INCLUDES}
  INCLUDE_PRIVATE
    ${CMAKE_BINARY_DIR}
  LINK_PUBLIC
    ${TPL_LIBRARIES}
  OPTION_PRIVATE
    ${XTRA_BUILD_FLAGS}
)

install(
  FILES
    ${PROJECT_SOURCE_DIR}/cmake/compilers.cmake
    ${PROJECT_SOURCE_DIR}/cmake/cuda.cmake
    ${PROJECT_SOURCE_DIR}/cmake/eigen3.cmake
    ${PROJECT_SOURCE_DIR}/cmake/eospac.cmake
    ${PROJECT_SOURCE_DIR}/cmake/FindEigen3.cmake
    ${PROJECT_SOURCE_DIR}/cmake/FindEOSPAC.cmake
    ${PROJECT_SOURCE_DIR}/cmake/FindHwloc.cmake
    ${PROJECT_SOURCE_DIR}/cmake/hdf5.cmake
    ${PROJECT_SOURCE_DIR}/cmake/kokkos.cmake
    ${PROJECT_SOURCE_DIR}/cmake/kokkoskernels.cmake
    ${PROJECT_SOURCE_DIR}/cmake/library.cmake
    ${PROJECT_SOURCE_DIR}/cmake/mathselect.cmake
    ${PROJECT_SOURCE_DIR}/cmake/portsofcall.cmake
    ${PROJECT_SOURCE_DIR}/cmake/subdirlist.cmake
    ${PROJECT_SOURCE_DIR}/cmake/variant.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

if(SINGULARITY_USE_HDF5)
  if(SINGULARITY_USE_EOSPAC AND SINGULARITY_BUILD_SESAME2SPINER)
    add_subdirectory(sesame2spiner)
  endif()
  if(SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER)
    add_subdirectory(stellarcollapse2spiner)
  endif()
endif()



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
