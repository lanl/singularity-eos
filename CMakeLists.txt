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
                ${CMAKE_CURRENT_SOURCE_DIR}/utils/cuda_compatibility.patch)


include(cuda)
include(hdf5)
include(kokkos)
include(kokkoskernels)
include(eigen3)
include(eospac)
include(portsofcall)

#----------------------------------------------------------------------------#
# Capture use settings.
#----------------------------------------------------------------------------#

# get_cmake_property(_variableNames VARIABLES)
# string (REGEX MATCHALL "(^|;)USE_[A-Za-z0-9_]*"
#   _matchedVars "${_variableNames}")

# foreach(_variableName ${_matchedVars})
#   set(SINGULARITY_${_variableName} ${${_variableName}})
# endforeach()

#----------------------------------------------------------------------------#
# Configure header.
#----------------------------------------------------------------------------#

configure_file(${PROJECT_SOURCE_DIR}/config/${PROJECT_NAME}-config.h.in
  ${CMAKE_BINARY_DIR}/flecsi-config.h @ONLY)

install(
  FILES ${CMAKE_BINARY_DIR}/${PROJECT_NAME}-config.h
  DESTINATION include
)

include(library)

add_library_target(${PROJECT_NAME} singularity-eos
  EXPORT_TARGET
    ${PROJECT_NAME}Targets
  VERSION
    ${PROJECT_VERSION}
  NAMESPACE
    ${PROJECT_NAME}
  DEFINE_PRIVATE
    ${TPL_DEFINES}
  INCLUDE_PUBLIC
    ${TPL_INCLUDES}
  INCLUDE_PRIVATE
    ${CMAKE_BINARY_DIR}
  LINK_PUBLIC
    ${TPL_LIBRARIES}
)

#----------------------------------------------------------------------------#
# Extract all project options so they can be exported to the
# FleCSIConfig.cmake file.
#----------------------------------------------------------------------------#

get_cmake_property(_variableNames VARIABLES)
string (REGEX MATCHALL "(^|;)SINGULARITY_[A-Za-z0-9_]*"
  _matchedVars "${_variableNames}")

foreach(_variableName ${_matchedVars})
  string(REGEX REPLACE "SINGULARITY" "singularity-eos" _tmpvariableName ${_variableName})
  set(SINGULARITY_CONFIG_CODE
    "${SINGULARITY_CONFIG_CODE}\nset(${_tmpvariableName} \"${${_variableName}}\")")
endforeach()

configure_file(${PROJECT_SOURCE_DIR}/config/${PROJECT_NAME}Config.cmake.in
${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}.cmake @ONLY)

install(
FILES ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}.cmake
DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIBDIR}/cmake/${PROJECT_NAME}
)

# cmakes
configure_file(${PROJECT_SOURCE_DIR}/config/${PROJECT_NAME}CMakeConfig.cmake.in
${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}CMakeConfig.cmake @ONLY)

install(
FILES ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}CMakeConfig.cmake
DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIBDIR}/cmake/${PROJECT_NAME}CMake
)

install(
FILES ${PROJECT_SOURCE_DIR}/cmake/cuda.cmake
FILES ${PROJECT_SOURCE_DIR}/cmake/eigen3.cmake
FILES ${PROJECT_SOURCE_DIR}/cmake/eospac.cmake
FILES ${PROJECT_SOURCE_DIR}/cmake/FindEigen3.cmake
FILES ${PROJECT_SOURCE_DIR}/cmake/FindEOSPAC.cmake
FILES ${PROJECT_SOURCE_DIR}/cmake/FindHwloc.cmake
FILES ${PROJECT_SOURCE_DIR}/cmake/FindPortsofCall.cmake
FILES ${PROJECT_SOURCE_DIR}/cmake/hdf5.cmake
FILES ${PROJECT_SOURCE_DIR}/cmake/kokkos.cmake
FILES ${PROJECT_SOURCE_DIR}/cmake/kokkoskernels.cmake
FILES ${PROJECT_SOURCE_DIR}/cmake/library.cmake
FILES ${PROJECT_SOURCE_DIR}/cmake/portsofcall.cmake
FILES ${PROJECT_SOURCE_DIR}/cmake/subdirlist.cmake
FILES ${PROJECT_SOURCE_DIR}/cmake/variant.cmake
DESTINATION ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/cmake
)


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
