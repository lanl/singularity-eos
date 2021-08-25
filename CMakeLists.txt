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


add_library(${PROJECT_NAME}::flags INTERFACE IMPORTED GLOBAL)

include(ProcessConfiguration)
include(SetupCompilers)
include(SetupFlags)

add_subdirectory(${PROJECT_SOURCE_DIR}/singularity-eos)

include(SetupInstall)

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
