#------------------------------------------------------------------------------
# © 2026. Triad National Security, LLC. All rights reserved.  This
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
#------------------------------------------------------------------------------

# This file was made in part with generative AI


# Create interface library (header-only)
add_library(singularity-utils INTERFACE)
add_library(singularity-utils::singularity-utils ALIAS singularity-utils)

# Set include directories
target_include_directories(singularity-utils
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

# singularity-utils depends on ports-of-call and spiner
# These should already be found by the parent singularity-eos project
target_link_libraries(singularity-utils
  INTERFACE
    ports-of-call::ports-of-call
)

# Propate SPINER defines
if (SINGULARITY_USE_SPINER)
  target_link_libraries(singularity-utils INTERFACE
      spiner::spiner
  )
  target_compile_definitions(singularity-utils
    INTERFACE
      $<$<BOOL:${SINGULARITY_USE_SPINER}>:SINGULARITY_USE_SPINER>
  )
endif()

# Propagate HDF5 usage definitions
# These must match the settings used by singularity-eos
if (SINGULARITY_USE_SPINER_WITH_HDF5)
  target_compile_definitions(singularity-utils
    INTERFACE
      $<$<BOOL:${SINGULARITY_USE_SPINER_WITH_HDF5}>:SINGULARITY_USE_SPINER_WITH_HDF5>
  )
endif()


# Propagate compile definitions for FastMath configuration
# These must match the settings used by singularity-eos
if (SINGULARITY_USE_TRUE_LOG_GRIDDING)
  target_compile_definitions(singularity-utils
    INTERFACE
      $<$<BOOL:${SINGULARITY_USE_TRUE_LOG_GRIDDING}>:SINGULARITY_USE_TRUE_LOG_GRIDDING>
  )
endif()

if (SINGULARITY_USE_SINGLE_LOGS)
  target_compile_definitions(singularity-utils
    INTERFACE
      $<$<BOOL:${SINGULARITY_USE_SINGLE_LOGS}>:SINGULARITY_USE_SINGLE_LOGS>
  )
endif()

if (SINGULARITY_NQT_PORTABLE)
  target_compile_definitions(singularity-utils
    INTERFACE
      $<$<BOOL:${SINGULARITY_NQT_PORTABLE}>:SINGULARITY_NQT_PORTABLE>
  )
endif()

if (SINGULARITY_NQT_ORDER_1)
  target_compile_definitions(singularity-utils
    INTERFACE
      $<$<BOOL:${SINGULARITY_NQT_ORDER_1}>:SINGULARITY_NQT_ORDER_1>
  )
endif()

# Version information
target_compile_definitions(singularity-utils
  INTERFACE
    SPINER_FORMAT_VERSION="${PROJECT_VERSION}"
    SPINER_FORMAT_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
    SPINER_FORMAT_VERSION_MINOR=${PROJECT_VERSION_MINOR}
    SPINER_FORMAT_VERSION_PATCH=${PROJECT_VERSION_PATCH}
)

# Installation rules
# singularity-utils is always built as part of singularity-eos via add_subdirectory,
# so we don't need the full package config machinery. Just install the headers.
if(NOT SINGULARITY_SUBMODULE_MODE)
  install(DIRECTORY singularity-utils/
    DESTINATION include/singularity-utils
  )
endif()
