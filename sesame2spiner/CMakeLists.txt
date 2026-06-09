#------------------------------------------------------------------------------
# © 2021-2025. Triad National Security, LLC. All rights reserved.  This
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


# ------------------------------------------------------------------------------
# sesame2spiner Library
# ------------------------------------------------------------------------------

add_library(sesame2spiner-lib
  sesame2spiner/io_eospac.cpp
  sesame2spiner/io_eospac.hpp
  sesame2spiner/generate_files.cpp
  sesame2spiner/generate_files.hpp
  sesame2spiner/parser.cpp
  sesame2spiner/parser.hpp
)

# Set library output name (without -lib suffix)
set_target_properties(sesame2spiner-lib PROPERTIES OUTPUT_NAME sesame2spiner)

target_include_directories(sesame2spiner-lib
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>
)

# sesame2spiner library depends on singularity-utils (not full singularity-eos)
target_link_libraries(sesame2spiner-lib
  PUBLIC
    singularity-utils::singularity-utils
    eospac-wrapper
    singularity-eos::singularity-eos_Common  # For HDF5
)

# Version information for library
target_compile_definitions(sesame2spiner-lib
  PUBLIC
  SESAME2SPINER_VERSION=\"${PROJECT_VERSION}\"
)

# Create alias for consistency
add_library(sesame2spiner::sesame2spiner ALIAS sesame2spiner-lib)

# ------------------------------------------------------------------------------
# sesame2spiner Executable
# ------------------------------------------------------------------------------

add_executable(sesame2spiner
  src/parse_cli.cpp
  src/parse_cli.hpp
  src/main.cpp
)

target_include_directories(sesame2spiner
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(sesame2spiner
  PRIVATE
    sesame2spiner-lib
)

# ------------------------------------------------------------------------------
# Installation
# ------------------------------------------------------------------------------

install(TARGETS sesame2spiner DESTINATION ${CMAKE_INSTALL_BINDIR})


install(DIRECTORY sesame2spiner/
  DESTINATION include/sesame2spiner
)

