#------------------------------------------------------------------------------
# © 2021-2026. Triad National Security, LLC. All rights reserved.  This
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


# sesame2spiner library
add_library(sesame2spiner-lib
  sesame2spiner/io_eospac.cpp
  sesame2spiner/io_eospac.hpp
  sesame2spiner/generate_files.cpp
  sesame2spiner/generate_files.hpp
  sesame2spiner/parser.cpp
  sesame2spiner/parser.hpp

)

set_target_properties(sesame2spiner-lib PROPERTIES OUTPUT_NAME sesame2spiner)

target_include_directories(sesame2spiner-lib
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(sesame2spiner-lib
  PUBLIC
    singularity-utils::singularity-utils
    eospac-wrapper
    singularity-eos::singularity-eos_Common  # For HDF5
)

target_compile_definitions(sesame2spiner-lib
  PUBLIC
  SESAME2SPINER_VERSION=\"${PROJECT_VERSION}\"
)

# Propate SPINER defines
if (DEFINED SINGULARITY_USE_SPINER)
  target_compile_definitions(sesame2spiner-lib
  PUBLIC
      $<$<BOOL:${SINGULARITY_USE_SPINER}>:SINGULARITY_USE_SPINER>
  )
endif()

# Propagate HDF5 usage definitions
# These must match the settings used by singularity-eos
if (DEFINED SINGULARITY_USE_SPINER_WITH_HDF5)
  target_compile_definitions(sesame2spiner-lib
  PUBLIC
      $<$<BOOL:${SINGULARITY_USE_SPINER_WITH_HDF5}>:SINGULARITY_USE_SPINER_WITH_HDF5>
  )
endif()


add_library(sesame2spiner::sesame2spiner ALIAS sesame2spiner-lib)

# sesame2spiner executable
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
  PUBLIC
    sesame2spiner-lib
)

# install

install(TARGETS sesame2spiner DESTINATION ${CMAKE_INSTALL_BINDIR})

install(TARGETS sesame2spiner-lib
  EXPORT sesame2spiner-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION include
)

install(DIRECTORY sesame2spiner/
  DESTINATION include/sesame2spiner
)

if(NOT SINGULARITY_SUBMODULE_MODE)
  install(EXPORT sesame2spiner-targets
    FILE sesame2spiner-targets.cmake
    NAMESPACE sesame2spiner::
    DESTINATION lib/cmake/sesame2spiner
  )
endif()


# TODO: Add tests for sesame2spiner here.
