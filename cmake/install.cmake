# ------------------------------------------------------------------------------#
# © 2021-2023. Triad National Security, LLC. All rights reserved.  This program
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

include(GNUInstallDirs)

# ----------------------------------------------------------------------------#
# Generate config file
# ----------------------------------------------------------------------------#

# get all available `SINGULARITY_` cmake variables set during configuration
get_cmake_property(_variableNames VARIABLES)
string(REGEX MATCHALL "(^|;)SINGULARITY_[A-Za-z0-9_]*" _matchedVars
             "${_variableNames}")

# use config template to generate the configuration of the build not sure why
# CMake doesn't do this automatically, but w/e
foreach(_variableName ${_matchedVars})
  set(SINGULARITY_EOS_CONFIG_CODE
      "${SINGULARITY_EOS_CONFIG_CODE}\nset(${_variableName} \"${${_variableName}}\")"
  )
endforeach()

# generate the configuration file NOTE: for some lunatic reason,
# `CMAKE_FILES_DIRECTORY` is set with a leading `/`
configure_file(
  ${PROJECT_SOURCE_DIR}/config/singularity-eosConfig.cmake.in
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/singularity-eosConfig.cmake @ONLY)

install(
  FILES ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/singularity-eosConfig.cmake
  DESTINATION
    ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/cmake/singularity-eos)

# ----------------------------------------------------------------------------#
# install cmake files
# ----------------------------------------------------------------------------#

install(DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/singularity-eos
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/cmake)

install(FILES ${PROJECT_SOURCE_DIR}/cmake/FindEOSPAC.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/cmake)

# ----------------------------------------------------------------------------#
# install export target
# ----------------------------------------------------------------------------#
install(
  TARGETS singularity-eos_Common
  EXPORT singularity-eos_Common
  DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(
  EXPORT singularity-eos_Common
  FILE singularity-eos_Common.cmake
  NAMESPACE "singularity-eos::"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/singularity-eos)

install(
  TARGETS singularity-eos_Interface
  EXPORT singularity-eos_Interface
  DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(
  EXPORT singularity-eos_Interface
  FILE singularity-eos_Interface.cmake
  NAMESPACE "singularity-eos::"
  COMPONENT singularity-eos_Interface
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/singularity-eos)

if(TARGET singularity-eos_Library)
  install(
    TARGETS singularity-eos_Library
    EXPORT singularity-eos_Library
    DESTINATION ${CMAKE_INSTALL_LIBDIR})

  install(
    EXPORT singularity-eos_Library
    FILE singularity-eos_Library.cmake
    NAMESPACE "singularity-eos::"
    COMPONENT singularity-eos_Library
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/singularity-eos)
endif()

# ----------------------------------------------------------------------------#
# Install headers
# ----------------------------------------------------------------------------#

# install singularity-eos headers
get_property(install_headers GLOBAL PROPERTY EOS_INSTALL_HEADERS)
foreach(src dst IN ZIP_LISTS eos_headers install_headers)
  get_filename_component(DIR ${dst} DIRECTORY)
  install(FILES ${src} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${DIR})
endforeach() # file

# Special sauce so generated file has proper include path
install(FILES ${CMAKE_BINARY_DIR}/generated/singularity-eos/eos/eos.hpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/singularity-eos/eos)

# install the fortran modules NB: cmake doesn't provide a clean way to handle
if(SINGULARITY_USE_FORTRAN)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran/
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/singularity-eos/eos)
endif()

# ----------------------------------------------------------------------------#
# local export
# ----------------------------------------------------------------------------#

# same as install step, but just places the file in the build tree. useful for
# downstream projects that use the source directly
export(
  EXPORT singularity-eos_Common
  FILE ${CMAKE_CURRENT_BINARY_DIR}/singularity-eos_Common.cmake
  NAMESPACE singularity-eos::)

export(
  EXPORT singularity-eos_Interface
  FILE ${CMAKE_CURRENT_BINARY_DIR}/singularity-eos_Interface.cmake
  NAMESPACE singularity-eos::)

if(TARGET singularity-eos_Library)
  export(
    EXPORT singularity-eos_Library
    FILE ${CMAKE_CURRENT_BINARY_DIR}/singularity-eos_Library.cmake
    NAMESPACE singularity-eos::)
endif()

# ----------------------------------------------------------------------------#
# Data files
# ----------------------------------------------------------------------------#

# install data files needed for various eos models
if(SINGULARITY_USE_HELMHOLTZ)
  # install data files
  install(DIRECTORY ${PROJECT_SOURCE_DIR}/data/helmholtz
          DESTINATION ${CMAKE_INSTALL_DATADIR}/singularity-eos/data)
endif()
