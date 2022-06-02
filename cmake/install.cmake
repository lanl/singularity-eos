#------------------------------------------------------------------------------
# placeholder
#------------------------------------------------------------------------------

include(GNUInstallDirs)

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
foreach(_depName ${SINGULARITY_DEP_PKGS})
  set(_components "")
  if(SINGULARITY_DEP_PKGS_${_depName})
    set(_components "COMPONENTS ${SINGULARITY_DEP_PKGS_${_depName}}")
  endif()
  set(SINGULARITY_CONFIG_DEPENDENCIES
    "${SINGULARITY_CONFIG_DEPENDENCIES}\nfind_package(${_depName} ${_components} REQUIRED)")
endforeach()

# generate the configuration file
# NOTE: for some lunatic reason, `CMAKE_FILES_DIRECTORY` is set with a leading `/`
configure_file(${PROJECT_SOURCE_DIR}/config/singularity-eosConfig.cmake.in
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/singularity-eosConfig.cmake @ONLY)

#----------------------------------------------------------------------------#
# Install library + attach to export
#----------------------------------------------------------------------------#

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
endforeach() # file

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

# MAUNEYC
# NB: extra step b/c kokkos doesn't do this? seems weird
# to the best of my reading, this is just bookkeeping after
# configuration step, but without it, there is a conflict between
# export sets. I'm likely doing something incorrect.
# if we are using Kokkos as a submodule
if(SINGULARITY_USE_KOKKOS AND NOT Kokkos_FOUND)
  export(
    EXPORT KokkosTargets 
    FILE ${CMAKE_CURRENT_BINARY_DIR}/cmake/KokkosTargets.cmake 
    NAMESPACE Kokkos::
  )
endif() # USE_KOKKOS AND NOT Kokkos_FOUND

if(SINGULARITY_USE_KOKKOSKERNELS AND NOT KokkosKernels_FOUND)
  export(
    EXPORT KokkosKernelsTargets 
    FILE ${CMAKE_CURRENT_BINARY_DIR}/cmake/KokkosKernelsTargets.cmake 
    NAMESPACE Kokkos::
  )
endif() # USE_KOKKOS AND NOT Kokkos_FOUND

# apparently, this doesn't do anything anymore,
# see https://cmake.org/cmake/help/latest/policy/CMP0090.html
# to re-enable, uncomment 
## set(CMAKE_EXPORT_PACKAGE_REGISTRY "ON")
export(PACKAGE singularity-eos)

