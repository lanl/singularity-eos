#------------------------------------------------------------------------------
# placeholder
#------------------------------------------------------------------------------

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
  singularity_msg(STATUS "NOTE: we export `KokkosTargets` here, because the `Kokkos` package was imported through `add_subdirectory()`, and that does not export these on it's own")
  export(
    EXPORT KokkosTargets 
    FILE ${CMAKE_CURRENT_BINARY_DIR}/cmake/KokkosTargets.cmake 
    NAMESPACE Kokkos::
  )
endif() # USE_KOKKOS AND NOT Kokkos_FOUND

# apparently, this doesn't do anything anymore,
# see https://cmake.org/cmake/help/latest/policy/CMP0090.html
# to re-enable, uncomment 
## set(CMAKE_EXPORT_PACKAGE_REGISTRY "ON")
export(PACKAGE singularity-eos)

