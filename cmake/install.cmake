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

