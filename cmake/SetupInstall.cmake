include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(INSTALL_CONFIG_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/singularity-eos)

install(TARGETS eos
  EXPORT singularity-eos_targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/singularity-eos/include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT singularity-eos_targets
  FILE
    singularity-eosTargets.cmake
  NAMESPACE
    singularity-eos::
  DESTINATION
    ${INSTALL_CONFIG_DIR}
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/singularity-eosConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
  ${CMAKE_CURRENT_LIST_DIR}/singularity-eosConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/singularity-eosConfig.cmake
  INSTALL_DESTINATION ${INSTALL_CONFIG_DIR}
)

#Install the config, configversion and custom find modules
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/singularity-eosConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/singularity-eosConfigVersion.cmake
  DESTINATION ${INSTALL_CONFIG_DIR}
)


export(
  EXPORT singularity-eos_targets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/singularity-eosTargets.cmake
  NAMESPACE singularity-eos::
)

#Register package in the User Package Registry
export(PACKAGE singularity-eos)
