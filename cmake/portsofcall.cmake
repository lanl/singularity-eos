# no options
include(GNUInstallDirs)
add_library(PortsofCall INTERFACE IMPORTED)

target_include_directories(PortsofCall INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utils>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/ports-of-call>
)

list(APPEND TPL_LIBRARIES PortsofCall)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/utils/ports-of-call
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ports-of-call)


