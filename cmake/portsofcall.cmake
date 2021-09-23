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

# when portsofcall gets a proper CMake, we can forego this
# find_path(PortsofCall_INCLUDE_DIR portability.hpp
#     HINTS ${PROJECT_SOURCE_DIR}/utils/ports-of-call
#     DOC "Directory where Ports-of-Call header files are located"
# )


# include(FindPackageHandleStandardArgs)
# find_package_handle_standard_args(PortsofCall
#     REQUIRED_VARS PortsofCall_INCLUDE_DIR
# )

# if(PortsofCall_FOUND)
#     set(PortsofCall_INCLUDE_DIRS ${PortsofCall_INCLUDE_DIR})
#     add_library(PortsofCall::PortsofCall INTERFACE IMPORTED)
#     set_target_properties(PortsofCall::PortsofCall PROPERTIES
# 		      INTERFACE_INCLUDE_DIRECTORIES "${PortsofCall_INCLUDE_DIR}")
#     list(APPEND TPL_LIBRARIES PortsofCall::PortsofCall)
# else()
#     message(FATAL_ERROR "For some odd reason, couldn't get Ports-of-Call")
# endif()


