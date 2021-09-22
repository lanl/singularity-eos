#------------------------------------------------------------------------------#
#  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
# /@@/////  /@@          @@////@@ @@////// /@@
# /@@       /@@  @@@@@  @@    // /@@       /@@
# /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
# /@@////   /@@/@@@@@@@/@@       ////////@@/@@
# /@@       /@@/@@//// //@@    @@       /@@/@@
# /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
# //       ///  //////   //////  ////////  //
#
# Copyright (c) 2016, Triad National Security, LLC
# All rights reserved
#------------------------------------------------------------------------------#

include(subdirlist)

if(NOT DEFINED LIBDIR)
  include(GNUInstallDirs)
  set(LIBDIR "${CMAKE_INSTALL_LIBDIR}")
endif(NOT DEFINED LIBDIR)

option(BUILD_SHARED_LIBS "Build shared libs" ON)
mark_as_advanced(BUILD_SHARED_LIBS)

function(add_library_target target directory)

  #----------------------------------------------------------------------------#
  # Setup argument options
  #----------------------------------------------------------------------------#

  set(options)
  set(one_value_args
    VERSION
    SOVERSION
    EXPORT_TARGET
    NAMESPACE
  )
  set(multi_value_args
    DEFINE_PUBLIC
    DEFINE_PRIVATE
    INCLUDE_PUBLIC
    INCLUDE_PRIVATE
    LINK_PUBLIC
    LINK_PRIVATE
  )

  cmake_parse_arguments(lib "${options}" "${one_value_args}"
    "${multi_value_args}" ${ARGN})

  #----------------------------------------------------------------------------#
  # Handle unspecified optional arguments
  #----------------------------------------------------------------------------#

  if(lib_VERSION AND NOT lib_SOVERSION)
    string(REPLACE "." ";" _VLIST ${lib_VERSION})
    list(GET _VLIST 0 _SOVERSION)
    set(lib_SOVERSION ${_SOVERSION})
  endif()

  #----------------------------------------------------------------------------#
  # Add target to list
  #----------------------------------------------------------------------------#

  message(STATUS
    "Adding library target ${target} with source directory ${directory}")

  #----------------------------------------------------------------------------#
  # Add any top-level source or headers files
  #----------------------------------------------------------------------------#

  set(_SOURCEDIR ${CMAKE_CURRENT_SOURCE_DIR}/${directory})

  if(EXISTS ${_SOURCEDIR}/library.cmake)
    include(${_SOURCEDIR}/library.cmake)
  endif()

  foreach(_HEADER ${${directory}_HEADERS})
    if(NOT EXISTS ${_SOURCEDIR}/${_HEADER})
      message(FATAL_ERROR
        "Header '${_HEADER}' from ${directory}_HEADERS does not exist.")
    endif()

    list(APPEND _INSTALL_HEADERS ${_HEADER})
    list(APPEND _HEADERS ${_SOURCEDIR}/${_HEADER})
  endforeach()

  foreach(_SOURCE ${${directory}_SOURCES})
    list(APPEND _SOURCES ${_SOURCEDIR}/${_SOURCE})
  endforeach()

  #----------------------------------------------------------------------------#
  # Add subdirectories
  #
  # This uses a glob, i.e., all sub-directories will be added at this level.
  # This is not true for levels below this one.  This allows some flexibility
  # while keeping the generic case as simple as possible.
  #----------------------------------------------------------------------------#

  make_subdirlist(_SUBDIRECTORIES ${_SOURCEDIR} False)

  #----------------------------------------------------------------------------#
  # Add subdirectory files
  #
  # This loop adds header and source files for each listed sub-directory
  # to the main header and source file lists.
  #----------------------------------------------------------------------------#

  foreach(_SUBDIR ${_SUBDIRECTORIES})

    if(NOT EXISTS ${_SOURCEDIR}/${_SUBDIR}/CMakeLists.txt)
      continue()
    endif()

    message(STATUS "Adding source subdirectory '${_SUBDIR}' to ${target}")

    unset(${_SUBDIR}_HEADERS)
    unset(${_SUBDIR}_SOURCES)

    add_subdirectory(${directory}/${_SUBDIR})
    list(APPEND _SUBDIRS ${_SOURCEDIR}/${_SUBDIR})

    foreach(_HEADER ${${_SUBDIR}_HEADERS})
      if(NOT EXISTS ${_SOURCEDIR}/${_SUBDIR}/${_HEADER})
        message(FATAL_ERROR
          "Header '${_HEADER}' from ${_SUBDIR}_HEADERS does not exist.")
      endif()

      list(APPEND _INSTALL_HEADERS ${_SUBDIR}/${_HEADER})
      list(APPEND _HEADERS ${_SOURCEDIR}/${_SUBDIR}/${_HEADER})
    endforeach()

    foreach(_SOURCE ${${_SUBDIR}_SOURCES})
      list(APPEND _SOURCES ${_SOURCEDIR}/${_SUBDIR}/${_SOURCE})
    endforeach()

  endforeach(_SUBDIR)

  #----------------------------------------------------------------------------#
  # Add the actual build target
  #----------------------------------------------------------------------------#

  if(_SOURCES)
    add_library(${target} ${_SOURCES} ${_HEADERS})
  else()
    add_library(${target} INTERFACE)
  endif()

  #----------------------------------------------------------------------------#
  # Create an alias for local builds
  #----------------------------------------------------------------------------#

  add_library(${lib_NAMESPACE}::${target} ALIAS ${target})

  #----------------------------------------------------------------------------#
  # Add compile defines
  #----------------------------------------------------------------------------#

  target_compile_definitions(${target}
    PUBLIC
      ${lib_DEFINE_PUBLIC}
    PRIVATE
      ${lib_DEFINE_PRIVATE}
  )

  #----------------------------------------------------------------------------#
  # Add the correct includes
  #----------------------------------------------------------------------------#

  target_include_directories(${target}
    PUBLIC
      $<INSTALL_INTERFACE:include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )

  target_include_directories(${target}
    SYSTEM
    PUBLIC
      ${lib_INCLUDE_PUBLIC}
    PRIVATE
      ${lib_INCLUDE_PRIVATE}
  )

  #----------------------------------------------------------------------------#
  # Add link libraries
  #----------------------------------------------------------------------------#

  target_link_libraries(${target}
    PUBLIC
      ${lib_LINK_PUBLIC}
    PRIVATE
      ${lib_LINK_PRIVATE}
  )

  #----------------------------------------------------------------------------#
  # Install
  #----------------------------------------------------------------------------#

  if(lib_EXPORT_TARGET)
    install(TARGETS ${target}
      EXPORT
        ${lib_EXPORT_TARGET}
      DESTINATION
        ${LIBDIR}
    )

    install(EXPORT ${lib_EXPORT_TARGET}
      FILE
        ${lib_EXPORT_TARGET}.cmake
      NAMESPACE
        "${lib_NAMESPACE}::"
      DESTINATION
        ${LIBDIR}/cmake/${target}
    )
  else()
    install(TARGETS ${target} DESTINATION ${LIBDIR})
  endif()

  #----------------------------------------------------------------------------#
  # Header install
  #----------------------------------------------------------------------------#

  foreach(file ${_INSTALL_HEADERS})
    get_filename_component(DIR ${file} DIRECTORY)
    install(FILES ${directory}/${file}
      DESTINATION include/${directory}/${DIR})
  endforeach()

  #----------------------------------------------------------------------------#
  # Version
  #----------------------------------------------------------------------------#

  if(lib_VERSION)
    set_target_properties(${target}
      PROPERTIES
        VERSION ${lib_VERSION}
        SOVERSION ${lib_SOVERSION}
    )

    include(CMakePackageConfigHelpers)
    write_basic_package_version_file(
      ${CMAKE_CURRENT_BINARY_DIR}/${target}ConfigVersion.cmake
      VERSION ${lib_VERSION}
      COMPATIBILITY AnyNewerVersion
    )
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${target}ConfigVersion.cmake
      DESTINATION
        ${LIBDIR}/cmake/${target}
    )
  endif()

  #----------------------------------------------------------------------------#
  # Local export
  #----------------------------------------------------------------------------#

  export(EXPORT ${lib_EXPORT_TARGET}
    FILE
      ${CMAKE_CURRENT_BINARY_DIR}/${target}Targets.cmake
    NAMESPACE
      ${lib_NAMESPACE}::
  )

  export(PACKAGE ${target})

endfunction()
