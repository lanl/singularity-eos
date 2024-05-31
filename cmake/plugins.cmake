#------------------------------------------------------------------------------#
# © 2024. Triad National Security, LLC. All rights reserved.  This
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
#------------------------------------------------------------------------------#

set_property(GLOBAL PROPERTY EOS_SRCS "")
set_property(GLOBAL PROPERTY EOS_HEADERS "")
set_property(GLOBAL PROPERTY EOS_INSTALL_HEADERS "")
set_property(GLOBAL PROPERTY PLUGIN_TESTS "")
set_property(GLOBAL PROPERTY PLUGIN_INCLUDE_PATHS "")

function(register_headers)
  set(keyword_args PLUGIN)
  cmake_parse_arguments(ARG "" "${keyword_args}" "" ${ARGN})
  set(variadic_args ${ARG_UNPARSED_ARGUMENTS})

  get_property(eos_headers GLOBAL PROPERTY EOS_HEADERS)
  get_property(install_headers GLOBAL PROPERTY EOS_INSTALL_HEADERS)

  foreach(arg IN LISTS variadic_args)
    file(RELATIVE_PATH relative_path ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
    if (ARG_PLUGIN)
      list(APPEND eos_headers ${relative_path}/${ARG_PLUGIN}/${arg})
      list(APPEND install_headers ${ARG_PLUGIN}/${arg})
    else()
      list(APPEND eos_headers ${relative_path}/${arg})
      list(APPEND install_headers ${relative_path}/${arg})
    endif()
  endforeach()
  set_property(GLOBAL PROPERTY EOS_HEADERS "${eos_headers}")
  set_property(GLOBAL PROPERTY EOS_INSTALL_HEADERS "${install_headers}")
endfunction()

function(register_srcs)
  get_property(eos_srcs GLOBAL PROPERTY EOS_SRCS)
  foreach(arg IN LISTS ARGN)
    file(RELATIVE_PATH relative_path ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
    list(APPEND eos_srcs ${relative_path}/${arg})
  endforeach()
  set_property(GLOBAL PROPERTY EOS_SRCS "${eos_srcs}")
endfunction()

function(register_tests)
  get_property(plugin_tests GLOBAL PROPERTY PLUGIN_TESTS)
  foreach(arg IN LISTS ARGN)
    list(APPEND plugin_tests ${CMAKE_CURRENT_SOURCE_DIR}/${arg})
  endforeach()
  set_property(GLOBAL PROPERTY PLUGIN_TESTS "${plugin_tests}")
endfunction()

macro(export_plugin)
  get_property(plugin_include_paths GLOBAL PROPERTY PLUGIN_INCLUDE_PATHS)
  list(APPEND plugin_include_paths ${CMAKE_CURRENT_SOURCE_DIR})
  set_property(GLOBAL PROPERTY PLUGIN_INCLUDE_PATHS "${plugin_include_paths}")
endmacro()
