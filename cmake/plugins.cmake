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

set(EOS_HEADERS "")
set(_install_headers "")
function(register_headers)
  foreach(arg IN LISTS ARGN)
    file(RELATIVE_PATH relative_path ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
    list(APPEND EOS_HEADERS ${relative_path}/${arg})
    list(APPEND _install_headers ${arg})
  endforeach()
  set(EOS_HEADERS ${EOS_HEADERS} PARENT_SCOPE)
  set(_install_headers ${_install_headers} PARENT_SCOPE)
endfunction()

set(EOS_SRCS "")
function(register_srcs)
  foreach(arg IN LISTS ARGN)
    file(RELATIVE_PATH relative_path ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
    list(APPEND EOS_SRCS ${relative_path}/${arg})
  endforeach()
  set(EOS_SRCS ${EOS_SRCS} PARENT_SCOPE)
endfunction()

set(PLUGIN_TESTS "")
function(register_tests)
  foreach(arg IN LISTS ARGN)
    list(APPEND PLUGIN_TESTS ${CMAKE_CURRENT_SOURCE_DIR}/${arg})
  endforeach()
  set(PLUGIN_TESTS ${PLUGIN_TESTS} PARENT_SCOPE)
endfunction()

macro(export_plugin)
  set(EOS_HEADERS ${EOS_HEADERS} PARENT_SCOPE)
  set(_install_headers ${_install_headers} PARENT_SCOPE)
  set(EOS_SRCS ${EOS_SRCS} PARENT_SCOPE)
  set(PLUGIN_TESTS ${PLUGIN_TESTS} PARENT_SCOPE)
endmacro()
