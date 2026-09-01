# ------------------------------------------------------------------------------#
# © 2026. Triad National Security, LLC. All rights reserved.  This program
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

macro(singularity_import_pybind11)
  if(NOT TARGET pybind11::headers)
    message(STATUS "Using pybind11 submodule")
    add_subdirectory(utils/pybind11)
  endif()
endmacro()

macro(singularity_find_pybind11)
  message(STATUS "Searching for system pybind11")
  find_package(pybind11 REQUIRED)
endmacro()
