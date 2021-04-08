#------------------------------------------------------------------------------#
# © 2021. Triad National Security, LLC. All rights reserved.  This
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

find_library(HWLOC_LIBRARY NAMES hwloc)
find_path(HWLOC_INCLUDE_DIR hwloc.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Hwloc DEFAULT_MSG HWLOC_LIBRARY HWLOC_INCLUDE_DIR)

add_library(Hwloc::hwloc UNKNOWN IMPORTED)
set_target_properties(Hwloc::hwloc PROPERTIES
    IMPORTED_LOCATION ${HWLOC_LIBRARY}
    INTERFACE_INCLUDE_DIRECTORIES ${HWLOC_INCLUDE_DIR})
