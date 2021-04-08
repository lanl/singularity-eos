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

# Eigen3
SET(Eigen3_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/utils/eigen")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Eigen3 DEFAULT_MSG)

# Copy the results to the output variables.
SET(Eigen3_INCLUDE_DIRS ${Eigen3_INCLUDE_DIR})

add_library(Eigen3::Eigen INTERFACE IMPORTED)
set_target_properties(Eigen3::Eigen PROPERTIES
		      INTERFACE_INCLUDE_DIRECTORIES "${Eigen3_INCLUDE_DIR}")
