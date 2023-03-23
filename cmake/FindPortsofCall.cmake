#------------------------------------------------------------------------------#
# © 2021-2023. Triad National Security, LLC. All rights reserved.  This
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

# Ports of Call
SET(PortsofCall_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/utils/ports-of-call")

# handle the QUIETLY and REQUIRED arguments and set EOSPAC_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PortsofCall DEFAULT_MSG
                                  )#PortsofCall_INCLUDE_DIR)

# Copy the results to the output variables.
SET(PortsofCall_INCLUDE_DIRS ${PortsofCall_INCLUDE_DIR})

add_library(PortsofCall::PortsofCall INTERFACE IMPORTED)
set_target_properties(PortsofCall::PortsofCall PROPERTIES
		      INTERFACE_INCLUDE_DIRECTORIES "${PortsofCall_INCLUDE_DIR}")
