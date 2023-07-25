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

# Find the native EOSPAC headers and libraries.
#
#  EOSPAC_INCLUDE_DIRS - where to find eos_Interface.h, etc.
#  EOSPAC_LIBRARIES    - List of libraries when using eospac6.
#  EOSPAC_FOUND        - True if eospac found.

#TODO: Add EOSPAC_MODULES (possibly part of include dirs) and clarify which interface we need

# Look for the header file.
FIND_PATH(EOSPAC_INCLUDE_DIR NAMES eos_Interface.h)

# Look for the library.
FIND_LIBRARY(EOSPAC_LIBRARY NAMES eospac6 libeospac6 libeospac6gpu.a)

# handle the QUIETLY and REQUIRED arguments and set EOSPAC_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(EOSPAC DEFAULT_MSG EOSPAC_LIBRARY EOSPAC_INCLUDE_DIR)

# Copy the results to the output variables.
SET(EOSPAC_LIBRARIES ${EOSPAC_LIBRARY})
SET(EOSPAC_INCLUDE_DIRS ${EOSPAC_INCLUDE_DIR})

MARK_AS_ADVANCED(EOSPAC_INCLUDE_DIR EOSPAC_LIBRARY)

if(EOSPAC_INCLUDE_DIR AND EOSPAC_LIBRARY)
	add_library(EOSPAC::eospac UNKNOWN IMPORTED)
	set_target_properties(EOSPAC::eospac PROPERTIES
		IMPORTED_LOCATION "${EOSPAC_LIBRARY}"
		INTERFACE_INCLUDE_DIRECTORIES "${EOSPAC_INCLUDE_DIR}")
endif()
