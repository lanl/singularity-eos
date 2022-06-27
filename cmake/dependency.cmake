#------------------------------------------------------------------------------#
# © 2021-2022. Triad National Security, LLC. All rights reserved.  This
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

#------------------------------------------------------------------------------
# including manually written submodule configs
#------------------------------------------------------------------------------
include(cmake/submodule_configs.cmake)

#------------------------------------------------------------------------------
# singularity_import_dependency
#------------------------------------------------------------------------------
#   optional arguments:
#     - SUBMODULE_ONLY:     only use the submodule, do not use `find_package` 
#                           or other
#   keyword arguments, single-value:
#     - PKG:                the name of the package (e.g. `find_package()`)
#     - SUBDIR (optional):  directory of submodule where to find package 
#                           source tree
#     - TARGET:             the target to import(/export later mayber).
#                           CANNOT be used with `TARGETS`
#   keyword arguments, multi-values:
#     - TARGETS:            list of targets to import(/export later mayber).
#                           CANNOT be used with `TARGET`
#     - COMPONENTS:         the components of PKG to require
#
#   the function takes a target/list of targets and package names from the caller, 
#   and trys to find the package, either using `find_package()` or using 
#   `add_subdirectory()`.
#
#   users may ask for `find_package()` to search a particular directory path
#   by setting the CMake variable `SINGULARITY_${PKG_UPPERCASE}_INSTALL_DIR` 
#   prior to configuration. By default, only some packages are given then option 
#   (see ${PROJECT_SOURCE_DIR}/CMakeList.txt). For example, to use the `Kokkos` 
#   package:
#   
#     $> cmake -DSINGULARITY_USE_KOKKOS=ON \
#       -DSINGULARITY_KOKKOS_INSTALL_DIR=<...>
#
#   If CMake cannot find the expected configuration files in the provided path,
#   an error will be generated and processing will stop.
#
#   users may also specify a source directory to "attach" to the build tree with
#   `SINGULARITY_${PKG_UPPERCASE}_IMPORT_DIR`. This is not recommended.
#
#     $> cmake -DSINGULARITY_USE_KOKKOSKERNELS=ON \
#       -DSINGULARITY_KOKKOSKERNELS_IMPORT_DIR=<...>
#
#   if both `SINGULARITY_${PKG_UPPERCASE}_INSTALL_DIR` and 
#   `SINGULARITY_${PKG_UPPERCASE}_IMPORT_DIR` are both given, an error will be
#   generated and processing will stop.
#
#   should `${PKG}` contain dashes ("-"), ${PKG_UPPERCASE} should use underscores
#   ("_"), e.g. for `ports-of-call` use `-DSINGULARITY_PORTS_OF_CALL_IMPORT_DIR=<...>`
#
#   some packages (e.g. `HDF5`, `MPI`) are are designed to provide `COMPONENTS` 
#   as seperate targets. These should always use the `TARGETS` and `COMPONENTS` 
#   keyword arguments in tandem, so that it is possible to export/install these
#   requirements later on.
#------------------------------------------------------------------------------
macro(singularity_import_dependency)
  set(options
    SUBMODULE_ONLY
  )
  set(one_value_args
    PKG
    TARGET
    SUBDIR
  )
  set(multi_value_args
    COMPONENTS
    TARGETS
  )
  cmake_parse_arguments(dep "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

  set(_SMSG_PREFIX "[IMPORT ${dep_PKG}]")

  # first, if this target is already defined, then ignore figuring out imports
  if(TARGET ${dep_TARGET})
    singularity_msg(STATUS "Detected ${dep_TARGET} present, no further processing for this dependency.")
  else()
    # check for user overrides
    # sets `SINGULARITY_IMPORT_USER_OVERRIDE_${dep_PKG}`
    # if that user requests import
    singularity_import_check_user_override(
      PKG ${dep_PKG}
      SUBDIR ${dep_SUBDIR}
    )

    # if only use submodule, or user requests import of external source
    if(dep_SUBMODULE_ONLY OR SINGULARITY_IMPORT_USER_OVERRIDE_${dep_PKG})
      singularity_import_submodule(
        PKG ${dep_PKG}
        SUBDIR ${dep_SUBDIR}
      )
    else()
      # the standard path of dependency importing.
      # check if pkg can be found in environment
      singularity_msg(STATUS "Attempting \"find_package(${dep_PKG})\"")
      find_package(${dep_PKG} QUIET COMPONENTS ${dep_COMPONENTS})

      # if we fail to find the package, try a submodule
      if (NOT ${dep_PKG}_FOUND)
        # does `add_subdirectory` with ${SUBDIR}
        singularity_import_submodule(
          PKG ${dep_PKG}
          SUBDIR ${dep_SUBDIR}
        )
      else()
        singularity_msg(STATUS "Found with find_package() [${${dep_PKG}_DIR}]")
      endif() # NOT FOUND
    endif() # SUBMODULE ONLY
  endif() # TARGET

  # if we made it hear, there should be valid imports available.
  # we record these to the global scope for later use and processing
  # TODO: split out public/private libs
  list(APPEND SINGULARITY_PUBLIC_LIBS ${dep_TARGETS})
  list(APPEND SINGULARITY_PUBLIC_LIBS ${dep_TARGET})
  list(APPEND SINGULARITY_DEP_PKGS ${dep_PKG})
  if(dep_COMPONENTS)
    list(APPEND SINGULARITY_DEP_PKGS_${dep_PKG} "${dep_COMPONENTS}")
  endif()
  unset(_SMSG_PREFIX)
endmacro() # singularity_import_dependency

#------------------------------------------------------------------------------
# Helper functions
#
# These 'cheat' scoping by assuming ${dep_*} is available, which can be 
# done becuase `macro` doesn't introduce a new scope, it just "injects"
# the body of the macro code at the callsite. As a result, they only can 
# be validly invoked within `singularity_import_dependency`.
# I do this for simplicity, passing around variables references is
# already a pain in CMake.
#------------------------------------------------------------------------------

macro(singularity_import_check_user_override)
  string(TOUPPER ${dep_PKG} dep_CAP)
  string(REPLACE "-" "_" dep_CAP "${dep_CAP}")

  if(SINGULARITY_${dep_CAP}_INSTALL_DIR AND SINGULARITY_${dep_CAP}_IMPORT_DIR)
    singularity_msg(FATAL_ERROR "Cannot set INSTALL_DIR and IMPORT_DIR for the same package. You may have a stale cache varaible.")
  endif()

  if(SINGULARITY_${dep_CAP}_INSTALL_DIR)
    find_path(${dep_PKG}_ROOT 
      NAMES "${dep_PKG}Config.cmake"
      PATHS "${SINGULARITY_${dep_CAP}_INSTALL_DIR}"
      PATH_SUFFIXES lib/cmake/${dep_PKG} lib64/cmake/${dep_PKG}
      NO_DEFAULT_PATH
    )
    if(NOT ${dep_PKG}_ROOT)
      singularity_msg(FATAL_ERROR "SINGULARITY_${dep_VARCASE}_INSTALL_DIR [${SINGULARITY_${dep_CAP}_INSTALL_DIR}] set, but did not find an installed CMake package. Use a valid install path or unset this variable.")
    else()
      singularity_msg(STATUS "located user-provided install, ${dep_PKG}_ROOT = ${${dep_PKG}_ROOT}")
    endif()
  endif()
  if(SINGULARITY_${dep_CAP}_IMPORT_DIR)
    singularity_msg(STATUS "SINGULARITY_${dep_CAP}_IMPORT_DIR is set. I will assume you know what you are doing. Overriding default path [${dep_SUBDIR}] to ${SINGULARITY_${dep_CAP}_IMPORT_DIR} and skipping find logic")
    set(dep_SUBDIR "${SINGULARITY_${dep_CAP}_IMPORT_DIR}")
    set(SINGULARITY_IMPORT_USER_OVERRIDE_${dep_PKG} ON)
  endif()
  
endmacro()

macro(singularity_import_submodule)
  if (NOT EXISTS ${dep_SUBDIR})      
    singularity_msg(FATAL_ERROR "${dep_PKG} requested, but the in-tree directory \"${dep_SUBDIR}\" does not exist")
  endif()
      
  if (NOT EXISTS ${dep_SUBDIR}/CMakeLists.txt)
    singularity_msg(WARNING "submodendency directory ${dep_SUBDIR} does not contain a CMakeLists.txt file. No target information about ${dep_PKG} will be available")
  endif()
      
  # if adding a subdirectory, set the config variables (if any) for the package
  singularity_cmake_config(${dep_PKG})
      
  add_subdirectory(${dep_SUBDIR} "${CMAKE_CURRENT_BINARY_DIR}/extern/${dep_PKG}")

  singularity_msg(STATUS "${dep_PKG} added from: ${dep_SUBDIR}")
endmacro()
