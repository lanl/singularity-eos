#------------------------------------------------------------------------------
# placeholder
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# including manually written submodule configs
#------------------------------------------------------------------------------
include(cmake/submodule_configs.cmake)

#------------------------------------------------------------------------------
# append_target_dependency
#------------------------------------------------------------------------------
#   positional arguments:
#     - targetlist(in/out): a list of targets that have been added by this function
#     - deplist(in/out): a list of dependencies this function has added
#   keyword arguments:
#     - PKG: the name of the package (e.g. `find_package()`)
#     - SUBDIR (optional): directory of submodule where to find package source tree
#     - TARGETS: list of targets to use from PKG
#     - COMPONENTS (optional): the components of PKG to require
#
#   the function takes a list of targets and package names from the caller, 
#   and trys to find the package, either using `find_package()` or using 
#   `add_subdirectory()`. The targets and package names are updated and passed
#   back to the caller. if any components were requested, these are returned as 
#   the list `${deplist}_${dep_PKG}`. for example:
#     append_target_dependency(TARS, DEPS
#       PKG MPI
#       TARGETS MPI::MPI MPI::CXX
#       COMPONENTS C CXX
#     )
#     # message("${DEPS_MPI}") will print "C;CXX"
#------------------------------------------------------------------------------
macro(append_target_dependency targetlist deplist)
  set(options)
  set(one_value_args
    PKG
    SUBDIR
  )
  set(multi_value_args
    TARGETS
    COMPONENTS
  )
  cmake_parse_arguments(dep "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

  # this is just temporary
  if(SINGULARITY_${dep_PKG}_INSTALL_DIR)
    list(APPEND CMAKE_PREFIX_PATH "${SINGULARITY_${dep_PKG}_INSTALL_DIR}")
  endif()

  if(NOT dep_SUBDIR)
    set(dep_SUBDIR ${PROJECT_SOURCE_DIR}/utils/${dep_PKG})
  endif()

  # TODO: this is a bit awkward, we want to check if target is already defined,
  # but in the general case a package may define multiple targets. currently,
  # we will bail if ANY of the targets we request already exist; that doesn't
  # feel robust, but it covers enough inputs for this project to not be an issue.
  set(_anyTargetsDefined 0)
  foreach(_targ ${dep_TARGETS})
    if(TARGET ${_targ})
      set(_anyTargetsDefined 1)
      break()
    endif()
  endforeach()

  if(NOT _anyTargetsDefined)
    find_package(${dep_PKG} QUIET COMPONENTS ${dep_COMPONENTS})
    if (NOT ${dep_PKG}_FOUND)
      if (NOT EXISTS ${dep_SUBDIR})
        message(FATAL_ERROR "dependency ${dep_PKG} requested, but cannot locate in system and the provided directory ${dep_SUBDIR} does not exist")
      endif()
      if (NOT EXISTS ${dep_SUBDIR}/CMakeLists.txt)
        message(WARNING "dependency directory ${dep_SUBDIR} does not contain a CMakeLists.txt file. No target information about ${dep_PKG} will be available")
      endif()
      # if adding a subdirectory, set the config variables (if any) for the package
      singularity_cmake_config(${dep_PKG})
      message(STATUS "invoking \"add_subdirectory(${dep_SUBDIR}\", output supressed")
      set(MESSAGE_QUIET ON)
      add_subdirectory(${dep_SUBDIR})
      unset(MESSAGE_QUIET)
      message(STATUS "dependency ${dep_PKG} added from in-tree: ${dep_SUBDIR}")
    else()
      message(STATUS "dependency ${dep_PKG} located: ${${dep_PKG}_DIR}")
    endif()
  endif()

  list(APPEND ${targetlist} ${dep_TARGETS})
  list(APPEND ${deplist} ${dep_PKG})

  if(dep_COMPONENTS)
    set(${deplist}_${dep_PKG} "${dep_COMPONENTS}")
  endif()

endmacro()

