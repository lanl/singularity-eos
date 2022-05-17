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
macro(import_dependency)
  set(options
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

  if(dep_TARGET AND dep_TARGETS)
    message(FATAL_ERROR "Can only import a single target, or multiple")
  endif()

  if(dep_TARGETS AND NOT dep_COMPONENTS)
    message(FATAL_ERROR "Need a target for each component")
  endif()

  if(TARGET ${dep_TARGET})
    singularity_msg(STATUS "Detected ${dep_TARGET} present, no further processing for this dependency.")
  else()
    singularity_msg(STATUS "Attempting \"find_package(${dep_PKG})\"")
    find_package(${dep_PKG} QUIET COMPONENTS ${dep_COMPONENTS})
    
    if (NOT ${dep_PKG}_FOUND)
    
      singularity_msg(STATUS "\"find_package(${dep_PKG})\" returned not-found. Trying to locate submodule at")
      
      if (NOT EXISTS ${dep_SUBDIR})
      
        message(FATAL_ERROR "${dep_PKG} requested, but cannot the in-tree directory ${dep_SUBDIR} does not exist")
      endif()
      
      if (NOT EXISTS ${dep_SUBDIR}/CMakeLists.txt)
      
        singularity_msg(WARNING "dependency directory ${dep_SUBDIR} does not contain a CMakeLists.txt file. No target information about ${dep_PKG} will be available")
      endif()
      
      # if adding a subdirectory, set the config variables (if any) for the package
      singularity_cmake_config(${dep_PKG})
      
      singularity_msg(STATUS "invoking \"add_subdirectory(${dep_SUBDIR}\", output supressed")
      set(MESSAGE_QUIET ON)
      add_subdirectory(${dep_SUBDIR})
      unset(MESSAGE_QUIET)

      #      get_target_property(_is_imported ${dep_TARGET} IMPORTED)
      #      set(_lib_type "INTERFACE")
      #      if(_is_imported)
      #        set(_lib_type "IMPORTED")
      #      endif()

      singularity_msg(STATUS "dependency ${dep_PKG} added from in-tree: ${dep_SUBDIR}")

      # anything with `add_subdirectory`...umm interface lib

      # get de-aliased target
      get_target_property(_taraliased ${dep_TARGET} ALIASED_TARGET)
      # link libs added to export
      get_target_property(_intflink ${dep_TARGET} INTERFACE_LINK_LIBRARIES)

      if(_taraliased)
        list(APPEND SINGULARITY_EXPORT_TARGETS ${_taraliased})
      else()
        list(APPEND SINGULARITY_EXPORT_TARGETS ${dep_TARGET})
      endif()
      if(_intflink)
        list(APPEND SINGULARITY_EXPORT_TARGETS ${_intflink})
      endif()

    else()
      singularity_msg(STATUS "dependency ${dep_PKG} located: ${${dep_PKG}_DIR}")
    endif()
  endif()

  if(dep_TARGETS)
    list(APPEND SINGULARITY_PUBLIC_LIBS ${dep_TARGETS})
  else()
    list(APPEND SINGULARITY_PUBLIC_LIBS ${dep_TARGET})
  endif()

  list(APPEND SINGULARITY_DEP_PKGS ${dep_PKG})
  if(dep_COMPONENTS)
    list(APPEND SINGULARITY_DEP_PKGS_${dep_PKG} ${${dep_COMPONETNS}})
  endif()

endmacro()

