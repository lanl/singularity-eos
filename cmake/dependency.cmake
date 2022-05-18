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
macro(singularity_import_dependency)
  set(options
    USER_INSTALL
    SUBMODULE_ONLY
    NO_SUBMODULE
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

  if(TARGET ${dep_TARGET})
    singularity_msg(STATUS "Detected ${dep_TARGET} present, no further processing for this dependency.")
  else()
    if(dep_SUBMODULE_ONLY)
      singularity_import_submodule(
        PKG ${dep_PKG}
        SUBDIR ${dep_SUBDIR}
      )
    else()
      if(dep_USER_INSTALL)
        singularity_import_user_install(
          PKG ${dep_PKG}
        )
      endif()
      singularity_import_system(
        PKG ${dep_PKG}
        TARGETS ${dep_TARGETS}
        COMPONENTS ${dep_COMPONENTS}
      )
      if (NOT ${dep_PKG}_FOUND)
        if(dep_NO_SUBMODULE)
          message(FATAL_ERROR "Could not locate ${dep_PKG} outside of project, and `singularity_import_dependency()` was called with `NO_SUBMODULE` option")
        endif()
        singularity_import_submodule(
          PKG ${dep_PKG}
          SUBDIR ${dep_SUBDIR}
        )
      else()
        singularity_msg(STATUS "Found with find_package() [${${dep_PKG}_DIR}]")
      endif()
    endif()
  endif()

  list(APPEND SINGULARITY_PUBLIC_LIBS ${dep_TARGETS})
  list(APPEND SINGULARITY_PUBLIC_LIBS ${dep_TARGET})

  if(dep_TARGET)
    #print_target_properties(${dep_TARGET})
  endif()
  list(APPEND SINGULARITY_DEP_PKGS ${dep_PKG})
  if(dep_COMPONENTS)
    list(APPEND SINGULARITY_DEP_PKGS_${dep_PKG} "${dep_COMPONENTS}")
  endif()

endmacro()

macro(singularity_import_user_install)
  set(options
  )
  set(one_value_args
    PKG
  )
  set(multi_value_args
  )
  cmake_parse_arguments(ui "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

  string(TOUPPER ${ui_PKG} ui_VARCASE)

  if(SINGULARITY_${ui_VARCASE}_INSTALL_DIR)
    find_path(ui_INSTALLCMAKE 
      NAMES "${ui_PKG}Config.cmake"
      PATHS "${SINGULARITY_${ui_VARCASE}_INSTALL_DIR}"
    )
    if(ui_INSTALLCMAKE-NOTFOUND)
      message(FATAL_ERROR "Could not find \"${ui_PKG}Config.cmake\" in \"SINGULARITY_${ui_VARCASE}_INSTALL_DIR}")
    endif()
    set(${ui_PKG}_ROOT "${ui_INSTALLCMAKE}")
  endif()
  
endmacro()

macro(singularity_import_submodule)
  set(options
  )
  set(one_value_args
    PKG
    SUBDIR
  )
  set(multi_value_args
  )
  cmake_parse_arguments(submod "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

  if (NOT EXISTS ${submod_SUBDIR})      
    message(FATAL_ERROR "${submod_PKG} requested, but cannot the in-tree directory ${submod_SUBDIR} does not exist")
  endif()
      
  if (NOT EXISTS ${submod_SUBDIR}/CMakeLists.txt)
    singularity_msg(WARNING "submodendency directory ${submod_SUBDIR} does not contain a CMakeLists.txt file. No target information about ${submod_PKG} will be available")
  endif()
      
  # if adding a subdirectory, set the config variables (if any) for the package
  singularity_cmake_config(${submod_PKG})
      
  singularity_msg(STATUS "invoking \"add_subdirectory(${submod_SUBDIR}\", output supressed")
  #  set(MESSAGE_QUIET ON)
  add_subdirectory(${submod_SUBDIR})
  #unset(MESSAGE_QUIET)

  singularity_msg(STATUS "${submod_PKG} added from in-tree: ${submod_SUBDIR}")
endmacro()

macro(singularity_import_system)
  set(options
  )
  set(one_value_args
    PKG
  )
  set(multi_value_args
    COMPONENTS
  )
  cmake_parse_arguments(sysinstall "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

    singularity_msg(STATUS "Attempting \"find_package(${sysinstall_PKG})\"")
    find_package(${sysinstall_PKG} QUIET COMPONENTS ${sysinstall_COMPONENTS})
endmacro()

