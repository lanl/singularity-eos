#------------------------------------------------------------------------------
# placeholder
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# including manually written submodule configs
#------------------------------------------------------------------------------
include(cmake/submodule_configs.cmake)

#------------------------------------------------------------------------------
# singularity_import_dependency
#------------------------------------------------------------------------------
#   optional arguments:
#     - USER_INSTALL:       look for a cmake configuration file located in
#                           `SINGULARITY_${PKG_UPPERCASE}_INSTALL_DIR`. The package
#                           search will continue if not found.
#     - SUBMODULE_ONLY:     only use the submodule, do not use `find_package` 
#                           or other
#     - NO_SUBMODULE:       run the different package import modes, except 
#                           for looking for a submodule
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
#     $> cmake -DSINGULARITY_USE_KOKKOS=ON -DSINGULARITY_KOKKOS_INSTALL_DIR=<...>
#
#   NB: this does procedure will NOT fail if `Kokkos` is not found in the provided 
#   directory, and will continue to the next step to search either a submodule or
#   through `find_package()`, depending on the other options provided in the call.
#
#   some packages (e.g. `HDF5`, `MPI`) are are designed to provide `COMPONENTS` 
#   as seperate targets. These should always use the `TARGETS` and `COMPONENTS` 
#   keyword arguments in tandem, so that it is possible to export/install these
#   requirements later on.
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

  # first, if this target is already defined, then ignore figuring out imports
  if(TARGET ${dep_TARGET})
    singularity_msg(STATUS "[IMPORT] Detected ${dep_TARGET} present, no further processing for this dependency.")
  else()
    # only use a submodule provided by `SUBDIR`
    if(dep_SUBMODULE_ONLY)
      singularity_import_submodule(
        PKG ${dep_PKG}
        SUBDIR ${dep_SUBDIR}
      )
    else()
      # if option set, try to find the a pre-existing configuration of `PKG`
      # NB: if successful, then `${PKG}_ROOT` is set, and `find_package()` will
      #     first pick up the configuration in `${PKG}_ROOT`.
      if(dep_USER_INSTALL)
        singularity_import_user_install(
          PKG ${dep_PKG}
        )
      endif() # dep_USER_INSTALL

      # proceed with simple `find_package()` search
      singularity_import_system(
        PKG ${dep_PKG}
        TARGETS ${dep_TARGETS}
        COMPONENTS ${dep_COMPONENTS}
      )
      # if we fail to find the package, try a submodule
      if (NOT ${dep_PKG}_FOUND)
        # if `NO_SUBMODULE` set, emit an error and stop
        if(dep_NO_SUBMODULE)
          message(FATAL_ERROR "[IMPORT] Could not locate ${dep_PKG} outside of project, and `singularity_import_dependency()` was called with `NO_SUBMODULE` option")
        endif()
        # does `add_subdirectory` with ${SUBDIR}
        singularity_import_submodule(
          PKG ${dep_PKG}
          SUBDIR ${dep_SUBDIR}
        )
      else()
        singularity_msg(STATUS "[IMPORT] Found with find_package() [${${dep_PKG}_DIR}]")
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
endmacro() # singularity_import_dependency

#------------------------------------------------------------------------------
# Helper functions
#------------------------------------------------------------------------------

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
    find_path(${ui_PKG}_ROOT 
      NAMES "${ui_PKG}Config.cmake"
      PATHS "${SINGULARITY_${ui_VARCASE}_INSTALL_DIR}"
      PATH_SUFFIXES lib/cmake/${ui_PKG} lib64/cmake/${ui_PKG}
    )
    if(NOT ${ui_PKG}_ROOT)
      singularity_msg(WARNING "[IMPORT:USER] SINGULARITY_${ui_PKG}_INSTALL_DIR [${SINGULARITY_${ui_PKG}_INSTALL_DIR}] set, but did not find \"${ui_PKG}Config.cmake\"")
    else()
      singularity_msg(STATUS "[IMPORT:USER] located cmake install, ${ui_PKG}_ROOT = ${${ui_PKG}_ROOT}")
    endif()
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
    singularity_msg(WARNING "[IMPORT:SUBMODULE] submodendency directory ${submod_SUBDIR} does not contain a CMakeLists.txt file. No target information about ${submod_PKG} will be available")
  endif()
      
  # if adding a subdirectory, set the config variables (if any) for the package
  singularity_cmake_config(${submod_PKG})
      
  #set(MESSAGE_QUIET ON)
  add_subdirectory(${submod_SUBDIR})
  #unset(MESSAGE_QUIET)

  singularity_msg(STATUS "[IMPORT:SUBMODULE] ${submod_PKG} added from in-tree: ${submod_SUBDIR}")
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

    singularity_msg(STATUS "[IMPORT:SYSTEM] Attempting \"find_package(${sysinstall_PKG})\"")
    find_package(${sysinstall_PKG} QUIET COMPONENTS ${sysinstall_COMPONENTS})
endmacro()

