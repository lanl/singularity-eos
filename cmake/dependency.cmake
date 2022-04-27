
include(cmake/submodule_configs.cmake)

function(append_target_dependency targetlist deplist)

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

  if(NOT dep_SUBDIR)
    set(dep_SUBDIR ${PROJECT_SOURCE_DIR}/utils/${dep_PKG})
  endif()

  set(_anyTargetsDefined 0)
  foreach(_tar ${dep_TARGETS})
    if(TARGET ${_tar})
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
      singularity_cmake_config(${dep_PKG})
      add_subdirectory(${dep_SUBDIR})
    endif()
  endif()

  list(APPEND ${targetlist} ${dep_TARGETS})
  list(APPEND ${deplist} ${dep_PKG})

  # recover variables in caller-scope
  set(${targetlist} "${${targetlist}}" PARENT_SCOPE)
  set(${deplist} "${${deplist}}" PARENT_SCOPE)
  # if components were requested, recover these as well
  if(dep_COMPONENTS)
    set(${deplist}_${dep_PKG} "${dep_COMPONENTS}" PARENT_SCOPE)
  endif()

endfunction()

