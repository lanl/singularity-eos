

function(add_target_dependency target)
  #argument options

  set(options)
  set(one_value_args
    PKG
    NAMESPACE
    SUBDIR
  )

  set(multi_value_args
    PKG_CACHE_OPTS
  )

  cmake_parse_arguments(dep "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

  if (NOT TARGET ${dep_NAMESPACE}::${dep_PKG})
    find_package(${dep_PKG} QUIET)
    if (NOT ${dep_PKG}_FOUND)
      foreach(_OPT ${dep_PKG_CACHE_OPTS})
        set(${_OPT} OFF CACHE BOOL "" FORCE)
      endforeach()
      add_subdirectory(${dep_SUBDIR})
    endif()
  endif()
  target_link_libraries(${target} PUBLIC ${dep_NAMESPACE}::${dep_PKG})

endfunction()

