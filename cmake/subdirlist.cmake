function(make_subdirlist result directory recursive)
  if(${recursive})
    file(GLOB_RECURSE _CHILDREN RELATIVE ${directory} ${directory}/*)
  else()
    file(GLOB _CHILDREN RELATIVE ${directory} ${directory}/*)
  endif()

  foreach(_CHILD ${_CHILDREN})
    if(NOT IS_DIRECTORY ${directory}/${_CHILD})
      get_filename_component(_DIR ${_CHILD} PATH)
    else()
      set(_DIR ${_CHILD})
    endif()

    list(APPEND _DIRLIST ${_DIR})
  endforeach()

  if(_DIRLIST)
    list(REMOVE_DUPLICATES _DIRLIST)
  endif()

  set(${result} ${_DIRLIST} PARENT_SCOPE)
endfunction()
