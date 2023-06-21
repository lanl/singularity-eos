macro(singularity_import_ports_of_call)
  if(NOT TARGET ports-of-call::ports-of-call)
    add_subdirectory(utils/ports-of-call)
  endif()
endmacro()

macro(singularity_find_ports_of_call)
  find_package(ports-of-call REQUIRED)
endmacro()

macro(singularity_enable_ports_of_call target)
  target_link_libraries(${target} PUBLIC ports-of-call::ports-of-call)
endmacro()
