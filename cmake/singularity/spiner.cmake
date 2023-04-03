macro(singularity_import_spiner)
  if(NOT TARGET spiner::spiner)
    add_subdirectory(utils/spiner)
  endif()
endmacro()

macro(singularity_find_spiner)
    find_package(spiner REQUIRED)
endmacro()

macro(singularity_enable_spiner target)
  target_link_libraries(${target} PUBLIC spiner::spiner)
endmacro()
