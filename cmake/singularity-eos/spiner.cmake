macro(singularity_import_spiner)
  if(NOT TARGET spiner::spiner)
    if(SINGULARITY_SPINER_USE_HDF5)
      set(SPINER_USE_HDF ON CACHE BOOL "" FORCE)
    endif()
    add_subdirectory(utils/spiner)
  endif()
endmacro()

macro(singularity_find_spiner)
    find_package(spiner REQUIRED)
endmacro()

macro(singularity_enable_spiner target)
  if(SINGULARITY_USE_SPINER_HDF5 AND NOT SPINER_USE_HDF)
    message(WARNING "\"SINGULARITY_USE_SPINER_HDF5\" is ON, 
      and a \"spiner\" target is already defined, but 
      \"SPINER_USE_HDF\" is either unset or set to OFF.
      The configuration will continue, but you will likely
      experience issues building and running.")
  endif()
  #TODO
  # This is redundent; `spiner::spiner` should already have this,
  # but need to double-check that both the submodule and the 
  # exported system config of `spiner` includes it.
  # There isn't any harm in potentially doing this twice, but 
  # a future cleanup should resolve this
  target_compile_definitions(${target} PUBLIC SPINER_USE_HDF)
  target_link_libraries(${target} PUBLIC spiner::spiner)
endmacro()
