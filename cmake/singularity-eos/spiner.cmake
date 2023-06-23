macro(_singularity_set_spiner_opt singularity_opt spiner_opt)
  if(${singularity_opt})
    set(${spiner_opt}
        ON
        CACHE BOOL "" FORCE)
  endif()
endmacro()

macro(singularity_import_spiner)
  if(NOT TARGET spiner::spiner)
    _singularity_set_spiner_opt(SINGULARITY_USE_SPINER_WITH_HDF5 SPINER_USE_HDF)
    _singularity_set_spiner_opt(SINGULARITY_USE_CUDA SPINER_USE_CUDA)
    _singularity_set_spiner_opt(SINGULARITY_USE_KOKKOS SPINER_USE_KOKKOS)
    add_subdirectory(utils/spiner)
  endif()
endmacro()

macro(singularity_find_spiner)
  find_package(spiner REQUIRED)
endmacro()

macro(_singularity_chk_spiner_opt singularity_opt spiner_opt)
  if(${singularity_opt} AND NOT ${spiner_opt})
    message(
      WARNING
        "\"${singularity_opt}\" is ON,
      and a \"spiner\" target is already defined, but
      \"${spiner_opt}\" is either unset or set to OFF.
      The configuration will continue, but you will likely
      experience issues building and running.")
  endif()
endmacro()

macro(singularity_enable_spiner target)
  _singularity_chk_spiner_opt(SINGULARITY_USE_SPINER_WITH_HDF5 SPINER_USE_HDF)
  _singularity_chk_spiner_opt(SINGULARITY_USE_CUDA SPINER_USE_CUDA)
  _singularity_chk_spiner_opt(SINGULARITY_USE_KOKKOS SPINER_USE_KOKKOS)

  target_compile_definitions(
    ${target}
    PUBLIC
      $<$<BOOL:${SINGULARITY_USE_SPINER}>:SINGULARITY_USE_SPINER>
      $<$<BOOL:${SINGULARITY_USE_SPINER_WITH_HDF5}>:SINGULARITY_USE_SPINER_WITH_HDF5>
  )
  target_link_libraries(${target} PUBLIC spiner::spiner)
endmacro()
