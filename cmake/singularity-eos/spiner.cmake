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
    # TODO I dunno if these are used
    _singularity_set_spiner_opt(SINGULARITY_USE_KOKKOS SPINER_USE_KOKKOS)
    _singularity_set_spiner_opt(SINGULARITY_USE_CUDA SPINER_USE_CUDA)
    add_subdirectory(utils/spiner)
  endif()
endmacro()

macro(singularity_find_spiner)
  find_package(spiner 1.6 REQUIRED)
endmacro()

macro(_singularity_chk_spiner_opt singularity_opt spiner_opt)
  if(${singularity_opt} AND NOT ${spiner_opt})
    message(
      FATAL_ERROR
        "\"${singularity_opt}\" is ON,
        but \"${spiner_opt}\" is OFF.
      Please ensure your \"spiner\" installation is configured with
      \"-D${spiner_opt}=ON\"")
  endif()
endmacro()

macro(singularity_enable_spiner target)
  _singularity_chk_spiner_opt(SINGULARITY_USE_SPINER_WITH_HDF5 SPINER_USE_HDF)
  # _singularity_chk_spiner_opt(SINGULARITY_USE_CUDA SPINER_USE_CUDA)
  # _singularity_chk_spiner_opt(SINGULARITY_USE_KOKKOS SPINER_USE_KOKKOS)

  target_compile_definitions(
    ${target}
    PUBLIC
      $<$<BOOL:${SINGULARITY_USE_SPINER}>:SINGULARITY_USE_SPINER>
      $<$<BOOL:${SINGULARITY_USE_SPINER_WITH_HDF5}>:SINGULARITY_USE_SPINER_WITH_HDF5>
  )
  target_link_libraries(${target} PUBLIC spiner::spiner)
endmacro()
