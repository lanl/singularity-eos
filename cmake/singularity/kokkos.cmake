macro(singularity_enable_kokkos target)

  if(SINGULARITY_CUSTOM_KOKKOS)
    message(WARNING "You have selected to use a custom source of Kokkos."
    " This is not recommended, and developer support for this feature is not secured.")

    _singularity_handle_custom_kokkos_tpl(
      NAME "Kokkos"
    )

    if(SINGULARITY_USE_KOKKOSKERNELS)
      message(STATUS "\"KokkosKernels\" requested with custom \"Kokkos\", directories and (optional) cache"
        " files will be queried.")
      _singularity_handle_custom_kokkos_tpl(
        NAME "KokkosKernels"
      )
    endif()
    
  else()
    if(NOT TARGET Kokkos::kokkos)
      find_package(Kokkos COMPONENTS separable_compilation REQUIRED)
      if(SINGULARITY_USE_KOKKOSKERNELS)
        if(NOT TARGET Kokkos::kokkoskernels)
          find_package(KokkosKernels REQUIRED)
        endif()
      endif()
    endif()
  endif()

  target_link_libraries(${target} PUBLIC Kokkos::kokkos)

  target_compile_definitions(${target} PUBLIC PORTABILITY_STRATEGY_KOKKOS)

  target_compile_definitions(${target} PUBLIC SINGULARITY_USE_KOKKOSKERNELS)
  #TODO: shouldn't be needed
  target_compile_definitions(${target} PUBLIC SPINER_USE_KOKKOS)

  if(SINGULARITY_USE_KOKKOSKERNELS AND SINGULARITY_BUILD_CLOSURE)
    target_link_libraries(${target} PUBLIC Kokkos::kokkoskernels)
  endif()
endmacro()

macro(_singularity_handle_custom_kokkos_tpl)
  set(options)
  set(one_value_args
    NAME
  )
  set(multi_value_args)

  cmake_parse_arguments(KTPL "${options}" "${one_value_args}"
    "${multi_value_args}" ${ARGN})

  if(NOT ("${KTPL_NAME}" STREQUAL "Kokkos" OR "${KTPL_NAME}" STREQUAL "KokkosKernels") )
    message(FATAL_ERROR "Custom TPL only supported for \"Kokkos\" and \"KokkosKernels\" (got \"${KTPL_NAME}\")")
  endif()

  string(TOUPPER ${KTPL_NAME} tplUp)
  string(TOLOWER ${KTPL_NAME} tplLo)
    
  # requires source dir
  if(NOT SINGULARITY_CUSTOM_${tplUp}_DIR)
    message(FATAL_ERROR "Selected custom ${KTPL_NAME} source, but SINGULARITY_CUSTOM_${tplUp}_DIR not set."
     "Set and reconfigure, or disable custom source.")
  else()
    if(NOT EXISTS ${SINGULARITY_CUSTOM_${tplUp}_DIR})
      message(FATAL_ERROR "Custom ${KTPL_NAME} source dir ${SINGULARITY_CUSTOM_${tplUP}_DIR} does not exist.")
    endif()
  endif()
  if(NOT SINGULARITY_CUSTOM_${tplUp}_CACHE)
    message(WARNING "Selected custom ${KTPL_NAME} source, but SINGULARITY_CUSTOM_${tplUp}_CACHE not set."
      " CMake will configure with default options."
      " See ${PROJECT_SOURCE_DIR}/utils/scripts/caches for sample caches.")
    set(SINGULARITY_CUSTOM_${tplUp}_CACHE "${PROJECT_SOURCE_DIR}/utils/scripts/caches/${tplLo}_default.cmake")
  endif() 
  if(NOT EXISTS ${SINGULARITY_CUSTOM_${tplUp}_CACHE})
    message(FATAL "Custom kokkos CMake cache file ${SINGULARITY_CUSTOM_${tplUp}_CACHE} does not exist.")
  endif()
  include(${SINGULARITY_CUSTOM_${tplUp}_CACHE})

  if(NOT TARGET Kokkos::${tplLo})
    add_subdirectory(${SINGULARITY_CUSTOM_${tplUp}_DIR} ${PROJECT_BINARY_DIR}/_deps/${tplLo}-build)
  else()
    message(WARNING "${KTPL_NAME} target already present, not using custom source directory.")
  endif()
 
endmacro()
