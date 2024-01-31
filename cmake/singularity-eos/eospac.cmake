macro(singularity_enable_eospac target)
  if(NOT TARGET EOSPAC::eospac)
    find_package(EOSPAC REQUIRED)
  endif()
  target_link_libraries(${target} INTERFACE EOSPAC::eospac)
  target_compile_definitions(${target} INTERFACE SINGULARITY_USE_EOSPAC)

  add_subdirectory(${PROJECT_SOURCE_DIR}/eospac-wrapper)
  target_link_libraries(${target} INTERFACE eospac-wrapper)
endmacro()


