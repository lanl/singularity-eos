macro(singularity_enable_eospac target)
  if(NOT TARGET EOSPAC::eospac)
    find_package(EOSPAC REQUIRED)
  endif()
  target_link_libraries(${target} PUBLIC EOSPAC::eospac)
  target_compile_definitions(${target} PUBLIC SINGULARITY_USE_EOSPAC)

  add_subdirectory(${PROJECT_SOURCE_DIR}/eospac-wrapper)
  target_link_libraries(${target} PUBLIC eospac-wrapper)
endmacro()


