macro(singularity_enable_eospac target)
  find_package(EOSPAC REQUIRED)
  target_link_libraries(${target} PUBLIC EOSPAC::eospac)
  target_compile_defines(${target} PUBLIC SINGULARITY_USE_EOSPAC)

  add_subdirectory(${PROJECT_SOURCE_DIR}/eospac-wrapper)
  target_link_libraries(${target} PUBLIC eospac-wrapper)
endmacro()


