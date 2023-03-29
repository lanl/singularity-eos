macro(singularity_enable_mpark_variant target)
  
  find_package(mpark_variant REQUIRED)

  target_link_libraries(${target} PUBLIC mpark_variant)
endmacro()
