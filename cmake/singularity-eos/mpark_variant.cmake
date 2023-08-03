macro(singularity_import_mpark_variant)
  # Patches variant to be compatible with cuda Assumes "patch" is present on
  # system
  if(SINGULARITY_PATCH_MPARK_VARIANT)
    message(STATUS "Patching variant")
    execute_process(
      COMMAND
        patch -N -s -V never
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/variant/include/mpark/variant.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/cuda_compatibility.patch)
  endif()
  if(NOT TARGET mpark_variant)
    add_subdirectory(utils/variant)
  endif()
endmacro()

macro(singularity_find_mpark_variant)
  find_package(mpark_variant REQUIRED)
  # TODO Check that imported `mpark_variant` has been correctly patched This may
  # require some additions upstream
endmacro()

macro(singularity_enable_mpark_variant target)
  target_link_libraries(${target} PUBLIC mpark_variant)
endmacro()
