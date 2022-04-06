#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "singularity-eos::singularity-eos" for configuration "Debug"
set_property(TARGET singularity-eos::singularity-eos APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(singularity-eos::singularity-eos PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libsingularity-eos.so.1.5.0"
  IMPORTED_SONAME_DEBUG "libsingularity-eos.so.1"
  )

list(APPEND _IMPORT_CHECK_TARGETS singularity-eos::singularity-eos )
list(APPEND _IMPORT_CHECK_FILES_FOR_singularity-eos::singularity-eos "${_IMPORT_PREFIX}/lib/libsingularity-eos.so.1.5.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
