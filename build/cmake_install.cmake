# Install script for directory: /home/mauneyc/devel/singularity_local/singularity-eos

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ports-of-call" TYPE DIRECTORY FILES "/home/mauneyc/devel/singularity_local/singularity-eos/utils/ports-of-call")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/mauneyc/devel/singularity_local/singularity-eos/build/utils/variant/cmake_install.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/fast-math" TYPE DIRECTORY FILES "/home/mauneyc/devel/singularity_local/singularity-eos/utils/fast-math")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/root-finding-1d" TYPE DIRECTORY FILES "/home/mauneyc/devel/singularity_local/singularity-eos/utils/root-finding-1d")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/spiner" TYPE DIRECTORY FILES "/home/mauneyc/devel/singularity_local/singularity-eos/utils/spiner")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/mauneyc/devel/singularity_local/singularity-eos/build/singularity-eos/base/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/mauneyc/devel/singularity_local/singularity-eos/build/singularity-eos/closure/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/mauneyc/devel/singularity_local/singularity-eos/build/singularity-eos/eos/cmake_install.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsingularity-eos.so.1.5.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsingularity-eos.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/mauneyc/devel/singularity_local/singularity-eos/build/libsingularity-eos.so.1.5.0"
    "/home/mauneyc/devel/singularity_local/singularity-eos/build/libsingularity-eos.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsingularity-eos.so.1.5.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsingularity-eos.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsingularity-eos.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsingularity-eos.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsingularity-eos.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/mauneyc/devel/singularity_local/singularity-eos/build/libsingularity-eos.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsingularity-eos.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsingularity-eos.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsingularity-eos.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/singularity-eos/singularity-eosTargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/singularity-eos/singularity-eosTargets.cmake"
         "/home/mauneyc/devel/singularity_local/singularity-eos/build/CMakeFiles/Export/lib/cmake/singularity-eos/singularity-eosTargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/singularity-eos/singularity-eosTargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/singularity-eos/singularity-eosTargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/singularity-eos" TYPE FILE FILES "/home/mauneyc/devel/singularity_local/singularity-eos/build/CMakeFiles/Export/lib/cmake/singularity-eos/singularity-eosTargets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/singularity-eos" TYPE FILE FILES "/home/mauneyc/devel/singularity_local/singularity-eos/build/CMakeFiles/Export/lib/cmake/singularity-eos/singularity-eosTargets-debug.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/singularity-eos/closure" TYPE FILE FILES "/home/mauneyc/devel/singularity_local/singularity-eos/singularity-eos/closure/mixed_cell_models.hpp")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/singularity-eos/eos" TYPE FILE FILES "/home/mauneyc/devel/singularity_local/singularity-eos/singularity-eos/eos/eos_builder.hpp")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/singularity-eos/eos" TYPE FILE FILES "/home/mauneyc/devel/singularity_local/singularity-eos/singularity-eos/eos/eos_variant.hpp")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/singularity-eos/eos" TYPE FILE FILES "/home/mauneyc/devel/singularity_local/singularity-eos/singularity-eos/eos/eos.hpp")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/singularity-eos/eos/modifiers" TYPE FILE FILES "/home/mauneyc/devel/singularity_local/singularity-eos/singularity-eos/eos/modifiers/eos_unitsystem.hpp")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/singularity-eos/eos/modifiers" TYPE FILE FILES "/home/mauneyc/devel/singularity_local/singularity-eos/singularity-eos/eos/modifiers/scaled_eos.hpp")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/singularity-eos/eos/modifiers" TYPE FILE FILES "/home/mauneyc/devel/singularity_local/singularity-eos/singularity-eos/eos/modifiers/shifted_eos.hpp")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/singularity-eos/eos/modifiers" TYPE FILE FILES "/home/mauneyc/devel/singularity_local/singularity-eos/singularity-eos/eos/modifiers/relativistic_eos.hpp")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/singularity-eos" TYPE FILE FILES "/home/mauneyc/devel/singularity_local/singularity-eos/build/singularity-eosConfig.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/singularity-eos" TYPE FILE FILES "/home/mauneyc/devel/singularity_local/singularity-eos/build/singularity-eosConfigVersion.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/singularity-eos" TYPE FILE FILES
    "/home/mauneyc/devel/singularity_local/singularity-eos/cmake/compilers.cmake"
    "/home/mauneyc/devel/singularity_local/singularity-eos/cmake/cuda.cmake"
    "/home/mauneyc/devel/singularity_local/singularity-eos/cmake/eigen3.cmake"
    "/home/mauneyc/devel/singularity_local/singularity-eos/cmake/eospac.cmake"
    "/home/mauneyc/devel/singularity_local/singularity-eos/cmake/FindEigen3.cmake"
    "/home/mauneyc/devel/singularity_local/singularity-eos/cmake/FindEOSPAC.cmake"
    "/home/mauneyc/devel/singularity_local/singularity-eos/cmake/FindHwloc.cmake"
    "/home/mauneyc/devel/singularity_local/singularity-eos/cmake/hdf5.cmake"
    "/home/mauneyc/devel/singularity_local/singularity-eos/cmake/kokkos.cmake"
    "/home/mauneyc/devel/singularity_local/singularity-eos/cmake/kokkoskernels.cmake"
    "/home/mauneyc/devel/singularity_local/singularity-eos/cmake/library.cmake"
    "/home/mauneyc/devel/singularity_local/singularity-eos/cmake/mathselect.cmake"
    "/home/mauneyc/devel/singularity_local/singularity-eos/cmake/portsofcall.cmake"
    "/home/mauneyc/devel/singularity_local/singularity-eos/cmake/subdirlist.cmake"
    "/home/mauneyc/devel/singularity_local/singularity-eos/cmake/variant.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/mauneyc/devel/singularity_local/singularity-eos/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
