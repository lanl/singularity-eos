macro(singularity_import_eigen)
  if(NOT TARGET Eigen3::Eigen)
    set(_storBT ${BUILD_TESTING})
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

    set(EIGEN_TEST_CXX11 OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)

    add_subdirectory(utils/eigen)
    
    set(BUILD_TESTING ${_storBT} CACHE BOOL "" FORCE)
    unset(_storBT)
  endif()
endmacro()

macro(singularity_find_eigen)
    find_package(Eigen3 REQUIRED)
endmacro()

macro(singularity_enable_eigen target)
  target_link_libraries(${target} PUBLIC Eigen3::Eigen)
endmacro()
