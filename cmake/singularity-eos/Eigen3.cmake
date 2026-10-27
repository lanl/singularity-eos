macro(singularity_find_eigen)
  if(NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 QUIET)
  endif()

  if(NOT TARGET Eigen3::Eigen)
    if(NOT SINGULARITY_FETCH_EIGEN)
      message(FATAL_ERROR
        "Could not find Eigen3. Install an Eigen3 package or configure with "
        "-DSINGULARITY_FETCH_EIGEN=ON to download it.")
    endif()

    include(FetchContent)
    set(_storBT ${BUILD_TESTING})
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

    set(EIGEN_TEST_CXX11 OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)

    FetchContent_Declare(
      eigen
      GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
      GIT_TAG 2859db0220cd8644c1a75a3bf04f62f551f73f22)
    FetchContent_MakeAvailable(eigen)

    set(BUILD_TESTING ${_storBT} CACHE BOOL "" FORCE)
    unset(_storBT)
  endif()

  if(NOT TARGET Eigen3::Eigen)
    message(FATAL_ERROR "Eigen3 did not provide the Eigen3::Eigen target.")
  endif()
endmacro()

macro(singularity_enable_eigen target)
  target_link_libraries(${target} INTERFACE Eigen3::Eigen)
endmacro()
