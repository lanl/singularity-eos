#------------------------------------------------------------------------------
# © 2021. Triad National Security, LLC. All rights reserved.  This
# program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by Triad
# National Security, LLC for the U.S.  Department of Energy/National
# Nuclear Security Administration. All rights in the program are
# reserved by Triad National Security, LLC, and the U.S. Department of
# Energy/National Nuclear Security Administration. The Government is
# granted for itself and others acting on its behalf a nonexclusive,
# paid-up, irrevocable worldwide license in this material to reproduce,
# prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.
#------------------------------------------------------------------------------

if (NOT TARGET Catch2::Catch2)
  find_package(Catch2 QUIET)
endif()

if(NOT Catch2_FOUND)
  message(STATUS "Fetching Catch2 as needed")
  # idiomatic FetchContent
  include(FetchContent)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.13.7)
  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/contrib)
endif()

add_executable(eos_unit_tests
                 eos_unit_tests.cpp)

if(SINGULARITY_USE_HDF5)
  if(SINGULARITY_TEST_SESAME)
    target_compile_definitions(eos_unit_tests PRIVATE SINGULARITY_TEST_SESAME)
  endif()
  if(SINGULARITY_TEST_STELLAR_COLLAPSE)
    target_compile_definitions(eos_unit_tests PRIVATE SINGULARITY_TEST_STELLAR_COLLAPSE)
  endif()
endif()

target_link_libraries(eos_unit_tests 
  PRIVATE
    Catch2::Catch2
    singularity-eos::singularity-eos
)

if(SINGULARITY_USE_EOSPAC AND SINGULARITY_USE_HDF5 AND SINGULARITY_TEST_SESAME)
  add_executable(compare_to_eospac
                 compare_to_eospac.cpp)
  target_link_libraries(compare_to_eospac
    PRIVATE
      Catch2::Catch2
      singularity-eos::singularity-eos
  )
endif()

include(Catch)
catch_discover_tests(eos_unit_tests)

if (SINGULARITY_BUILD_CLOSURE)
  add_executable(test_pte
                 test_pte.cpp)
  target_link_libraries(test_pte PRIVATE
                        Catch2::Catch2
                        singularity-eos::singularity-eos)
  add_test(pte test_pte)
endif()

if (SINGULARITY_USE_HDF5 AND SINGULARITY_TEST_STELLAR_COLLAPSE)
  add_executable(profile_stellar_collapse
                 profile_stellar_collapse.cpp)
  target_link_libraries(profile_stellar_collapse singularity-eos::singularity-eos)

  message(STATUS "Generating stellar collapse table for tests using Python script")
  include(FindPythonInterp)
  if (HDF5_IS_PARALLEL)
    execute_process(COMMAND mpirun -n 1 ${PYTHON_EXECUTABLE}
      ${CMAKE_SOURCE_DIR}/utils/scripts/make_tabulated_ideal_sc.py
      -o ${PROJECT_BINARY_DIR}/stellar_collapse_ideal.h5)
  else()
    execute_process(COMMAND ${PYTHON_EXECUTABLE}
      ${CMAKE_SOURCE_DIR}/utils/scripts/make_tabulated_ideal_sc.py
      -o ${PROJECT_BINARY_DIR}/stellar_collapse_ideal.h5)
  endif()
endif()
