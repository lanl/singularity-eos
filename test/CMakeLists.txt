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

add_executable(eos_unit_tests
                 eos_unit_tests.cpp)
if(SINGULARITY_USE_HDF5 AND SINGULARITY_TEST_SESAME)
  target_compile_definitions(eos_unit_tests PRIVATE SINGULARITY_TEST_SESAME)
endif()
if(SINGULARITY_USE_HDF5 AND SINGULARITY_TEST_STELLAR_COLLAPSE)
  target_compile_definitions(eos_unit_tests PRIVATE SINGULARITY_TEST_STELLAR_COLLAPSE)
endif()
target_link_libraries(eos_unit_tests Catch2::Catch2 eos singularity-eos::libs singularity-eos::flags)

if(SINGULARITY_USE_EOSPAC AND SINGULARITY_USE_HDF5 AND SINGULARITY_TEST_SESAME)
  add_executable(compare_to_eospac
                 compare_to_eospac.cpp)
  target_link_libraries(compare_to_eospac Catch2::Catch2 eos singularity-eos::libs singularity-eos::flags)
  target_include_directories(compare_to_eospac PUBLIC
                             ${CMAKE_SOURCE_DIR}/utils/spiner
                             ${CMAKE_SOURCE_DIR}/utils/ports-of-call
                             ${CMAKE_SOURCE_DIR}/utils/sesame2spiner)
endif()
catch_discover_tests(eos_unit_tests)

if (SINGULARITY_BUILD_CLOSURE)
  add_executable(test_pte
                 test_pte.cpp)
  target_link_libraries(test_pte PRIVATE
                        Catch2::Catch2
                        eos
                        singularity-eos::libs
                        singularity-eos::flags)
  add_test(pte test_pte)
endif()

if (SINGULARITY_USE_HDF5 AND SINGULARITY_TEST_STELLAR_COLLAPSE)
  add_executable(profile_stellar_collapse
                 profile_stellar_collapse.cpp)
  target_link_libraries(profile_stellar_collapse eos singularity-eos::libs singularity-eos::flags)

  message(STATUS "Generating stellar collapse table for tests using Python script")
  include(FindPythonInterp)
  if (HDF5_IS_PARALLEL)
    execute_process(COMMAND mpirun -n 1 ${PYTHON_EXECUTABLE}
      ${CMAKE_SOURCE_DIR}/utils/scripts/make_tabulated_gamma.py
      -o ${PROJECT_BINARY_DIR}/stellar_collapse_ideal.h5)
  else()
    execute_process(COMMAND ${PYTHON_EXECUTABLE}
      ${CMAKE_SOURCE_DIR}/utils/scripts/make_tabulated_gamma.py
      -o ${PROJECT_BINARY_DIR}/stellar_collapse_ideal.h5)
  endif()
endif()
