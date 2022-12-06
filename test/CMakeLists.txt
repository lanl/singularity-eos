#------------------------------------------------------------------------------
# © 2021-2022. Triad National Security, LLC. All rights reserved.  This
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
  catch2_define.cpp
  eos_unit_test_helpers.hpp
  test_eos_unit.cpp
  test_eos_gruneisen.cpp
  test_eos_snlgrun.cpp
  test_eos_vinet.cpp
)

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
include(Catch)
catch_discover_tests(eos_unit_tests PROPERTIES TIMEOUT 60)

if(SINGULARITY_USE_EOSPAC AND SINGULARITY_USE_HDF5 AND SINGULARITY_TEST_SESAME)
  add_executable(compare_to_eospac
                 compare_to_eospac.cpp)
  target_link_libraries(compare_to_eospac
    PRIVATE
      Catch2::Catch2
      singularity-eos::singularity-eos
  )
endif()

if (SINGULARITY_BUILD_CLOSURE)
  add_executable(test_pte
                 test_pte.cpp)
  target_link_libraries(test_pte PRIVATE
                        Catch2::Catch2
                        singularity-eos::singularity-eos)
  add_test(pte test_pte)
endif()

if (SINGULARITY_USE_HDF5 AND SINGULARITY_TEST_STELLAR_COLLAPSE)
  set(SINGULARITY_GOLDFILE_URL
    https://github.com/lanl/singularity-eos/releases/download/${SINGULARITY_GOLDFILES_VERSION}/goldfiles.tar.gz)
  message(STATUS "Attempting to download gold files for regression tests. "
    "File is located at ${SINGULARITY_GOLDFILE_URL} "
    "If you do not want to use the gold file, you can set -DSINGULARITY_TEST_STELLAR_COLLAPSE=OFF."
    )
  include(FetchContent)
  FetchContent_Declare(goldfiles
    URL ${SINGULARITY_GOLDFILE_URL}
    URL_HASH SHA256=${SINGULARITY_GOLDFILE_HASH}
    SOURCE_DIR goldfiles
    )
  FetchContent_MakeAvailable(goldfiles)

  add_executable(profile_stellar_collapse
                 profile_stellar_collapse.cpp)
  target_link_libraries(profile_stellar_collapse singularity-eos::singularity-eos)
endif()

if (SINGULARITY_BUILD_PYTHON AND SINGULARITY_TEST_PYTHON)
  find_package(Python COMPONENTS Interpreter Development REQUIRED)
  set(PYTHON_TEST_ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH})

  add_test(NAME PythonBindings
           COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python_bindings.py -v)
  set_tests_properties(PythonBindings PROPERTIES ENVIRONMENT "${PYTHON_TEST_ENVIRONMENT}")
endif()
