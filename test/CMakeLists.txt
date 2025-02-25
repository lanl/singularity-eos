# ------------------------------------------------------------------------------
# © 2021-2024. Triad National Security, LLC. All rights reserved.  This program
# was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
# National Laboratory (LANL), which is operated by Triad National Security, LLC
# for the U.S.  Department of Energy/National Nuclear Security Administration.
# All rights in the program are reserved by Triad National Security, LLC, and
# the U.S. Department of Energy/National Nuclear Security Administration. The
# Government is granted for itself and others acting on its behalf a
# nonexclusive, paid-up, irrevocable worldwide license in this material to
# reproduce, prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.
# ------------------------------------------------------------------------------

add_executable(
  eos_analytic_unit_tests
  catch2_define.cpp
  eos_unit_test_helpers.hpp
  test_eos_ideal.cpp
  test_eos_gruneisen.cpp
  test_eos_sap_polynomial.cpp
  test_eos_noble_abel.cpp
  test_eos_stiff.cpp
  test_eos_carnahan_starling.cpp
  test_eos_vinet.cpp
  test_eos_mgusup.cpp
  test_eos_powermg.cpp
  $<$<BOOL:${SINGULARITY_BUILD_CLOSURE}>:test_pte_ideal.cpp>
  test_eos_zsplit.cpp
  )

add_executable(
  eos_infrastructure_tests
  catch2_define.cpp
  eos_unit_test_helpers.hpp
  #test_eos_modifiers.cpp
  test_eos_modifiers_minimal.cpp
  test_eos_vector.cpp
  test_math_utils.cpp
  test_variadic_utils.cpp
  test_bounds.cpp
  test_indexable_types.cpp
  )

add_executable(
  eos_tabulated_unit_tests
  catch2_define.cpp
  eos_unit_test_helpers.hpp
  test_spiner_tricks.cpp
  test_eos_helmholtz.cpp
  test_eos_tabulated.cpp
  test_eos_stellar_collapse.cpp
  )

if (SINGULARITY_BUILD_CLOSURE)
  add_executable(
    closure_unit_tests
    catch2_define.cpp
    test_closure_pte.cpp
    test_kpt_models.cpp
    )
endif()

get_property(plugin_tests GLOBAL PROPERTY PLUGIN_TESTS)
if (plugin_tests)
  add_executable(
    eos_plugin_tests
    catch2_define.cpp
    eos_unit_test_helpers.hpp
    ${plugin_tests}
    )
endif()

if(SINGULARITY_TEST_HELMHOLTZ)
    configure_file(${PROJECT_SOURCE_DIR}/data/helmholtz/helm_table.dat ${CMAKE_BINARY_DIR}/data/helmholtz/helm_table.dat COPYONLY)
  target_compile_definitions(eos_tabulated_unit_tests PRIVATE SINGULARITY_TEST_HELMHOLTZ SINGULARITY_USE_HELMHOLTZ)
endif()

if(SINGULARITY_TEST_SESAME)
  target_compile_definitions(eos_tabulated_unit_tests PRIVATE SINGULARITY_TEST_SESAME)
endif()
if(SINGULARITY_TEST_STELLAR_COLLAPSE)
  target_compile_definitions(eos_tabulated_unit_tests
                             PRIVATE SINGULARITY_TEST_STELLAR_COLLAPSE)
endif()

target_link_libraries(eos_analytic_unit_tests PRIVATE Catch2::Catch2
  singularity-eos::singularity-eos)
target_link_libraries(eos_infrastructure_tests PRIVATE Catch2::Catch2
  singularity-eos::singularity-eos)
target_link_libraries(eos_tabulated_unit_tests PRIVATE Catch2::Catch2
  singularity-eos::singularity-eos)
if (plugin_tests)
  target_link_libraries(eos_plugin_tests PRIVATE Catch2::Catch2
    singularity-eos::singularity-eos)
endif()
include(Catch)
catch_discover_tests(eos_analytic_unit_tests PROPERTIES TIMEOUT 60)
catch_discover_tests(eos_infrastructure_tests PROPERTIES TIMEOUT 60)
catch_discover_tests(eos_tabulated_unit_tests PROPERTIES TIMEOUT 120)
if (plugin_tests)
  catch_discover_tests(eos_plugin_tests PROPERTIES TIMEOUT 60)
endif()

if(SINGULARITY_USE_EOSPAC AND SINGULARITY_TEST_SESAME)
  add_executable(compare_to_eospac compare_to_eospac.cpp)
  target_link_libraries(compare_to_eospac
                        PRIVATE Catch2::Catch2 singularity-eos::singularity-eos)
endif()

if(SINGULARITY_BUILD_CLOSURE)
  if (SINGULARITY_TEST_SESAME)
    target_compile_definitions(closure_unit_tests PRIVATE SINGULARITY_TEST_SESAME)
  endif()
  target_link_libraries(closure_unit_tests PRIVATE Catch2::Catch2
    singularity-eos::singularity-eos)
  catch_discover_tests(closure_unit_tests PROPERTIES TIMEOUT 120)

  if(SINGULARITY_USE_SPINER)
    add_executable(test_pte test_pte.cpp)
    target_link_libraries(test_pte PRIVATE Catch2::Catch2
                                         singularity-eos::singularity-eos)
    target_include_directories(test_pte PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(pte test_pte)

    add_executable(test_pte_2phase test_pte_2phase.cpp)
    target_link_libraries(test_pte_2phase PRIVATE Catch2::Catch2
                                         singularity-eos::singularity-eos)
    target_include_directories(test_pte_2phase PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(pte_2phase test_pte_2phase)

    if(SINGULARITY_USE_EOSPAC)
      add_executable(test_pte_3phase test_pte_3phase.cpp)
      target_link_libraries(test_pte_3phase PRIVATE Catch2::Catch2
                                           singularity-eos::singularity-eos)
      target_include_directories(test_pte_3phase PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
      add_test(pte_3phase test_pte_3phase)
    endif()
  endif()
  # disable fortran tests when using OpenMP due to segfaults
  if(SINGULARITY_USE_KOKKOS AND SINGULARITY_BUILD_FORTRAN_BACKEND AND NOT Kokkos_ENABLE_OPENMP)
    add_executable(test_get_sg_eos test_get_sg_eos.cpp)
    target_link_libraries(test_get_sg_eos PRIVATE Catch2::Catch2
                                         singularity-eos::singularity-eos)
    target_include_directories(test_get_sg_eos PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(get_sg_eos test_get_sg_eos)
  endif()
endif()

if(SINGULARITY_TEST_STELLAR_COLLAPSE)
  set(SINGULARITY_GOLDFILE_URL
      "https://github.com/lanl/singularity-eos/releases/download/${SINGULARITY_GOLDFILES_VERSION}/goldfiles.tar.gz"
      CACHE STRING "URL to goldfiles.tar.gz"
  )
  mark_as_advanced(SINGULARITY_GOLDFILE_URL)
  message(
    STATUS
      "Attempting to download gold files for regression tests. "
      "File is located at ${SINGULARITY_GOLDFILE_URL} "
      "If you do not want to use the gold file, you can set -DSINGULARITY_TEST_STELLAR_COLLAPSE=OFF."
  )
  include(FetchContent)
  FetchContent_Declare(
    goldfiles
    URL ${SINGULARITY_GOLDFILE_URL}
    URL_HASH SHA256=${SINGULARITY_GOLDFILE_HASH}
    SOURCE_DIR goldfiles)
  FetchContent_MakeAvailable(goldfiles)

  add_executable(profile_stellar_collapse profile_stellar_collapse.cpp)
  target_link_libraries(profile_stellar_collapse
                        singularity-eos::singularity-eos)
endif()

if(SINGULARITY_BUILD_PYTHON AND SINGULARITY_TEST_PYTHON)
  find_package(
    Python
    COMPONENTS Interpreter Development
    REQUIRED)
  set(PYTHON_TEST_ENVIRONMENT
      PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH})

  add_test(NAME PythonBindings
           COMMAND ${Python_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/python_bindings.py -v)
  set_tests_properties(PythonBindings PROPERTIES ENVIRONMENT
                                                 "${PYTHON_TEST_ENVIRONMENT}")
endif()

if(SINGULARITY_USE_SPINER)
  add_executable(profile_eos profile_eos.cpp)
  target_link_libraries(profile_eos singularity-eos::singularity-eos)
endif()

if(SINGULARITY_USE_FORTRAN)
  add_executable(ftn_interface test_f_iface.f90)
  target_link_libraries(ftn_interface singularity-eos::singularity-eos)
  set_target_properties(ftn_interface PROPERTIES LINKER_LANGUAGE Fortran)
  add_test(test_fortran_interface ftn_interface)
endif()

if(SINGULARITY_USE_EOSPAC
   AND SINGULARITY_TEST_SESAME
   AND NOT SINGULARITY_USE_CUDA)
  set(METHODS_TO_TEST
      TemperatureFromDensityInternalEnergy
      PressureFromDensityTemperature
      InternalEnergyFromDensityTemperature
      PressureFromDensityInternalEnergy
      SpecificHeatFromDensityTemperature
      SpecificHeatFromDensityInternalEnergy
      BulkModulusFromDensityTemperature
      BulkModulusFromDensityInternalEnergy
      GruneisenParamFromDensityTemperature
      GruneisenParamFromDensityInternalEnergy)

  foreach(MUT IN LISTS METHODS_TO_TEST)
    add_executable(profile_eospac_${MUT} profile_eospac.cpp)
    target_compile_definitions(profile_eospac_${MUT}
                               PUBLIC METHOD_UNDER_TEST=${MUT})
    target_link_libraries(profile_eospac_${MUT}
                          singularity-eos::singularity-eos)

    add_test(NAME Profile_EOSPAC_${MUT} COMMAND profile_eospac_${MUT} 100 16)
  endforeach()
endif()
