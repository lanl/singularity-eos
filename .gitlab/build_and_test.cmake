set(CTEST_SOURCE_DIRECTORY "$ENV{SOURCE_DIR}")
set(CTEST_BINARY_DIRECTORY "$ENV{BUILD_DIR}")
set(CTEST_PROJECT_NAME "$ENV{PROJECT_NAME}")


set(CTEST_SITE "$ENV{SYSTEM_NAME}")
if("$ENV{CI_PIPELINE_SOURCE}" STREQUAL "merge_request_event")
  set(CTEST_BUILD_NAME "$ENV{SPACK_ENV_NAME}, source=$ENV{CI_MERGE_REQUEST_SOURCE_BRANCH_NAME}, target=$ENV{CI_MERGE_REQUEST_TARGET_BRANCH_NAME}, MR=$ENV{CI_MERGE_REQUEST_IID}")
else()
  set(CTEST_BUILD_NAME "$ENV{SPACK_ENV_NAME}, Branch=$ENV{CI_COMMIT_REF_NAME}")
endif()

set(CTEST_SUBMIT_URL "$ENV{CDASH_SERVER_URL}/submit.php?project=${CTEST_PROJECT_NAME}")
set(CTEST_CURL_OPTIONS CURLOPT_SSL_VERIFYPEER_OFF CURLOPT_SSL_VERIFYHOST_OFF)
set(CTEST_NIGHTLY_START_TIME "22:00:00 MDT")
set(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE 1024000)
set(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE 1024000)

set(CTEST_UPDATE_COMMAND "git")
set(CTEST_GIT_UPDATE_CUSTOM "${CMAKE_COMMAND}" "-E" "echo" "Skipping git update (no-op).")
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_BUILD_CONFIGURATION "RelWithDebInfo")

set(CTEST_MODE "$ENV{CTEST_MODE}")

set(CTEST_OUTPUT_ON_FAILURE ON)
set(CTEST_USE_LAUNCHERS TRUE)

cmake_host_system_information(RESULT NUM_PHYSICAL_CORES QUERY NUMBER_OF_PHYSICAL_CORES)

if(${NUM_PHYSICAL_CORES} EQUAL 1)
  # workaround for GraceHopper
  include(ProcessorCount)
  ProcessorCount(NUM_PHYSICAL_CORES)
endif()

if(${CTEST_SCRIPT_ARG} MATCHES Configure)
  ctest_start(${CTEST_MODE})
else()
  ctest_start(${CTEST_MODE} APPEND)
endif()

if(${CTEST_SCRIPT_ARG} MATCHES Configure)
  ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}")
  ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE configure_error)

  if(configure_error)
    if(${CTEST_SCRIPT_ARG} MATCHES ReportErrors)
        ctest_submit()
    endif()
    message(FATAL_ERROR "configure error")
  endif()
endif()

if(${CTEST_SCRIPT_ARG} MATCHES Build)
  ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" PARALLEL_LEVEL ${NUM_PHYSICAL_CORES} RETURN_VALUE build_error)

  if(build_error)
    if(${CTEST_SCRIPT_ARG} MATCHES ReportErrors)
        ctest_submit()
    endif()
    message(FATAL_ERROR "build error")
  endif()
endif()

if(${CTEST_SCRIPT_ARG} MATCHES Test)
  ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" OUTPUT_JUNIT tests.xml RETURN_VALUE test_error)

  if(test_error)
    if(${CTEST_SCRIPT_ARG} MATCHES ReportErrors)
        ctest_submit()
    endif()
    message(FATAL_ERROR "test error")
  endif()
endif()

if(${CTEST_SCRIPT_ARG} MATCHES Submit)
  ctest_submit()
endif()
