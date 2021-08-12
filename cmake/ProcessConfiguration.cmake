# Feature Summary used to present the user
# with an easy-to-read summary of options
# and features of the build configuration
#include(FeatureSummary)
include(CMakeDependentOption)
include(ExternalConfig)

option (SINGULARITY_USE_INTERNAL_DEPS "Use submodules when available" OFF)

option (SINGULARITY_USE_KOKKOS "Use Kokkos for portability" OFF)
cmake_dependent_option (SINGULARITY_USE_CUDA "Enable cuda support" ON "SINGULARITY_USE_KOKKOS" OFF)
cmake_dependent_option (SINGULARITY_USE_KOKKOSKERNELS "Use kokkos-kernels for linear algebra" ON "SINGULARITY_USE_KOKKOS" OFF)

option (SINGULARITY_USE_EOSPAC "Pull in eospac" OFF)
cmake_dependent_option (SINGULARITY_INVERT_AT_SETUP "Use eospacs preinverted tables" OFF "SINGULARITY_USE_EOSPAC" OFF)

option (SINGULARITY_USE_HDF5 "Pull in hdf5" ON)
option (SINGULARITY_USE_FORTRAN "Enable fortran bindings" ON)

option (SINGULARITY_BUILD_SESAME2SPINER "Compile sesame2spiner" OFF)
option (SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER "Compile stellarcollapse2spiner" OFF)

option (SINGULARITY_BETTER_DEBUG_FLAGS
  "Better debug flags for singularity" ON)
option (SINGULARITY_HIDE_MORE_WARNINGS "hide more warnings" OFF)
option (SINGULARITY_BUILD_CLOSURE "Mixed cell closure" ON)

option (SINGULARITY_BUILD_TESTS "Compile tests" OFF)
option (SINGULARITY_TEST_SESAME "Test the Sesame table readers" OFF)
option (SINGULARITY_TEST_STELLAR_COLLAPSE "Test the stellar collapse table readers" OFF)

option (SINGULARITY_EOS_SKIP_EXTRAP "Turn off eospac extrap checks" OFF)
option (SINGULARITY_USE_SINGLE_LOGS "Use single precision logs. Can harm accuracy." OFF)
option (SINGULARITY_FMATH_USE_ORDER_4 "4th order interpolant for fast logs. Default is 7th order." OFF)
option (SINGULARITY_FMATH_USE_ORDER_5 "5th order interpolant for fast logs. Default is 7th order." OFF)

if(SINGULARITY_USE_CUDA)
    if(NOT TARGET CUDA::toolkit)
        find_package(CUDAToolkit REQUIRED)
    endif()
endif()


if(SINGULARITY_USE_KOKKOS)
    if(NOT TARGET Kokkos::kokkos)
        find_package(Kokkos QUIET)
        # build from
        if(NOT Kokkos_FOUND)
            singularity_config_internal_kokkos(SINGULARITY_USE_CUDA)
            singularity_select_dep(
                USE_INTERNAL ${SINGULARITY_USE_INTERNAL_DEPS}
                DEP Kokkos
                GITURL https://github.com/kokkos/kokkos.git
                GITBRANCH master
                SUBDIR kokkos
            )
        endif()
        message(STATUS "[Kokkos] Kokkos setup complete")
    endif()
    if(SINGULARITY_USE_KOKKOSKERNELS)
        if(NOT TARGET Kokkos::kokkoskernels)
            find_package(KokkosKernels QUIET)
            if(NOT KokkosKernels_FOUND)
                singularity_config_internal_kokkoskernels()
                singularity_select_dep(
                    USE_INTERNAL OFF
                    DEP KokkosKernels
                    GITURL https://github.com/kokkos/kokkos-kernels.git
                    GITBRANCH master
                    SUBDIR ""
                )
            endif()
            message(STATUS "[KokkosKernels] KokkosKernels setup complete")
        endif()
    endif()
else()
    if(SINGULARITY_USE_CUDA)
        message(FATAL_ERROR "CUDA must be used with Kokkos")
    endif()
endif()

if(NOT SINGULARITY_USE_KOKKOSKERNELS)
    if(NOT TARGET Eigen3::Eigen)
        find_package(Eigen3 QUIET)
        if(NOT Eigen3_FOUND)
            singularity_select_dep(
                USE_INTERNAL ${SINGULARITY_USE_INTERNAL_DEPS}
                DEP Eigen3
                GITURL https://gitlab.com/libeigen/eigen.git
                GITBRANCH  master
                SUBDIR "eigen"
            )
        endif()
    endif()
    message(STATUS "[Eigen3] Setup complete")
endif()

if(SINGULARITY_USE_EOSPAC)
    find_package(EOSPAC REQUIRED)
endif()

if(SINGULARITY_BUILD_SESAME2SPINER)
    find_package(json QUIET)
    if(NOT json_FOUND)
        singularity_config_internal_json()
        singularity_select_dep(
            USE_INTERNAL ${SINGULARITY_USE_INTERNAL_DEPS}
            DEP json
            GITURL https://github.com/nlohmann/json.git
            GITBRANCH  develop
            SUBDIR "json"
        )
    endif()
    message(STATUS "[json] Setup complete")
endif()

find_package(PortsofCall REQUIRED)

if(SINGULARITY_BUILD_TESTS)
    if(NOT TARGET Catch2::Catch2)
        find_package(Catch2 QUIET)
        if(NOT Catch2_FOUND)
            singularity_select_dep(
                USE_INTERNAL ${SINGULARITY_USE_INTERNAL_DEPS}
                DEP Catch2
                GITURL https://github.com/catchorg/Catch2.git
                GITBRANCH devel
                SUBDIR "catch2"
            )
        endif()
    endif()
    message(STATUS "[Catch2] Setup complete")
    # other
    #include(CTest)
    #include(${PROJECT_SOURCE_DIR}/utils/spiner/Catch2/contrib/Catch.cmake)
    #add_subdirectory(test)
endif()

if(SINGULARITY_USE_HDF5)
    find_package(HDF5 COMPONENTS C HL QUIET)

    if (HDF5_FOUND)
        add_library(${PROJECT_NAME}::hdf5 INTERFACE IMPORTED)
        set_target_properties(${PROJECT_NAME}::hdf5
        PROPERTIES
            INTERFACE_LINK_LIBRARIES "${HDF5_LIBRARIES};${HDF5_HL_LIBRARIES}"
            INTERFACE_COMPILE_DEFINITIONS "SINGULARITY_USE_HDF5"
            INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIRS}"
        )
        if(HDF5_IS_PARALLEL)
            # TODO: this shouldn't be necessary, tho
            # i'm leaving it in for now
            if(NOT TARGET MPI::MPI_CXX)
                find_package(MPI COMPONENTS CXX QUIET)
                if(NOT MPI_FOUND)
                    message(FATAL_ERROR "The HDF5 package requires MPI")
                endif()
            endif()
        endif()
    endif()
endif()

if (SINGULARITY_USE_HDF5 AND SINGULARITY_BUILD_STELLARCOLLAPSE2SPINER)
   add_subdirectory(${PROJECT_SOURCE_DIR}/stellarcollapse2spiner)
   install(TARGETS stellarcollapse2spiner DESTINATION bin)
endif()