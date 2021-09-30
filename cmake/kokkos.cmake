include(CMakeDependentOption)

option (SINGULARITY_USE_KOKKOS "Use Kokkos for portability" OFF)
mark_as_advanced(ENABLE_KOKKSKERNELS)

if(SINGULARITY_USE_KOKKOS)

    find_package(Kokkos REQUIRED)
    if(NOT Kokkos_FOUND)
        message(WARNING "Kokkos requested, but not found. Reverting to internal")
        set(Kokkos_ENABLE_CUDA ON CACHE BOOL "" FORCE)
        set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
        set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "" FORCE)
        set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE ON CACHE BOOL "" FORCE)
        add_subdirectory(${PROJECT_SOURCE_DIR}/utils/kokkos)
    endif()

    get_target_property(KOKKOS_COMPILE_OPTIONS Kokkos::kokkoscore
        INTERFACE_COMPILE_OPTIONS)

    list(APPEND TPL_LIBRARIES Kokkos::kokkos)
    list(APPEND TPL_DEFINES PORTABILITY_STRATEGY_KOKKOS)
    set(PORTABILITY_STRATEGY_KOKKOS ON)


endif()