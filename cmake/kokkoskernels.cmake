include(CMakeDependentOption)

option(SINGULARITY_USE_KOKKOSKERNELS "Use kokkos kernels for linear algebra" OFF)
mark_as_advanced(SINGULARITY_USE_KOKKOSKERNELS)

if(SINGULARITY_USE_KOKKOSKERNELS)
    if(NOT SINGULARITY_USE_KOKKOS)
        message(FATAL_ERROR "Kokkoskernels requires Kokkos")
    endif()

    if(SINGULARITY_KOKKOSKERNELS_SUB_DIR)
        list(APPEND GPU_LINALG
            BLAS MKL LAPACK CUBLAS CUSPARSE MAMGA)
        list(APPEND GPU_MEMSET
            COMPLEX_DOUBLE COMPLEX_FLOAT DOUBLE FLOAT
            EXECSPACE_OPENMP EXECSPACE_SERIAL EXECSPACE_THREADS EXECSPACE_CUDA
            LAYOUTLEFT LAYOUTRIGHT MEMSPACE_HOSTSPACE
            OFFSET_INT OFFSET_SIZE_T
            ORDINAL_INT ORDINAL_INT
        )
        foreach(_KKTPL ${GPU_LINALG})
            set(KokkosKernels_ENABLE_TPL_${_KKTPL} OFF CACHE BOOL "" FORCE)
        endforeach()
        foreach(_KKMS ${GPU_MEMSET})
            set(KokkosKernels_INST_${_KKTPL} OFF CACHE BOOL "" FORCE)
        endforeach()

        add_subdirectory(${SINGULARITY_KOKKOSKERNELS_SUB_DIR}
                        ${PROJECT_BINARY_DIR}/utils/kokkos-kernels
        )
    else()
        if(SINGULARITY_KOKKOSKERNELS_INSTALL_DIR)
            set(KokkosKernels_ROOT ${SINGULARITY_KOKKOSKERNELS_INSTALL_DIR})
        endif()
        find_package(KokkosKernels REQUIRED)
        if(NOT KokkosKernels_FOUND)
            message(FATAL_ERROR "Kokkoskernels requested, but not found.\n \
            Use the CMake options SINGULARITY_KOKKOSKERNELS_SUB_DIR or SINGULARITY_KOKKOSKERNELS_INSTALL_DIR\n \
            if you know the source path.")
        endif()
    endif()
    list(APPEND TPL_LIBRARIES Kokkos::kokkoskernels)
    list(APPEND TPL_DEFINES SINGULARITY_USE_KOKKOSKERNELS)
endif()