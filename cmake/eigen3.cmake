option(SINGULARITY_USE_EIGEN "Use Eigen linear algebra" OFF)

if(SINGULARITY_USE_EIGEN)
    if(SINGULARITY_USE_KOKKOSKERNELS)
        message(FATAL_ERROR "Cannot use Eigen and Kokkoskernels together")
    endif()
    find_package(Eigen3 REQUIRED)
    list(APPEND TPL_LIBRARIES Eigen3::Eigen)
endif()