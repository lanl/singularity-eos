include(CMakeDependentOption)
option(SINGULARITY_USE_CUDA "Enable CUDA" OFF)
cmake_dependent_option(SINGULARITY_HIDE_MORE_WARNINGS "Disable flags the suppress warnings" ON "SINGULARITY_USE_CUDA" OFF)
cmake_dependent_option(SINGULARITY_BETTER_DEBUG_FLAGS "Disable flags the suppress warnings" ON "SINGULARITY_USE_CUDA; CMAKE_BUILD_TYPE MATCHES DEBUG" OFF)

mark_as_advanced(SINGULARITY_USE_CUDA)

if(SINGULARITY_USE_CUDA)
    find_package(CUDAToolkit REQUIRED)
    if(NOT CUDAToolkit_FOUND)
        message(FATAL_ERROR "CUDA requested, but not found")
    endif()

    list(APPEND XTRA_BUILD_FLAGS "$<${build_release}:-use_fast_math>")
    list(APPEND XTRA_BUILD_FLAGS "$<${with_cxx}:--expt-relaxed-constexpr>")
    if(NOT SINGULARITY_HIDE_MORE_WARNINGS)
        list(APPEND XTRA_BUILD_FLAGS "$<${with_cxx}:-Xcudafe>")
        list(APPEND XTRA_BUILD_FLAGS "$<${with_cxx}:--diag_suppress=esa_on_defaulted_function_ignored>")
    endif()
    if(SINGULARITY_BETTER_DEBUG_FLAGS)
        list(APPEND XTRA_BUILD_FLAGS "$<AND<${build_debug},${with_cxx}>:-G;-lineinfo")
    endif()

endif()