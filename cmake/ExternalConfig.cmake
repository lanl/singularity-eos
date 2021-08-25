#------------------------------------------------------------------------------#
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
#------------------------------------------------------------------------------#

# for submodules
include(FetchContent)

# macro to automate dependency resolution
# arguments:
#   USE_INTERNAL    : use in-tree version if not immediately found
#   DEP             : the dependency name (as in `find_package(${DEP}`))
#   TARGET          : the cmake target (if any) defined by the dependency
#   GITURL          : the git url of the dependency, for remote builds
#   GITTAG          : the tag to fetch
#   SUBDIR          : the internal directory path, if USE_INTERNAL
# ----------
# this macro is ordered as:
#   1.) check if target exists.
#       |->if so, assume everything is okay and leave
#   2.) check if `find_package` succeeds
#       |->if so, assume `find_package` did it's job, and leave
#   3.) no trace of package in system
#       |-> if USE_INTERNAL, select as SUBDIR the source directory
#       |-> else select GITURL:GITTAG as the source directory and fetch it
#   4.) (implicit) users are expected to have an environment that isn't pathological
#
macro(singularity_select_dep)
    set(options USE_INTERNAL)
    set(oneValueArgs DEP TARGET GITURL GITTAG SUBDIR)
    set(multiValueArgs OTHERS)
    cmake_parse_arguments(SD "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    if(NOT SD_TARGET)
        singularity_message(STATUS "[${SD_DEP}] Target not defined, attempting to find in environment")
        find_package(${SD_DEP} QUIET)
        if(NOT ${SD_DEP}_FOUND)
            singularity_message(STATUS "[${SD_DEP}] find_package couldn't locate")
            if(SD_USE_INTERNAL)
                set(UTIL_PKG "${PROJECT_SOURCE_DIR}/utils/${SD_SUBDIR}")
                singularity_message(STATUS "[${SD_DEP}] USE_INTERNAL flag set, using ${UTIL_PKG} ")
                FetchContent_Declare(
                    ${SD_DEP}
                    FETCHCONTENT_SOURCE_DIR ${UTIL_PKG}
                )
                set(${SD_DEP}_IS_INTERNAL "1")
            else()
                singularity_message(STATUS "[${SD_DEP}] fetching source from ${SD_GITURL}")
                FetchContent_Declare(
                    ${SD_DEP}
                    GIT_REPOSITORY ${SD_GITURL}
                    GIT_TAG        ${SD_GITTAG}
                )
                set(${SD_DEP}_IS_REMOTE "1")
            endif()
            FetchContent_MakeAvailable(${SD_DEP})
        else()
            singularity_message(STATUS "[${SD_DEP}] found with find_package")
            set(${SD_DEP}_IS_EXTERNAL "1")
        endif()
    else()
        singularity_message(STATUS "[${SD_DEP}] Target already defined, skipping package search")
    endif()

endmacro()

#######################################
## internal configurations ############
#######################################
## these are set before calls to    ###
## `singularity_select_dep` to      ###
## prepare to a possible internal   ###
## source directory                 ###
#######################################

macro(singularity_config_internal_kokkos WITH_CUDA)
    set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
    if(WITH_CUDA)
        set(Kokkos_ENABLE_CUDA ON CACHE BOOL "" FORCE)
        set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "" FORCE)
        set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE ON CACHE BOOL "" FORCE)
    endif()
endmacro()

macro(singularity_config_internal_kokkoskernels)
    # Disable TPLs
    set(KokkosKernels_ENABLE_TPL_BLAS OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_ENABLE_TPL_MKL OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_ENABLE_TPL_LAPACK OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_ENABLE_TPL_CUBLAS OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_ENABLE_TPL_CUSPARSE OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_ENABLE_TPL_MAGMA OFF CACHE BOOL "" FORCE)
    # Disable ETIs
    set(KokkosKernels_INST_COMPLEX_DOUBLE OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_INST_COMPLEX_FLOAT OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_INST_DOUBLE OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_INST_EXECSPACE_OPENMP OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_INST_EXECSPACE_SERIAL OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_INST_EXECSPACE_THREADS OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_INST_EXECSPACE_CUDA OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_INST_FLOAT OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_INST_LAYOUTLEFT OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_INST_LAYOUTRIGHT OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_INST_MEMSPACE_HOSTSPACE OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_INST_OFFSET_INT OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_INST_OFFSET_SIZE_T OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_INST_ORDINAL_INT OFF CACHE BOOL "" FORCE)
    set(KokkosKernels_INST_ORDINAL_INT OFF CACHE BOOL "" FORCE)
endmacro()

macro(singularity_config_internal_json)
    set(JSON_BuildTests OFF CACHE INTERNAL "")
    set(JSON_Install OFF CACHE INTERNAL "")
endmacro()