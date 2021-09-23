option (SINGULARITY_USE_SINGLE_LOGS "Use single precision logs. Can harm accuracy." OFF)

set(SINGULARITY_USE_FMATH_ORDERS 7 5 0)
if(NOT SINGULARITY_USE_FMATH_ORDER)
    list(GET SINGULARITY_USE_FMATH_ORDERS 0 SINGULARITY_USE_FMATH_ORDER)
endif()

set(SINGULARITY_USE_FMATH_ORDER "${SINGULARITY_USE_FMATH_ORDER}" CACHE STRING "order interpolant for fast logs. Default is 7th order")
set_property(CACHE SINGULARITY_USE_FMATH_ORDER PROPERTY STRINGS ${SINGULARITY_FMATH_ORDERS})

mark_as_advanced(SINGULARITY_USE_FMATH_ORDER)

if (SINGULARITY_USE_SINGLE_LOGS)
    list(APPEND TPL_DEFINES SINGULARITY_USE_SINGLE_LOGS)
endif()

list(APPEND TPL_DEFINES SINGULARITY_USE_FMATH_ORDER=${SINGULARITY_USE_FMATH_ORDER})

include(GNUInstallDirs)

# kinda winging it here
function(add_small_library libname)
    add_library(${libname} INTERFACE IMPORTED)

    target_include_directories(${libname} INTERFACE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utils>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${libname}>
    )

    install(DIRECTORY ${PROJECT_SOURCE_DIR}/utils/${libname}
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${libname})
endfunction()

add_small_library(fast-math)
add_small_library(root-finding-1d)
add_small_library(spiner)