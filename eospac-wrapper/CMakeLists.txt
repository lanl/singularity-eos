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

option(SINGULARITY_EOSPAC_SKIP_EXTRAP "Skip extrapolation checks" OFF)

# Fix for #228, tho should probably be part of singularity-eos lib
list(APPEND EOSWRAPPER_HEADERS eospac_wrapper.hpp)
list(APPEND EOSWRAPPER_SRCS eospac_wrapper.cpp)

add_library(eospac-wrapper
  ${EOSWRAPPER_SRCS} ${EOSWRAPPER_HEADERS}
)

target_include_directories(eospac-wrapper
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

if (SINGULARITY_EOSPAC_SKIP_EXTRAP)
  target_compile_definitions(eospac-wrapper
    PRIVATE
      SINGULARITY_EOSPAC_SKIP_EXTRAP
  )
endif() #SINGULARITY_EOSPAC_SKIP_EXTRAP

target_link_libraries(eospac-wrapper 
  PUBLIC
    EOSPAC::eospac
)

# install header file
install(
  FILES ${EOSWRAPPER_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/eospac-wrapper/
)

install(
  TARGETS eospac-wrapper EXPORT eospac-wrapperTargets DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(
  EXPORT eospac-wrapperTargets
  FILE eospac-wrapperTargets.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/singularity-eos
  )

export(
  EXPORT eospac-wrapperTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/cmake/eospac-wrapperTargets.cmake
  )
