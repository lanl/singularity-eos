# ------------------------------------------------------------------------------#
# © 2021-2026. Triad National Security, LLC. All rights reserved.  This program
# was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
# National Laboratory (LANL), which is operated by Triad National Security, LLC
# for the U.S.  Department of Energy/National Nuclear Security Administration.
# All rights in the program are reserved by Triad National Security, LLC, and
# the U.S. Department of Energy/National Nuclear Security Administration. The
# Government is granted for itself and others acting on its behalf a
# nonexclusive, paid-up, irrevocable worldwide license in this material to
# reproduce, prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.
# ------------------------------------------------------------------------------#

pybind11_add_module(singularity_eos
  module.cpp
  modifier_shifted.cpp
  modifier_scaled.cpp
  modifier_bilinear_ramp.cpp
)
target_link_libraries(singularity_eos PRIVATE singularity-eos_Interface)

include(GNUInstallDirs)

set(PYTHON_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages")

install(TARGETS singularity_eos
  COMPONENT python
  LIBRARY DESTINATION "${PYTHON_INSTALL_DIR}")
