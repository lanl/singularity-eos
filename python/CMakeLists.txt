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
