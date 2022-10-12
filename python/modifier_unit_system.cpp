//------------------------------------------------------------------------------
// © 2021-2022. Triad National Security, LLC. All rights reserved.  This
// program was produced under U.S. Government contract 89233218CNA000001
// for Los Alamos National Laboratory (LANL), which is operated by Triad
// National Security, LLC for the U.S.  Department of Energy/National
// Nuclear Security Administration. All rights in the program are
// reserved by Triad National Security, LLC, and the U.S. Department of
// Energy/National Nuclear Security Administration. The Government is
// granted for itself and others acting on its behalf a nonexclusive,
// paid-up, irrevocable worldwide license in this material to reproduce,
// prepare derivative works, distribute copies to the public, perform
// publicly and display publicly, and to permit others to do so.
//------------------------------------------------------------------------------
// clang-format off
#include "module.hpp"

template<typename T>
void unit_system_eos_class(py::module_ & m) {
  // define UnitSystem utility function
  m.def("UnitSystem", [](typename T::BaseType eos, eos_units_init::ThermalUnitsInit, const Real rho_unit, const Real sie_unit, const Real temp_unit){
    return T(std::move(eos), eos_units_init::ThermalUnitsInit(), rho_unit, sie_unit, temp_unit);
  }, py::arg("eos"), py::arg("units"), py::arg("rho_unit"), py::arg("sie_unit"), py::arg("temp_unit"));
  m.def("UnitSystem", [](typename T::BaseType eos, eos_units_init::LengthTimeUnitsInit, const Real time_unit, const Real mass_unit, const Real length_unit, const Real temp_unit){
    return T(std::move(eos), eos_units_init::LengthTimeUnitsInit(), time_unit, mass_unit, length_unit, temp_unit);
  }, py::arg("eos"), py::arg("units"), py::arg("time_unit"), py::arg("mass_unit"), py::arg("length_unit"), py::arg("temp_unit"));
  m.def("UnitSystem", [](typename T::BaseType eos, const Real rho_unit, const Real sie_unit, const Real temp_unit){
    return T(std::move(eos), rho_unit, sie_unit, temp_unit);
  }, py::arg("eos"), py::arg("rho_unit"), py::arg("sie_unit"), py::arg("temp_unit"));

  // define UnitSystem class
  eos_class<T>(m, T::EosPyType())
    .def(
      py::init<typename T::BaseType, eos_units_init::ThermalUnitsInit, Real, Real, Real>(),
      py::arg("eos"), py::arg("units"), py::arg("rho_unit"), py::arg("sie_unit"), py::arg("temp_unit")
    )
    .def(
      py::init<typename T::BaseType, eos_units_init::LengthTimeUnitsInit, Real, Real, Real, Real>(),
      py::arg("eos"), py::arg("units"), py::arg("time_unit"), py::arg("mass_unit"), py::arg("length_unit"), py::arg("temp_unit")
    )
    .def(
      py::init<typename T::BaseType, Real, Real, Real>(),
      py::arg("eos"), py::arg("rho_unit"), py::arg("sie_unit"), py::arg("temp_unit")
    );
}

template<typename... Ts>
void create_unit_system(py::module_ &m, tl<Ts...>) {
  // C++14 workaround, since we don't have C++17 fold expressions
  auto l = { (unit_system_eos_class<Ts>(m), 0)... };
}

void create_unit_system_eos_classes(py::module_ & m) {
  create_unit_system(m, singularity::unit_system);
}
