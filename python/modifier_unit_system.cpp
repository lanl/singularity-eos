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
void unit_system_helper(py::module_ & m) {
  m.def("UnitSystem", [](T eos, eos_units_init::ThermalUnitsInit, const Real rho_unit, const Real sie_unit, const Real temp_unit){
    return UnitSystem<T>(std::move(eos), eos_units_init::ThermalUnitsInit(), rho_unit, sie_unit, temp_unit);
  }, py::arg("eos"), py::arg("units"), py::arg("rho_unit"), py::arg("sie_unit"), py::arg("temp_unit"));
  m.def("UnitSystem", [](T eos, eos_units_init::LengthTimeUnitsInit, const Real time_unit, const Real mass_unit, const Real length_unit, const Real temp_unit){
    return UnitSystem<T>(std::move(eos), eos_units_init::LengthTimeUnitsInit(), time_unit, mass_unit, length_unit, temp_unit);
  }, py::arg("eos"), py::arg("units"), py::arg("time_unit"), py::arg("mass_unit"), py::arg("length_unit"), py::arg("temp_unit"));
  m.def("UnitSystem", [](T eos, const Real rho_unit, const Real sie_unit, const Real temp_unit){
    return UnitSystem<T>(std::move(eos), rho_unit, sie_unit, temp_unit);
  }, py::arg("eos"), py::arg("rho_unit"), py::arg("sie_unit"), py::arg("temp_unit"));
}

template<typename T>
void unit_system_class_helper(py::module_ & m, std::string name) {
  eos_class<UnitSystem<T>>(m, std::string("UnitSystem") + name)
    .def(
      py::init<T, eos_units_init::ThermalUnitsInit, Real, Real, Real>(),
      py::arg("eos"), py::arg("units"), py::arg("rho_unit"), py::arg("sie_unit"), py::arg("temp_unit")
    )
    .def(
      py::init<T, eos_units_init::LengthTimeUnitsInit, Real, Real, Real, Real>(),
      py::arg("eos"), py::arg("units"), py::arg("time_unit"), py::arg("mass_unit"), py::arg("length_unit"), py::arg("temp_unit")
    )
    .def(
      py::init<T, Real, Real, Real>(),
      py::arg("eos"), py::arg("rho_unit"), py::arg("sie_unit"), py::arg("temp_unit")
    );
}

template<typename T>
void unit_system_eos_class(py::module_ & m, const char * name) {
  // define UnitSystem utility function
  unit_system_helper<T>(m);
  unit_system_helper<ShiftedEOS<T>>(m);
  unit_system_helper<ScaledEOS<T>>(m);
  unit_system_helper<ScaledEOS<ShiftedEOS<T>>>(m);

  unit_system_class_helper<T>(m, name);
  unit_system_class_helper<ShiftedEOS<T>>(m, std::string("Shifted") + name);
  unit_system_class_helper<ScaledEOS<T>>(m, std::string("Scaled") + name);
  unit_system_class_helper<ScaledEOS<ShiftedEOS<T>>>(m, std::string("ScaledShifted") + name);

  unit_system_helper<BilinearRampEOS<T>>(m);
  unit_system_helper<BilinearRampEOS<ShiftedEOS<T>>>(m);
  unit_system_helper<BilinearRampEOS<ScaledEOS<T>>>(m);
  unit_system_helper<BilinearRampEOS<ScaledEOS<ShiftedEOS<T>>>>(m);

  unit_system_class_helper<BilinearRampEOS<T>>(m, std::string("BilinearRamp") + name);
  unit_system_class_helper<BilinearRampEOS<ShiftedEOS<T>>>(m, std::string("BilinearRampShifted") + name);
  unit_system_class_helper<BilinearRampEOS<ScaledEOS<T>>>(m, std::string("BilinearRampScaled") + name);
  unit_system_class_helper<BilinearRampEOS<ScaledEOS<ShiftedEOS<T>>>>(m, std::string("BilinearRampScaledShifted") + name);

  unit_system_helper<RelativisticEOS<T>>(m);
  unit_system_helper<RelativisticEOS<ShiftedEOS<T>>>(m);
  unit_system_helper<RelativisticEOS<ScaledEOS<T>>>(m);
  unit_system_helper<RelativisticEOS<ScaledEOS<ShiftedEOS<T>>>>(m);

  unit_system_class_helper<RelativisticEOS<T>>(m, std::string("Relativistic") + name);
  unit_system_class_helper<RelativisticEOS<ShiftedEOS<T>>>(m, std::string("RelativisticShifted") + name);
  unit_system_class_helper<RelativisticEOS<ScaledEOS<T>>>(m, std::string("RelativisticScaled") + name);
  unit_system_class_helper<RelativisticEOS<ScaledEOS<ShiftedEOS<T>>>>(m, std::string("RelativisticScaledShifted") + name);

  unit_system_helper<RelativisticEOS<BilinearRampEOS<T>>>(m);
  unit_system_helper<RelativisticEOS<BilinearRampEOS<ShiftedEOS<T>>>>(m);
  unit_system_helper<RelativisticEOS<BilinearRampEOS<ScaledEOS<T>>>>(m);
  unit_system_helper<RelativisticEOS<BilinearRampEOS<ScaledEOS<ShiftedEOS<T>>>>>(m);

  unit_system_class_helper<RelativisticEOS<BilinearRampEOS<T>>>(m, std::string("RelativisticBilinearRamp") + name);
  unit_system_class_helper<RelativisticEOS<BilinearRampEOS<ShiftedEOS<T>>>>(m, std::string("RelativisticBilinearRampShifted") + name);
  unit_system_class_helper<RelativisticEOS<BilinearRampEOS<ScaledEOS<T>>>>(m, std::string("RelativisticBilinearRampScaled") + name);
  unit_system_class_helper<RelativisticEOS<BilinearRampEOS<ScaledEOS<ShiftedEOS<T>>>>>(m, std::string("RelativisticBilinearRampScaledShifted") + name);
}


void create_unit_system_eos_classes(py::module_ & m) {		       
  unit_system_eos_class<IdealGas>(m, "IdealGas");

#ifdef SPINER_USE_HDF
  unit_system_eos_class<SpinerEOSDependsRhoT>(m, "SpinerEOSDependsRhoT");
  unit_system_eos_class<SpinerEOSDependsRhoSie>(m, "SpinerEOSDependsRhoSie");
  unit_system_eos_class<StellarCollapse>(m, "StellarCollapse");
#endif
}
