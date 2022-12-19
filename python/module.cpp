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

PYBIND11_MODULE(singularity_eos, m) {
  py::class_<EOSState>(m, "EOSState")
    .def(py::init())
    .def_readwrite("density", &EOSState::density)
    .def_readwrite("specific_internal_energy", &EOSState::specific_internal_energy)
    .def_readwrite("pressure", &EOSState::pressure)
    .def_readwrite("temperature", &EOSState::temperature)
    .def_readwrite("specific_heat", &EOSState::specific_heat)
    .def_readwrite("bulk_modulus", &EOSState::bulk_modulus)
    .def("__repr__", &EOSState::to_string);

  eos_class<IdealGas>(m, "IdealGas")
    .def(py::init())
    .def(
      py::init<Real, Real>(),
      py::arg("gm1"), py::arg("Cv")
    );

  eos_class<Gruneisen>(m, "Gruneisen")
    .def(py::init())
    .def(
      py::init<Real, Real, Real, Real, Real, Real, Real, Real, Real, Real>(),
      py::arg("C0"), py::arg("s1"), py::arg("s2"), py::arg("s3"), py::arg("G0"),
      py::arg("b"), py::arg("rho0"), py::arg("T0"), py::arg("P0"), py::arg("Cv")
    );

  eos_class<JWL>(m, "JWL")
    .def(py::init())
    .def(
      py::init<Real, Real, Real, Real, Real, Real, Real>(),
      py::arg("A"), py::arg("B"), py::arg("R1"), py::arg("R2"),
      py::arg("w"), py::arg("rho0"), py::arg("Cv")
    );

  eos_class<DavisReactants>(m, "DavisReactants")
    .def(py::init())
    .def(
      py::init<Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real>(),
      py::arg("rho0"), py::arg("e0"), py::arg("P0"), py::arg("T0"),
      py::arg("A"), py::arg("B"), py::arg("C"), py::arg("G0"), py::arg("Z"),
      py::arg("alpha"), py::arg("Cv0")
    );

  eos_class<DavisProducts>(m, "DavisProducts")
    .def(py::init())
    .def(
      py::init<Real, Real, Real, Real, Real, Real, Real, Real>(),
      py::arg("a"), py::arg("b"), py::arg("k"), py::arg("n"), py::arg("vc"),
      py::arg("pc"), py::arg("Cv"), py::arg("E0")
    );

#ifdef SPINER_USE_HDF
  eos_class<SpinerEOSDependsRhoT>(m, "SpinerEOSDependsRhoT")
    .def(py::init())
    .def(py::init<const std::string&,int,bool>(), py::arg("filename"), py::arg("matid"), py::arg("reproduciblity_mode")=false)
    .def(py::init<const std::string&,const std::string&,bool>(), py::arg("filename"), py::arg("materialName"), py::arg("reproduciblity_mode")=false)
    .def_property_readonly("matid", &SpinerEOSDependsRhoT::matid)
    .def_property_readonly("lRhoOffset", &SpinerEOSDependsRhoT::lRhoOffset)
    .def_property_readonly("lTOffset", &SpinerEOSDependsRhoT::lTOffset)
    .def_property_readonly("rhoMin", &SpinerEOSDependsRhoT::rhoMin)
    .def_property_readonly("rhoMax", &SpinerEOSDependsRhoT::rhoMax);

  eos_class<SpinerEOSDependsRhoSie>(m, "SpinerEOSDependsRhoSie")
    .def(py::init())
    .def(py::init<const std::string&,int,bool>(), py::arg("filename"), py::arg("matid"), py::arg("reproduciblity_mode")=false)
    .def(py::init<const std::string&,const std::string&,bool>(), py::arg("filename"), py::arg("materialName"), py::arg("reproduciblity_mode")=false)
    .def_property_readonly("matid", &SpinerEOSDependsRhoSie::matid)
    .def_property_readonly("lRhoOffset", &SpinerEOSDependsRhoSie::lRhoOffset)
    .def_property_readonly("lTOffset", &SpinerEOSDependsRhoSie::lTOffset)
    .def_property_readonly("lEOffset", &SpinerEOSDependsRhoSie::lEOffset)
    .def_property_readonly("rhoMin", &SpinerEOSDependsRhoSie::rhoMin)
    .def_property_readonly("rhoMax", &SpinerEOSDependsRhoSie::rhoMax)
    .def_property_readonly("TMin", &SpinerEOSDependsRhoSie::TMin)
    .def_property_readonly("TMax", &SpinerEOSDependsRhoSie::TMax)
    .def_property_readonly("sieMin", &SpinerEOSDependsRhoSie::sieMin)
    .def_property_readonly("sieMax", &SpinerEOSDependsRhoSie::sieMax);

  eos_class<StellarCollapse>(m, "StellarCollapse")
    .def(py::init())
    .def(py::init<const std::string&, bool, bool>(), py::arg("filename"), py::arg("use_sp5")=false, py::arg("filter_bmod")=true)
    .def("Save", &StellarCollapse::Save, py::arg("filename"))
    .def_property_readonly("lRhoOffset", &StellarCollapse::lRhoOffset)
    .def_property_readonly("lTOffset", &StellarCollapse::lTOffset)
    .def_property_readonly("lEOffset", &StellarCollapse::lEOffset)
    .def_property_readonly("lRhoMin", &StellarCollapse::lRhoMin)
    .def_property_readonly("lRhoMax", &StellarCollapse::lRhoMax)
    .def_property_readonly("rhoMin", &StellarCollapse::rhoMin)
    .def_property_readonly("rhoMax", &StellarCollapse::rhoMax)
    .def_property_readonly("lTMin", &StellarCollapse::lTMin)
    .def_property_readonly("lTMax", &StellarCollapse::lTMax)
    .def_property_readonly("TMin", &StellarCollapse::TMin)
    .def_property_readonly("TMax", &StellarCollapse::TMax)
    .def_property_readonly("YeMin", &StellarCollapse::YeMin)
    .def_property_readonly("YeMax", &StellarCollapse::YeMax)
    .def_property_readonly("sieMin", &StellarCollapse::sieMin)
    .def_property_readonly("sieMax", &StellarCollapse::sieMax);
#endif

#ifdef SINGULARITY_USE_EOSPAC
  constexpr bool use_scratch = true;
  eos_class<EOSPAC, use_scratch>(m, "EOSPAC")
    .def(py::init())
    .def(py::init<int, bool>(), py::arg("matid"), py::arg("invert_at_setup")=false)

#endif

  create_shifted_eos_classes(m);
  create_scaled_eos_classes(m);
  create_bilinear_ramp_eos_classes(m);
  create_relativistic_eos_classes(m);
  create_unit_system_eos_classes(m);
  
  m.def("pAlpha2BilinearRampParams", [](const IdealGas &eos, const Real alpha0, const Real Pe, const Real Pc){
    Real r0, a, b, c;
    pAlpha2BilinearRampParams(eos, alpha0, Pe, Pc, r0, a, b, c);
    return py::make_tuple(r0, a, b, c);
  }, py::arg("eos"), py::arg("alpha0"), py::arg("Pe"), py::arg("Pc"));

  py::module thermalqs = m.def_submodule("thermalqs");
  thermalqs.attr("none") = pybind11::int_(thermalqs::none);
  thermalqs.attr("density") = pybind11::int_(thermalqs::density);
  thermalqs.attr("specific_internal_energy") = pybind11::int_(thermalqs::specific_internal_energy);
  thermalqs.attr("pressure") = pybind11::int_(thermalqs::pressure);
  thermalqs.attr("temperature") = pybind11::int_(thermalqs::temperature);
  thermalqs.attr("specific_heat") = pybind11::int_(thermalqs::specific_heat);
  thermalqs.attr("bulk_modulus") = pybind11::int_(thermalqs::bulk_modulus);
  thermalqs.attr("all_values") = pybind11::int_(thermalqs::all_values);

  py::module eos_units = m.def_submodule("eos_units");
  py::class_<eos_units_init::ThermalUnitsInit>(eos_units, "_ThermalUnits");
  py::class_<eos_units_init::LengthTimeUnitsInit>(eos_units, "_LengthTimeUnits");
  eos_units.attr("ThermalUnits") = eos_units_init::thermal_units_init_tag;
  eos_units.attr("LengthTimeUnits") = eos_units_init::length_time_units_init_tag;

  m.doc() = "Singularity EOS Python Bindings";
}
