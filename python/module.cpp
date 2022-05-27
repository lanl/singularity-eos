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
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <singularity-eos/eos/eos.hpp>
#include <map>
#include <string>
#include <ostream>
#include <sstream>
#include <limits>

namespace py = pybind11;
using namespace singularity;

// Helper function to convert lambda numpy array to double* buffer
// With std::optional we would add support for a default value of lambda=None
template<typename T, PORTABLE_FUNCTION Real(T::*Func)(const Real, const Real, Real*) const>
Real two_params(const T& self, const Real a, const Real b, py::array_t<Real> lambda) {
  return (self.*Func)(a, b, lambda.mutable_data());
}

class LambdaHelper {
  py::array_t<Real> & lambdas;
public:
  LambdaHelper(py::array_t<Real> & lambdas) : lambdas(lambdas) {}
  Real * operator[](const int i) const {
    return lambdas.mutable_data(i,0);
  }
};

class NoLambdaHelper {
public:
  Real * operator[](const int i) const {
    return nullptr;
  }
};

// so far didn't find a good way of working with template member function pointers
// to generalize this without the preprocessor.
#define EOS_VEC_FUNC_TMPL(func, a, b, out)                                        \
template<typename T>                                                              \
void func(const T & self, py::array_t<Real> a, py::array_t<Real> b,                \
          py::array_t<Real> out, const int num, py::array_t<Real> lambdas){       \
  py::buffer_info lambdas_info = lambdas.request();                               \
  if (lambdas_info.ndim != 2)                                                     \
      throw std::runtime_error("lambdas dimension must be 2!");                   \
                                                                                  \
  if(lambdas_info.shape[1] > 0) {                                                 \
    self.func(a.data(), b.data(), out.mutable_data(), num, LambdaHelper(lambdas)); \
  } else {                                                                        \
    self.func(a.data(), b.data(), out.mutable_data(), num, NoLambdaHelper());      \
  }                                                                               \
}

EOS_VEC_FUNC_TMPL(TemperatureFromDensityInternalEnergy, rhos, sies, temperatures)
EOS_VEC_FUNC_TMPL(InternalEnergyFromDensityTemperature, rhos, temperatures, sies)
EOS_VEC_FUNC_TMPL(PressureFromDensityTemperature, rhos, temperatures, pressures)
EOS_VEC_FUNC_TMPL(PressureFromDensityInternalEnergy, rhos, sies, pressures)
EOS_VEC_FUNC_TMPL(SpecificHeatFromDensityTemperature, rhos, temperatures, cvs)
EOS_VEC_FUNC_TMPL(SpecificHeatFromDensityInternalEnergy, rhos, sies, cvs)
EOS_VEC_FUNC_TMPL(BulkModulusFromDensityTemperature, rhos, temperatures, bmods)
EOS_VEC_FUNC_TMPL(BulkModulusFromDensityInternalEnergy, rhos, sies, bmods)
EOS_VEC_FUNC_TMPL(GruneisenParamFromDensityTemperature, rhos, temperatures, gm1s)
EOS_VEC_FUNC_TMPL(GruneisenParamFromDensityInternalEnergy, rhos, sies, gm1s)

struct EOSState {
  Real density;
  Real specific_internal_energy;
  Real pressure;
  Real temperature;
  Real specific_heat;
  Real bulk_modulus;

  std::string to_string() const {
    std::stringstream ss;
    ss << "density: " << density << std::endl;
    ss << "specific_internal_energy: " << specific_internal_energy << std::endl;
    ss << "pressure: " << pressure << std::endl;
    ss << "temperature: " << temperature << std::endl;
    ss << "specific_heat: " << specific_heat << std::endl;
    ss << "bulk_modulus: " << bulk_modulus;
    return ss.str();
  }
};

template<typename T>
py::class_<T> eos_class(py::module_ & m, const char * name) {
  return py::class_<T>(m, name)
    .def("TemperatureFromDensityInternalEnergy", &two_params<T, &T::TemperatureFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"), py::arg("lmbda"))
    .def("InternalEnergyFromDensityTemperature", &two_params<T, &T::InternalEnergyFromDensityTemperature>, py::arg("rho"), py::arg("temperature"), py::arg("lmbda"))
    .def("PressureFromDensityTemperature", &two_params<T, &T::PressureFromDensityTemperature>, py::arg("rho"), py::arg("temperature"), py::arg("lmbda"))
    .def("PressureFromDensityInternalEnergy", &two_params<T, &T::PressureFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"), py::arg("lmbda"))
    .def("SpecificHeatFromDensityTemperature", &two_params<T, &T::SpecificHeatFromDensityTemperature>, py::arg("rho"), py::arg("temperature"), py::arg("lmbda"))
    .def("SpecificHeatFromDensityInternalEnergy", &two_params<T, &T::SpecificHeatFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"), py::arg("lmbda"))
    .def("BulkModulusFromDensityTemperature", &two_params<T, &T::BulkModulusFromDensityTemperature>, py::arg("rho"), py::arg("temperature"), py::arg("lmbda"))
    .def("BulkModulusFromDensityInternalEnergy", &two_params<T, &T::BulkModulusFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"), py::arg("lmbda"))
    .def("GruneisenParamFromDensityTemperature", &two_params<T, &T::GruneisenParamFromDensityTemperature>, py::arg("rho"), py::arg("temperature"), py::arg("lmbda"))
    .def("GruneisenParamFromDensityInternalEnergy", &two_params<T, &T::GruneisenParamFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"), py::arg("lmbda"))


    // TODO .def("GetOnDevice")
    .def("FillEos", [](const T & self, const py::kwargs& kwargs) {
      unsigned long output = thermalqs::none;
      Real rho, temp, sie, press, cv, bmod;
      rho = temp = sie = press = cv = bmod = std::numeric_limits<double>::quiet_NaN();
      std::map<std::string, std::pair<unsigned long,Real*>> param_mapping {
        {"rho", {thermalqs::density, &rho}},
        {"sie", {thermalqs::specific_internal_energy, &sie}},
        {"press", {thermalqs::pressure, &press}},
        {"temp", {thermalqs::temperature, &temp}},
        {"cv", {thermalqs::specific_heat, &cv}},
        {"bmod", {thermalqs::bulk_modulus, &bmod}},
      };

      for(auto const & it : param_mapping) {
        auto param = it.first.c_str();
        if(kwargs.contains(param)) {
          auto param_bitmask = it.second.first;
          auto param_ptr = it.second.second;
          output |= param_bitmask;
          *param_ptr = kwargs[param].cast<Real>();
        }
      }

      // override auto-detection
      if(kwargs.contains("output")) {
        output = kwargs["output"].cast<unsigned long>();
      }

      if(kwargs.contains("lmbda")) {
        auto lambda = kwargs["lmbda"].cast<py::array_t<Real>>();
        self.FillEos(rho, temp, sie, press, cv, bmod, output, lambda.mutable_data());
      } else {
        self.FillEos(rho, temp, sie, press, cv, bmod, output, nullptr);
      }
      return EOSState {rho, sie, press, temp, cv, bmod};
    })
    // TODO .def("ValuesAtReferenceState")

    // Generic functions provided by the base class. These contain e.g. the vector
    // overloads that use the scalar versions declared here
    .def("TemperatureFromDensityInternalEnergy", &TemperatureFromDensityInternalEnergy<T>, py::arg("rhos"), py::arg("sies"), py::arg("temperatures"), py::arg("num"), py::arg("lmbdas"))
    .def("InternalEnergyFromDensityTemperature", &InternalEnergyFromDensityTemperature<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("sies"), py::arg("num"), py::arg("lmbdas"))
    .def("PressureFromDensityTemperature", &PressureFromDensityTemperature<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("pressures"), py::arg("num"), py::arg("lmbdas"))
    .def("PressureFromDensityInternalEnergy", &PressureFromDensityInternalEnergy<T>, py::arg("rhos"), py::arg("sies"), py::arg("pressures"), py::arg("num"), py::arg("lmbdas"))
    .def("SpecificHeatFromDensityTemperature", &SpecificHeatFromDensityTemperature<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("cvs"), py::arg("num"), py::arg("lmbdas"))
    .def("SpecificHeatFromDensityInternalEnergy", &SpecificHeatFromDensityInternalEnergy<T>, py::arg("rhos"), py::arg("sies"), py::arg("cvs"), py::arg("num"), py::arg("lmbdas"))
    .def("BulkModulusFromDensityTemperature", &BulkModulusFromDensityTemperature<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("bmods"), py::arg("num"), py::arg("lmbdas"))
    .def("BulkModulusFromDensityInternalEnergy", &BulkModulusFromDensityInternalEnergy<T>, py::arg("rhos"), py::arg("sies"), py::arg("bmods"), py::arg("num"), py::arg("lmbdas"))
    .def("GruneisenParamFromDensityTemperature", &GruneisenParamFromDensityTemperature<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("gm1s"), py::arg("num"), py::arg("lmbdas"))
    .def("GruneisenParamFromDensityInternalEnergy", &GruneisenParamFromDensityInternalEnergy<T>, py::arg("rhos"), py::arg("sies"), py::arg("gm1s"), py::arg("num"), py::arg("lmbdas"))

    .def("nlambda", &T::nlambda)
    .def_static("PreferredInput", &T::PreferredInput)
    .def("PrintParams", &T::PrintParams)
    // TODO .def("DensityEnergyFromPressureTemperature") reference semantics
    .def("Finalize", &T::Finalize)
    .def_static("EosType", &T::EosType);
}

PYBIND11_MODULE(singularity_eos, m) {
  py::class_<EOSState>(m, "EOSState")
    .def_readwrite("density", &EOSState::density)
    .def_readwrite("specific_internal_energy", &EOSState::specific_internal_energy)
    .def_readwrite("pressure", &EOSState::pressure)
    .def_readwrite("temperature", &EOSState::temperature)
    .def_readwrite("specific_heat", &EOSState::specific_heat)
    .def_readwrite("bulk_modulus", &EOSState::bulk_modulus)
    .def("__str__", &EOSState::to_string);

  eos_class<IdealGas>(m, "IdealGas")
    .def(
      py::init<Real, Real>(),
      py::arg("gm1"), py::arg("Cv")
    );

  eos_class<Gruneisen>(m, "Gruneisen")
    .def(
      py::init<Real, Real, Real, Real, Real, Real, Real, Real, Real, Real>(),
      py::arg("C0"), py::arg("s1"), py::arg("s2"), py::arg("s3"), py::arg("G0"),
      py::arg("b"), py::arg("rho0"), py::arg("T0"), py::arg("P0"), py::arg("Cv")
    );

  eos_class<JWL>(m, "JWL")
    .def(
      py::init<Real, Real, Real, Real, Real, Real, Real>(),
      py::arg("A"), py::arg("B"), py::arg("R1"), py::arg("R2"),
      py::arg("w"), py::arg("rho0"), py::arg("Cv")
    );

  eos_class<DavisReactants>(m, "DavisReactants")
    .def(
      py::init<Real, Real, Real, Real, Real, Real, Real, Real, Real, Real, Real>(),
      py::arg("rho0"), py::arg("e0"), py::arg("P0"), py::arg("T0"),
      py::arg("A"), py::arg("B"), py::arg("C"), py::arg("G0"), py::arg("Z"),
      py::arg("alpha"), py::arg("Cv0")
    );

  eos_class<DavisProducts>(m, "DavisProducts")
    .def(
      py::init<Real, Real, Real, Real, Real, Real, Real, Real>(),
      py::arg("a"), py::arg("b"), py::arg("k"), py::arg("n"), py::arg("vc"),
      py::arg("pc"), py::arg("Cv"), py::arg("E0")
    );

  py::module thermalqs = m.def_submodule("thermalqs");
  thermalqs.attr("none") = pybind11::int_(thermalqs::none);
  thermalqs.attr("density") = pybind11::int_(thermalqs::density);
  thermalqs.attr("specific_internal_energy") = pybind11::int_(thermalqs::specific_internal_energy);
  thermalqs.attr("pressure") = pybind11::int_(thermalqs::pressure);
  thermalqs.attr("temperature") = pybind11::int_(thermalqs::temperature);
  thermalqs.attr("specific_heat") = pybind11::int_(thermalqs::specific_heat);
  thermalqs.attr("bulk_modulus") = pybind11::int_(thermalqs::bulk_modulus);
  thermalqs.attr("all_values") = pybind11::int_(thermalqs::all_values);

  m.doc() = "Singularity EOS Python Bindings";
}
