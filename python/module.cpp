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

namespace py = pybind11;
using namespace singularity;

// Helper function to convert lambda numpy array to double* buffer
// With std::optional we would add support for a default value of lambda=None
template<typename T, PORTABLE_FUNCTION Real(T::*Func)(const Real, const Real, Real*) const>
Real two_params(const T& obj, const Real a, const Real b, py::array_t<Real> lambda) {
  return (obj.*Func)(a, b, lambda.mutable_data());
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
void func(const T & obj, py::array_t<Real> a, py::array_t<Real> b,                \
          py::array_t<Real> out, const int num, py::array_t<Real> lambdas){       \
  py::buffer_info lambdas_info = lambdas.request();                               \
  if (lambdas_info.ndim != 2)                                                     \
      throw std::runtime_error("lambdas dimension must be 2!");                   \
                                                                                  \
  if(lambdas_info.shape[1] > 0) {                                                 \
    obj.func(a.data(), b.data(), out.mutable_data(), num, LambdaHelper(lambdas)); \
  } else {                                                                        \
    obj.func(a.data(), b.data(), out.mutable_data(), num, NoLambdaHelper());      \
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

template<typename T>
py::class_<T> eos_class(py::module_ & m, const char * name) {

  return py::class_<T>(m, name)
    .def("TemperatureFromDensityInternalEnergy", &two_params<T, &T::TemperatureFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"), py::arg("lambda"))
    .def("InternalEnergyFromDensityTemperature", &two_params<T, &T::InternalEnergyFromDensityTemperature>, py::arg("rho"), py::arg("temperature"), py::arg("lambda"))
    .def("PressureFromDensityTemperature", &two_params<T, &T::PressureFromDensityTemperature>, py::arg("rho"), py::arg("temperature"), py::arg("lambda"))
    .def("PressureFromDensityInternalEnergy", &two_params<T, &T::PressureFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"), py::arg("lambda"))
    .def("SpecificHeatFromDensityTemperature", &two_params<T, &T::SpecificHeatFromDensityTemperature>, py::arg("rho"), py::arg("temperature"), py::arg("lambda"))
    .def("SpecificHeatFromDensityInternalEnergy", &two_params<T, &T::SpecificHeatFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"), py::arg("lambda"))
    .def("BulkModulusFromDensityTemperature", &two_params<T, &T::BulkModulusFromDensityTemperature>, py::arg("rho"), py::arg("temperature"), py::arg("lambda"))
    .def("BulkModulusFromDensityInternalEnergy", &two_params<T, &T::BulkModulusFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"), py::arg("lambda"))
    .def("GruneisenParamFromDensityTemperature", &two_params<T, &T::GruneisenParamFromDensityTemperature>, py::arg("rho"), py::arg("temperature"), py::arg("lambda"))
    .def("GruneisenParamFromDensityInternalEnergy", &two_params<T, &T::GruneisenParamFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"), py::arg("lambda"))


    // TODO .def("GetOnDevice")
    // TODO .def("FillEos")  what are the call semantics, since it uses references? return dict?
    // TODO .def("ValuesAtReferenceState")

    // Generic functions provided by the base class. These contain e.g. the vector
    // overloads that use the scalar versions declared here
    .def("TemperatureFromDensityInternalEnergy", &TemperatureFromDensityInternalEnergy<T>, py::arg("rhos"), py::arg("sies"), py::arg("temperatures"), py::arg("num"), py::arg("lambdas"))
    .def("InternalEnergyFromDensityTemperature", &InternalEnergyFromDensityTemperature<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("sies"), py::arg("num"), py::arg("lambdas"))
    .def("PressureFromDensityTemperature", &PressureFromDensityTemperature<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("pressures"), py::arg("num"), py::arg("lambdas"))
    .def("PressureFromDensityInternalEnergy", &PressureFromDensityInternalEnergy<T>, py::arg("rhos"), py::arg("sies"), py::arg("pressures"), py::arg("num"), py::arg("lambdas"))
    .def("SpecificHeatFromDensityTemperature", &SpecificHeatFromDensityTemperature<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("cvs"), py::arg("num"), py::arg("lambdas"))
    .def("SpecificHeatFromDensityInternalEnergy", &SpecificHeatFromDensityInternalEnergy<T>, py::arg("rhos"), py::arg("sies"), py::arg("cvs"), py::arg("num"), py::arg("lambdas"))
    .def("BulkModulusFromDensityTemperature", &BulkModulusFromDensityTemperature<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("bmods"), py::arg("num"), py::arg("lambdas"))
    .def("BulkModulusFromDensityInternalEnergy", &BulkModulusFromDensityInternalEnergy<T>, py::arg("rhos"), py::arg("sies"), py::arg("bmods"), py::arg("num"), py::arg("lambdas"))
    .def("GruneisenParamFromDensityTemperature", &GruneisenParamFromDensityTemperature<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("gm1s"), py::arg("num"), py::arg("lambdas"))
    .def("GruneisenParamFromDensityInternalEnergy", &GruneisenParamFromDensityInternalEnergy<T>, py::arg("rhos"), py::arg("sies"), py::arg("gm1s"), py::arg("num"), py::arg("lambdas"))

    .def("nlambda", &T::nlambda)
    .def_static("PreferredInput", &T::PreferredInput)
    .def("PrintParams", &T::PrintParams)
    // TODO .def("DensityEnergyFromPressureTemperature") reference semantics
    .def("Finalize", &T::Finalize)
    .def_static("EosType", &T::EosType);
}

PYBIND11_MODULE(singularity_eos, m) {
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
