//------------------------------------------------------------------------------
// © 2021-2023. Triad National Security, LLC. All rights reserved.  This
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
#include <cmath>

namespace py = pybind11;
using namespace singularity;

// Helper function to convert lambda numpy array to double* buffer
// With std::optional we would add support for a default value of lambda=None
template<typename T, PORTABLE_FUNCTION Real(T::*Func)(const Real, const Real, Real*&&) const>
Real two_params(const T& self, const Real a, const Real b, py::array_t<Real> lambda) {
  return (self.*Func)(a, b, lambda.mutable_data());
}

template<typename T, PORTABLE_FUNCTION Real(T::*Func)(const Real, const Real, Real*&&) const>
Real two_params_no_lambda(const T& self, const Real a, const Real b) {
  return (self.*Func)(a, b, static_cast<Real*>(nullptr));
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
void func(const T & self, py::array_t<Real> a, py::array_t<Real> b,               \
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
}                                                                                 \
                                                                                  \
template<typename T>                                                              \
void func####WithScratch(const T & self, py::array_t<Real> a, py::array_t<Real> b,\
          py::array_t<Real> out, py::array_t<Real> scratch, const int num, py::array_t<Real> lambdas){       \
  py::buffer_info lambdas_info = lambdas.request();                               \
  if (lambdas_info.ndim != 2)                                                     \
      throw std::runtime_error("lambdas dimension must be 2!");                   \
                                                                                  \
  if(lambdas_info.shape[1] > 0) {                                                 \
    self.func(a.data(), b.data(), out.mutable_data(), scratch.mutable_data(), num, LambdaHelper(lambdas)); \
  } else {                                                                        \
    self.func(a.data(), b.data(), out.mutable_data(), scratch.mutable_data(), num, NoLambdaHelper());      \
  }                                                                               \
}                                                                                 \
                                                                                  \
template<typename T>                                                              \
void func####NoLambda(const T & self, py::array_t<Real> a, py::array_t<Real> b,   \
          py::array_t<Real> out, const int num){                                  \
  self.func(a.data(), b.data(), out.mutable_data(), num, NoLambdaHelper());       \
}                                                                                 \
                                                                                  \
template<typename T>                                                              \
void func####NoLambdaWithScratch(const T & self, py::array_t<Real> a, py::array_t<Real> b,  \
          py::array_t<Real> out, py::array_t<Real> scratch, const int num){                 \
  self.func(a.data(), b.data(), out.mutable_data(), scratch.mutable_data(), num, NoLambdaHelper()); \
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
  Real dpde;
  Real dvdt;
  Real dpdr;
  Real dtdr;
  Real dtde;

  EOSState() :
    density(std::numeric_limits<Real>::quiet_NaN()),
    specific_internal_energy(std::numeric_limits<Real>::quiet_NaN()),
    pressure(std::numeric_limits<Real>::quiet_NaN()),
    temperature(std::numeric_limits<Real>::quiet_NaN()),
    specific_heat(std::numeric_limits<Real>::quiet_NaN()),
    bulk_modulus(std::numeric_limits<Real>::quiet_NaN()),
    dpde(std::numeric_limits<Real>::quiet_NaN()),
    dvdt(std::numeric_limits<Real>::quiet_NaN()),
    dpdr(std::numeric_limits<Real>::quiet_NaN()),
    dtdr(std::numeric_limits<Real>::quiet_NaN()),
    dtde(std::numeric_limits<Real>::quiet_NaN()) {
  }

  std::string to_string() const {
    std::stringstream ss;
    if(!std::isnan(density)) ss << "density: " << density << std::endl;
    if(!std::isnan(specific_internal_energy)) ss << "specific_internal_energy: " << specific_internal_energy << std::endl;
    if(!std::isnan(pressure)) ss << "pressure: " << pressure << std::endl;
    if(!std::isnan(temperature)) ss << "temperature: " << temperature << std::endl;
    if(!std::isnan(specific_heat)) ss << "specific_heat: " << specific_heat << std::endl;
    if(!std::isnan(bulk_modulus)) ss << "bulk_modulus: " << bulk_modulus << std::endl;
    if(!std::isnan(dpde)) ss << "dpde:" << dpde << std::endl;
    if(!std::isnan(dvdt)) ss << "dvdt:" << dvdt << std::endl;
    if(!std::isnan(dpdr)) ss << "dpdr:" << dpdr << std::endl;
    if(!std::isnan(dtdr)) ss << "dtdr:" << dtdr << std::endl;
    if(!std::isnan(dtde)) ss << "dtde:" << dtde << std::endl;
    return ss.str();
  }
};

template<typename T, bool use_scratch = false>
struct VectorFunctions {
  static void add(py::class_<T> & cls) {
    // Generic functions provided by the base class. These contain e.g. the vector
    // overloads that use the scalar versions declared here
    cls.def("TemperatureFromDensityInternalEnergy", &TemperatureFromDensityInternalEnergy<T>, py::arg("rhos"), py::arg("sies"), py::arg("temperatures"), py::arg("num"), py::arg("lmbdas"))
    .def("InternalEnergyFromDensityTemperature", &InternalEnergyFromDensityTemperature<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("sies"), py::arg("num"), py::arg("lmbdas"))
    .def("PressureFromDensityTemperature", &PressureFromDensityTemperature<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("pressures"), py::arg("num"), py::arg("lmbdas"))
    .def("PressureFromDensityInternalEnergy", &PressureFromDensityInternalEnergy<T>, py::arg("rhos"), py::arg("sies"), py::arg("pressures"), py::arg("num"), py::arg("lmbdas"))
    .def("SpecificHeatFromDensityTemperature", &SpecificHeatFromDensityTemperature<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("cvs"), py::arg("num"), py::arg("lmbdas"))
    .def("SpecificHeatFromDensityInternalEnergy", &SpecificHeatFromDensityInternalEnergy<T>, py::arg("rhos"), py::arg("sies"), py::arg("cvs"), py::arg("num"), py::arg("lmbdas"))
    .def("BulkModulusFromDensityTemperature", &BulkModulusFromDensityTemperature<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("bmods"), py::arg("num"), py::arg("lmbdas"))
    .def("BulkModulusFromDensityInternalEnergy", &BulkModulusFromDensityInternalEnergy<T>, py::arg("rhos"), py::arg("sies"), py::arg("bmods"), py::arg("num"), py::arg("lmbdas"))
    .def("GruneisenParamFromDensityTemperature", &GruneisenParamFromDensityTemperature<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("gm1s"), py::arg("num"), py::arg("lmbdas"))
    .def("GruneisenParamFromDensityInternalEnergy", &GruneisenParamFromDensityInternalEnergy<T>, py::arg("rhos"), py::arg("sies"), py::arg("gm1s"), py::arg("num"), py::arg("lmbdas"))

    .def("TemperatureFromDensityInternalEnergy", &TemperatureFromDensityInternalEnergyNoLambda<T>, py::arg("rhos"), py::arg("sies"), py::arg("temperatures"), py::arg("num"))
    .def("InternalEnergyFromDensityTemperature", &InternalEnergyFromDensityTemperatureNoLambda<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("sies"), py::arg("num"))
    .def("PressureFromDensityTemperature", &PressureFromDensityTemperatureNoLambda<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("pressures"), py::arg("num"))
    .def("PressureFromDensityInternalEnergy", &PressureFromDensityInternalEnergyNoLambda<T>, py::arg("rhos"), py::arg("sies"), py::arg("pressures"), py::arg("num"))
    .def("SpecificHeatFromDensityTemperature", &SpecificHeatFromDensityTemperatureNoLambda<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("cvs"), py::arg("num"))
    .def("SpecificHeatFromDensityInternalEnergy", &SpecificHeatFromDensityInternalEnergyNoLambda<T>, py::arg("rhos"), py::arg("sies"), py::arg("cvs"), py::arg("num"))
    .def("BulkModulusFromDensityTemperature", &BulkModulusFromDensityTemperatureNoLambda<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("bmods"), py::arg("num"))
    .def("BulkModulusFromDensityInternalEnergy", &BulkModulusFromDensityInternalEnergyNoLambda<T>, py::arg("rhos"), py::arg("sies"), py::arg("bmods"), py::arg("num"))
    .def("GruneisenParamFromDensityTemperature", &GruneisenParamFromDensityTemperatureNoLambda<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("gm1s"), py::arg("num"))
    .def("GruneisenParamFromDensityInternalEnergy", &GruneisenParamFromDensityInternalEnergyNoLambda<T>, py::arg("rhos"), py::arg("sies"), py::arg("gm1s"), py::arg("num"))


    .def("FillEos", [](const T & self, py::array_t<Real> rhos,
                       py::array_t<Real> temperatures, py::array_t<Real> sies,
                       py::array_t<Real> pressures, py::array_t<Real> cvs, py::array_t<Real> bmods,
                       const int num, const unsigned long output, py::array_t<Real> lambdas) {
      py::buffer_info lambdas_info = lambdas.request();
      if (lambdas_info.ndim != 2)
        throw std::runtime_error("lambdas dimension must be 2!");

      if(lambdas_info.shape[1] > 0) {
        self.FillEos(rhos.mutable_data(), temperatures.mutable_data(),
                     sies.mutable_data(), pressures.mutable_data(), cvs.mutable_data(),
                     bmods.mutable_data(), num, output, LambdaHelper(lambdas));
      } else {
        self.FillEos(rhos.mutable_data(), temperatures.mutable_data(),
                     sies.mutable_data(), pressures.mutable_data(), cvs.mutable_data(),
                     bmods.mutable_data(), num, output, NoLambdaHelper());
      }
    }, py::arg("rhos"), py::arg("temperatures"), py::arg("sies"),
       py::arg("pressures"), py::arg("cvs"), py::arg("bmods"), py::arg("num"),
       py::arg("output"), py::arg("lmbdas")
    )
    .def("FillEos", [](const T & self, py::array_t<Real> rhos,
                       py::array_t<Real> temperatures, py::array_t<Real> sies, py::array_t<Real>
                       pressures, py::array_t<Real> cvs, py::array_t<Real> bmods, const int num,
                       const unsigned long output) {
      self.FillEos(rhos.mutable_data(), temperatures.mutable_data(),
                   sies.mutable_data(), pressures.mutable_data(), cvs.mutable_data(),
                   bmods.mutable_data(), num, output, NoLambdaHelper());
    }, py::arg("rhos"), py::arg("temperatures"), py::arg("sies"), py::arg("pressures"), py::arg("cvs"), py::arg("bmods"), py::arg("num"), py::arg("output"));
  }
};

template<typename T>
struct VectorFunctions<T,true> {
  static void add(py::class_<T> & cls) {
    cls.def("TemperatureFromDensityInternalEnergy", &TemperatureFromDensityInternalEnergyWithScratch<T>, py::arg("rhos"), py::arg("sies"), py::arg("temperatures"), py::arg("scratch"), py::arg("num"), py::arg("lmbdas"))
    .def("InternalEnergyFromDensityTemperature", &InternalEnergyFromDensityTemperatureWithScratch<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("sies"), py::arg("scratch"), py::arg("num"), py::arg("lmbdas"))
    .def("PressureFromDensityTemperature", &PressureFromDensityTemperatureWithScratch<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("pressures"), py::arg("scratch"), py::arg("num"), py::arg("lmbdas"))
    .def("PressureFromDensityInternalEnergy", &PressureFromDensityInternalEnergyWithScratch<T>, py::arg("rhos"), py::arg("sies"), py::arg("pressures"), py::arg("scratch"), py::arg("num"), py::arg("lmbdas"))
    .def("SpecificHeatFromDensityTemperature", &SpecificHeatFromDensityTemperatureWithScratch<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("cvs"), py::arg("scratch"), py::arg("num"), py::arg("lmbdas"))
    .def("SpecificHeatFromDensityInternalEnergy", &SpecificHeatFromDensityInternalEnergyWithScratch<T>, py::arg("rhos"), py::arg("sies"), py::arg("cvs"), py::arg("scratch"), py::arg("num"), py::arg("lmbdas"))
    .def("BulkModulusFromDensityTemperature", &BulkModulusFromDensityTemperatureWithScratch<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("bmods"), py::arg("scratch"), py::arg("num"), py::arg("lmbdas"))
    .def("BulkModulusFromDensityInternalEnergy", &BulkModulusFromDensityInternalEnergyWithScratch<T>, py::arg("rhos"), py::arg("sies"), py::arg("bmods"), py::arg("scratch"), py::arg("num"), py::arg("lmbdas"))
    .def("GruneisenParamFromDensityTemperature", &GruneisenParamFromDensityTemperatureWithScratch<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("gm1s"), py::arg("scratch"), py::arg("num"), py::arg("lmbdas"))
    .def("GruneisenParamFromDensityInternalEnergy", &GruneisenParamFromDensityInternalEnergyWithScratch<T>, py::arg("rhos"), py::arg("sies"), py::arg("gm1s"), py::arg("scratch"), py::arg("num"), py::arg("lmbdas"))

    .def("TemperatureFromDensityInternalEnergy", &TemperatureFromDensityInternalEnergyNoLambdaWithScratch<T>, py::arg("rhos"), py::arg("sies"), py::arg("temperatures"), py::arg("scratch"), py::arg("num"))
    .def("InternalEnergyFromDensityTemperature", &InternalEnergyFromDensityTemperatureNoLambdaWithScratch<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("sies"), py::arg("scratch"), py::arg("num"))
    .def("PressureFromDensityTemperature", &PressureFromDensityTemperatureNoLambdaWithScratch<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("pressures"), py::arg("scratch"), py::arg("num"))
    .def("PressureFromDensityInternalEnergy", &PressureFromDensityInternalEnergyNoLambdaWithScratch<T>, py::arg("rhos"), py::arg("sies"), py::arg("pressures"), py::arg("scratch"), py::arg("num"))
    .def("SpecificHeatFromDensityTemperature", &SpecificHeatFromDensityTemperatureNoLambdaWithScratch<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("cvs"), py::arg("scratch"), py::arg("num"))
    .def("SpecificHeatFromDensityInternalEnergy", &SpecificHeatFromDensityInternalEnergyNoLambdaWithScratch<T>, py::arg("rhos"), py::arg("sies"), py::arg("cvs"), py::arg("scratch"), py::arg("num"))
    .def("BulkModulusFromDensityTemperature", &BulkModulusFromDensityTemperatureNoLambdaWithScratch<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("bmods"), py::arg("scratch"), py::arg("num"))
    .def("BulkModulusFromDensityInternalEnergy", &BulkModulusFromDensityInternalEnergyNoLambdaWithScratch<T>, py::arg("rhos"), py::arg("sies"), py::arg("bmods"), py::arg("scratch"), py::arg("num"))
    .def("GruneisenParamFromDensityTemperature", &GruneisenParamFromDensityTemperatureNoLambdaWithScratch<T>, py::arg("rhos"), py::arg("temperatures"), py::arg("gm1s"), py::arg("scratch"), py::arg("num"))
    .def("GruneisenParamFromDensityInternalEnergy", &GruneisenParamFromDensityInternalEnergyNoLambdaWithScratch<T>, py::arg("rhos"), py::arg("sies"), py::arg("gm1s"), py::arg("scratch"), py::arg("num"))

    .def("FillEos", [](const T & self, py::array_t<Real> rhos,
                       py::array_t<Real> temperatures, py::array_t<Real> sies,
                       py::array_t<Real> pressures, py::array_t<Real> cvs, py::array_t<Real> bmods,
                       py::array_t<Real> scratch,
                       const int num, const unsigned long output, py::array_t<Real> lambdas) {
      py::buffer_info lambdas_info = lambdas.request();
      if (lambdas_info.ndim != 2)
        throw std::runtime_error("lambdas dimension must be 2!");

      if(lambdas_info.shape[1] > 0) {
        self.FillEos(rhos.mutable_data(), temperatures.mutable_data(),
                     sies.mutable_data(), pressures.mutable_data(), cvs.mutable_data(),
                     bmods.mutable_data(), scratch.mutable_data(), num, output, LambdaHelper(lambdas));
      } else {
        self.FillEos(rhos.mutable_data(), temperatures.mutable_data(),
                     sies.mutable_data(), pressures.mutable_data(), cvs.mutable_data(),
                     bmods.mutable_data(), scratch.mutable_data(), num, output, NoLambdaHelper());
      }
    }, py::arg("rhos"), py::arg("temperatures"), py::arg("sies"),
       py::arg("pressures"), py::arg("cvs"), py::arg("bmods"), py::arg("scratch"), py::arg("num"),
       py::arg("output"), py::arg("lmbdas")
    )
    .def("FillEos", [](const T & self, py::array_t<Real> rhos,
                       py::array_t<Real> temperatures, py::array_t<Real> sies, py::array_t<Real>
                       pressures, py::array_t<Real> cvs, py::array_t<Real> bmods,
                       py::array_t<Real> scratch,
                       const int num, const unsigned long output) {
      self.FillEos(rhos.mutable_data(), temperatures.mutable_data(),
                   sies.mutable_data(), pressures.mutable_data(), cvs.mutable_data(),
                   bmods.mutable_data(), scratch.mutable_data(), num, output, NoLambdaHelper());
    }, py::arg("rhos"), py::arg("temperatures"), py::arg("sies"), py::arg("pressures"), py::arg("cvs"), py::arg("bmods"), py::arg("scratch"), py::arg("num"), py::arg("output"));
  }
};

template<typename T, bool use_scratch = false>
py::class_<T> eos_class(py::module_ & m, std::string name) {
  py::class_<T> cls(m, name.c_str());
  cls.def("TemperatureFromDensityInternalEnergy", &two_params<T, &T::TemperatureFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"), py::arg("lmbda"))
    .def("InternalEnergyFromDensityTemperature", &two_params<T, &T::InternalEnergyFromDensityTemperature>, py::arg("rho"), py::arg("temperature"), py::arg("lmbda"))
    .def("PressureFromDensityTemperature", &two_params<T, &T::PressureFromDensityTemperature>, py::arg("rho"), py::arg("temperature"), py::arg("lmbda"))
    .def("PressureFromDensityInternalEnergy", &two_params<T, &T::PressureFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"), py::arg("lmbda"))
    .def("SpecificHeatFromDensityTemperature", &two_params<T, &T::SpecificHeatFromDensityTemperature>, py::arg("rho"), py::arg("temperature"), py::arg("lmbda"))
    .def("SpecificHeatFromDensityInternalEnergy", &two_params<T, &T::SpecificHeatFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"), py::arg("lmbda"))
    .def("BulkModulusFromDensityTemperature", &two_params<T, &T::BulkModulusFromDensityTemperature>, py::arg("rho"), py::arg("temperature"), py::arg("lmbda"))
    .def("BulkModulusFromDensityInternalEnergy", &two_params<T, &T::BulkModulusFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"), py::arg("lmbda"))
    .def("GruneisenParamFromDensityTemperature", &two_params<T, &T::GruneisenParamFromDensityTemperature>, py::arg("rho"), py::arg("temperature"), py::arg("lmbda"))
    .def("GruneisenParamFromDensityInternalEnergy", &two_params<T, &T::GruneisenParamFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"), py::arg("lmbda"))

    .def("TemperatureFromDensityInternalEnergy", &two_params_no_lambda<T, &T::TemperatureFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"))
    .def("InternalEnergyFromDensityTemperature", &two_params_no_lambda<T, &T::InternalEnergyFromDensityTemperature>, py::arg("rho"), py::arg("temperature"))
    .def("PressureFromDensityTemperature", &two_params_no_lambda<T, &T::PressureFromDensityTemperature>, py::arg("rho"), py::arg("temperature"))
    .def("PressureFromDensityInternalEnergy", &two_params_no_lambda<T, &T::PressureFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"))
    .def("SpecificHeatFromDensityTemperature", &two_params_no_lambda<T, &T::SpecificHeatFromDensityTemperature>, py::arg("rho"), py::arg("temperature"))
    .def("SpecificHeatFromDensityInternalEnergy", &two_params_no_lambda<T, &T::SpecificHeatFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"))
    .def("BulkModulusFromDensityTemperature", &two_params_no_lambda<T, &T::BulkModulusFromDensityTemperature>, py::arg("rho"), py::arg("temperature"))
    .def("BulkModulusFromDensityInternalEnergy", &two_params_no_lambda<T, &T::BulkModulusFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"))
    .def("GruneisenParamFromDensityTemperature", &two_params_no_lambda<T, &T::GruneisenParamFromDensityTemperature>, py::arg("rho"), py::arg("temperature"))
    .def("GruneisenParamFromDensityInternalEnergy", &two_params_no_lambda<T, &T::GruneisenParamFromDensityInternalEnergy>, py::arg("rho"), py::arg("sie"))

    .def("FillEos", [](const T & self, const py::kwargs& kwargs) {
      unsigned long output = thermalqs::none;
      EOSState s;
      std::map<std::string, std::pair<unsigned long,Real*>> param_mapping {
        {"rho", {thermalqs::density, &s.density}},
        {"sie", {thermalqs::specific_internal_energy, &s.specific_internal_energy}},
        {"press", {thermalqs::pressure, &s.pressure}},
        {"temp", {thermalqs::temperature, &s.temperature}},
        {"cv", {thermalqs::specific_heat, &s.specific_heat}},
        {"bmod", {thermalqs::bulk_modulus, &s.bulk_modulus}},
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
        self.FillEos(s.density, s.temperature, s.specific_internal_energy, s.pressure, s.specific_heat, s.bulk_modulus, output, lambda.mutable_data());
      } else {
        self.FillEos(s.density, s.temperature, s.specific_internal_energy, s.pressure, s.specific_heat, s.bulk_modulus, output, nullptr);
      }
      return s;
    })
    .def("ValuesAtReferenceState", [](const T & self, py::array_t<Real> lambda){
      EOSState s;
      self.ValuesAtReferenceState(s.density, s.temperature, s.specific_internal_energy, s.pressure, s.specific_heat, s.bulk_modulus, s.dpde, s.dvdt, lambda.mutable_data());
      return s;
    })
    .def("ValuesAtReferenceState", [](const T & self){
      EOSState s;
      self.ValuesAtReferenceState(s.density, s.temperature, s.specific_internal_energy, s.pressure, s.specific_heat, s.bulk_modulus, s.dpde, s.dvdt, nullptr);
      return s;
    })


    .def("MinimumDensity", &T::MinimumDensity)
    .def("MinimumTemperature", &T::MinimumTemperature)
    .def_property_readonly("nlambda", &T::nlambda)
    .def_property_readonly_static("PreferredInput", [](py::object) { return T::PreferredInput(); })
    .def("PrintParams", &T::PrintParams)
    .def("DensityEnergyFromPressureTemperature", [](const T & self, const Real press, const Real temp, py::array_t<Real> lambda) {
      Real rho, sie;
      self.DensityEnergyFromPressureTemperature(press, temp, lambda.mutable_data(), rho, sie);
      return std::pair<Real, Real>(rho, sie);
    }, py::arg("press"), py::arg("temp"), py::arg("lmbda"))
    .def("DensityEnergyFromPressureTemperature", [](const T & self, const Real press, const Real temp) {
      Real rho, sie;
      self.DensityEnergyFromPressureTemperature(press, temp, nullptr, rho, sie);
      return std::pair<Real, Real>(rho, sie);
    }, py::arg("press"), py::arg("temp"))
    .def("Finalize", &T::Finalize)
    .def_property_readonly_static("EosType", [](py::object) { return T::EosType(); })
    .def_static("scratch_size", &T::scratch_size, py::arg("method_name"), py::arg("nelements"))
    .def_static("max_scratch_size", &T::max_scratch_size, py::arg("nelements"));

    VectorFunctions<T, use_scratch>::add(cls);

    return cls;
}

// functions to move modifier instantiations into their own compilation unit
// allows parallel compilation and reduces compiler memory footprint
void create_shifted_eos_classes(py::module_ & m);
void create_scaled_eos_classes(py::module_ & m);
void create_bilinear_ramp_eos_classes(py::module_ &m);
void create_relativistic_eos_classes(py::module_ &m);
void create_unit_system_eos_classes(py::module_ &m);
