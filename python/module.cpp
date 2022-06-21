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
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <singularity-eos/eos/eos.hpp>

namespace py = pybind11;
using namespace singularity;

// Helper function to convert lambda numpy array to double* buffer
// With std::optional we would add support for a default value of lambda=None
template<typename T, PORTABLE_FUNCTION Real(T::*Func)(const Real, const Real, Real*) const>
Real two_params(const T& obj, const Real a, const Real b, py::array_t<double> lambda) {
  py::buffer_info lambda_info = lambda.request();
  auto lambda_ptr = static_cast<double *>(lambda_info.ptr);
  return (obj.*Func)(a, b, lambda_ptr);
}

PYBIND11_MODULE(singularity_eos, m) {
  py::class_<IdealGas>(m, "IdealGas")
    .def(py::init<Real, Real>())
    .def("TemperatureFromDensityInternalEnergy", &two_params<IdealGas, &IdealGas::TemperatureFromDensityInternalEnergy>)
    .def("InternalEnergyFromDensityTemperature", &two_params<IdealGas, &IdealGas::InternalEnergyFromDensityTemperature>)
    .def("PressureFromDensityTemperature", &two_params<IdealGas, &IdealGas::PressureFromDensityTemperature>)
    .def("PressureFromDensityInternalEnergy", &two_params<IdealGas, &IdealGas::PressureFromDensityInternalEnergy>)
    .def("SpecificHeatFromDensityTemperature", &two_params<IdealGas, &IdealGas::SpecificHeatFromDensityTemperature>)
    .def("SpecificHeatFromDensityInternalEnergy", &two_params<IdealGas, &IdealGas::SpecificHeatFromDensityInternalEnergy>)
    .def("BulkModulusFromDensityTemperature", &two_params<IdealGas, &IdealGas::BulkModulusFromDensityTemperature>)
    .def("BulkModulusFromDensityInternalEnergy", &two_params<IdealGas, &IdealGas::BulkModulusFromDensityInternalEnergy>)
    .def("GruneisenParamFromDensityTemperature", &two_params<IdealGas, &IdealGas::GruneisenParamFromDensityTemperature>)
    .def("GruneisenParamFromDensityInternalEnergy", &two_params<IdealGas, &IdealGas::GruneisenParamFromDensityInternalEnergy>)

    // TODO .def("GetOnDevice")
    // TODO .def("FillEos")  what are the call semantics, since it uses references? return dict?
    // TODO .def("ValuesAtReferenceState")

    // Generic functions provided by the base class. These contain e.g. the vector
    // overloads that use the scalar versions declared here
    // TODO  SG_ADD_BASE_CLASS_USINGS(IdealGas)
    .def("nlambda", &IdealGas::nlambda)
    .def_static("PreferredInput", &IdealGas::PreferredInput)
    .def("PrintParams", &IdealGas::PrintParams)
    // TODO .def("DensityEnergyFromPressureTemperature") reference semantics
    .def("Finalize", &IdealGas::Finalize)
    .def_static("EosType", &IdealGas::EosType);

  m.doc() = "Singularity EOS Python Bindings";
}
