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

template <typename T>
void bilinear_ramp_helper(py::module_ &m) {
  m.def(
      "BilinearRamp",
      [](T eos, const Real r0, const Real a, const Real b, const Real c) {
        return BilinearRampEOS<T>(std::move(eos), r0, a, b, c);
      },
      py::arg("eos"), py::arg("r0"), py::arg("a"), py::arg("b"), py::arg("c"));
}

template <typename T>
void bilinear_ramp_class_helper(py::module_ &m, std::string name) {
  eos_class<BilinearRampEOS<T>>(m, std::string("BilinearRamp") + name)
      .def(py::init<T, Real, Real, Real, Real>(), py::arg("eos"), py::arg("r0"),
           py::arg("a"), py::arg("b"), py::arg("c"));
}

template <typename T>
void bilinear_ramp_eos_class(py::module_ &m, const char *name) {
  // define BilinearRamp utility function
  bilinear_ramp_helper<T>(m);
  bilinear_ramp_helper<ShiftedEOS<T>>(m);
  bilinear_ramp_helper<ScaledEOS<T>>(m);
  bilinear_ramp_helper<ScaledEOS<ShiftedEOS<T>>>(m);

  bilinear_ramp_class_helper<T>(m, name);
  bilinear_ramp_class_helper<ShiftedEOS<T>>(m, std::string("Shifted") + name);
  bilinear_ramp_class_helper<ScaledEOS<T>>(m, std::string("Scaled") + name);
  bilinear_ramp_class_helper<ScaledEOS<ShiftedEOS<T>>>(m, std::string("ScaledShifted") + name);
}

void create_bilinear_ramp_eos_classes(py::module_ &m) {
  bilinear_ramp_eos_class<IdealGas>(m, "IdealGas");

#ifdef SPINER_USE_HDF
  bilinear_ramp_eos_class<SpinerEOSDependsRhoT>(m, "SpinerEOSDependsRhoT");
  bilinear_ramp_eos_class<SpinerEOSDependsRhoSie>(m, "SpinerEOSDependsRhoSie");
  bilinear_ramp_eos_class<StellarCollapse>(m, "StellarCollapse");
#endif
}
