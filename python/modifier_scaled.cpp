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
#include "module.hpp"

template <typename T>
void scaled_eos_class(py::module_ &m, const char *name) {
  // define Scaled utility function
  m.def(
      "Scaled", [](T eos, Real scaled) { return ScaledEOS<T>(std::move(eos), scaled); },
      py::arg("eos"), py::arg("scaled"));
  m.def(
      "Scaled",
      [](ShiftedEOS<T> eos, Real scaled) {
        return ScaledEOS<ShiftedEOS<T>>(std::move(eos), scaled);
      },
      py::arg("eos"), py::arg("scaled"));

  // each scaled can also be shifted
  eos_class<ScaledEOS<ShiftedEOS<T>>>(m, std::string("ScaledShifted") + name)
      .def(py::init<ShiftedEOS<T>, Real>(), py::arg("eos"), py::arg("scale"));

  // define scaled class
  eos_class<ScaledEOS<T>>(m, std::string("Scaled") + name)
      .def(py::init<T, Real>(), py::arg("eos"), py::arg("scale"));
}

void create_scaled_eos_classes(py::module_ &m) {
  scaled_eos_class<IdealGas>(m, "IdealGas");
  scaled_eos_class<Gruneisen>(m, "Gruneisen");
  scaled_eos_class<JWL>(m, "JWL");
  scaled_eos_class<DavisReactants>(m, "DavisReactants");
  scaled_eos_class<DavisProducts>(m, "DavisProducts");

#ifdef SPINER_USE_HDF
  scaled_eos_class<SpinerEOSDependsRhoT>(m, "SpinerEOSDependsRhoT");
  scaled_eos_class<SpinerEOSDependsRhoSie>(m, "SpinerEOSDependsRhoSie");
  scaled_eos_class<StellarCollapse>(m, "StellarCollapse");
#endif

#ifdef SINGULARITY_USE_EOSPAC
  scaled_eos_class<EOSPAC>(m, "EOSPAC");
#endif
}
