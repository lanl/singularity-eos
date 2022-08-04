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
py::class_<T> shifted_eos_class(py::module_ &m, const char *name) {
  // define Shifted utility function
  m.def(
      "Shifted", [](T eos, Real shift) { return ShiftedEOS<T>(std::move(eos), shift); },
      py::arg("eos"), py::arg("shift"));

  // define shifted class
  return eos_class<ShiftedEOS<T>>(m, std::string("Shifted") + name)
      .def(py::init<T, Real>(), py::arg("eos"), py::arg("shift"));
}

void create_shifted_eos_classes(py::module_ &m) {
  shifted_eos_class<IdealGas>(m, "IdealGas");
  shifted_eos_class<Gruneisen>(m, "Gruneisen");
  shifted_eos_class<JWL>(m, "JWL");
  shifted_eos_class<DavisReactants>(m, "DavisReactants");
  shifted_eos_class<DavisProducts>(m, "DavisProducts");

#ifdef SPINER_USE_HDF
  shifted_eos_class<SpinerEOSDependsRhoT>(m, "SpinerEOSDependsRhoT");
  shifted_eos_class<SpinerEOSDependsRhoSie>(m, "SpinerEOSDependsRhoSie");
  shifted_eos_class<StellarCollapse>(m, "StellarCollapse");
#endif

#ifdef SINGULARITY_USE_EOSPAC
  shifted_eos_class<EOSPAC>(m, "EOSPAC");
#endif
}
