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
void relativistic_helper(py::module_ & m) {
  m.def("Relativistic", [](T eos, const Real cl){
    return RelativisticEOS<T>(std::move(eos), cl);
  }, py::arg("eos"), py::arg("cl"));
}

template<typename T>
void relativistic_class_helper(py::module_ & m, std::string name) {
  eos_class<RelativisticEOS<T>>(m, std::string("Relativistic") + name)
    .def(
      py::init<T, Real>(),
      py::arg("eos"), py::arg("cl")
    );
}

template<typename T>
void relativistic_eos_class(py::module_ & m, const char * name) {
  // define Relativistic utility function
  relativistic_helper<T>(m);
  relativistic_helper<ShiftedEOS<T>>(m);
  relativistic_helper<ScaledEOS<T>>(m);
  relativistic_helper<ScaledEOS<ShiftedEOS<T>>>(m);

  relativistic_class_helper<T>(m, name);
  relativistic_class_helper<ShiftedEOS<T>>(m, std::string("Shifted") + name);
  relativistic_class_helper<ScaledEOS<T>>(m, std::string("Scaled") + name);
  relativistic_class_helper<ScaledEOS<ShiftedEOS<T>>>(m, std::string("ScaledShifted") + name);

  relativistic_helper<BilinearRampEOS<T>>(m);
  relativistic_helper<BilinearRampEOS<ShiftedEOS<T>>>(m);
  relativistic_helper<BilinearRampEOS<ScaledEOS<T>>>(m);
  relativistic_helper<BilinearRampEOS<ScaledEOS<ShiftedEOS<T>>>>(m);

  relativistic_class_helper<BilinearRampEOS<T>>(m, std::string("Bilinear") + name);
  relativistic_class_helper<BilinearRampEOS<ShiftedEOS<T>>>(m, std::string("BilinearShifted") + name);
  relativistic_class_helper<BilinearRampEOS<ScaledEOS<T>>>(m, std::string("BilinearScaled") + name);
  relativistic_class_helper<BilinearRampEOS<ScaledEOS<ShiftedEOS<T>>>>(m, std::string("BilinearScaledShifted") + name);
}


void create_relativistic_eos_classes(py::module_ & m) {
  relativistic_eos_class<IdealGas>(m, "IdealGas");

#ifdef SPINER_USE_HDF
  relativistic_eos_class<SpinerEOSDependsRhoT>(m, "SpinerEOSDependsRhoT");
  relativistic_eos_class<SpinerEOSDependsRhoSie>(m, "SpinerEOSDependsRhoSie");
  relativistic_eos_class<StellarCollapse>(m, "StellarCollapse");
#endif
}
