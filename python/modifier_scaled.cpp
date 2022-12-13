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
void scaled_eos_class(py::module_ &m) {
  // define Scaled utility function
  m.def(
      "Scaled",
      [](typename T::BaseType eos, Real scaled) { return T(std::move(eos), scaled); },
      py::arg("eos"), py::arg("scaled"));

  // define scaled class
  eos_class<T>(m, T::EosPyType())
      .def(py::init<typename T::BaseType, Real>(), py::arg("eos"), py::arg("scale"));
}

template <typename... Ts>
void create_scaled(py::module_ &m, tl<Ts...>) {
  // C++14 workaround, since we don't have C++17 fold expressions
  auto l = {(scaled_eos_class<Ts>(m), 0)...};
}

void create_scaled_eos_classes(py::module_ &m) {
  create_scaled(m, singularity::scaled);
  create_scaled(m, singularity::scaled_of_shifted);
}
