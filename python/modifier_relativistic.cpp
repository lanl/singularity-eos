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
#include "module.hpp"

template<typename T>
void relativistic_eos_class(py::module_ & m) {
  // define Relativistic utility function
  m.def("Relativistic", [](typename T::BaseType eos, const Real cl){
    return T(std::move(eos), cl);
  }, py::arg("eos"), py::arg("cl"));

  eos_class<T>(m, T::EosPyType())
    .def(
      py::init<typename T::BaseType, Real>(),
      py::arg("eos"), py::arg("cl")
    );
}

template<typename... Ts>
void create_relativistic(py::module_ &m, tl<Ts...>) {
  // C++14 workaround, since we don't have C++17 fold expressions
  auto l = { (relativistic_eos_class<Ts>(m), 0)... };
}

void create_relativistic_eos_classes(py::module_ & m) {
  create_relativistic(m, singularity::relativistic);
}
