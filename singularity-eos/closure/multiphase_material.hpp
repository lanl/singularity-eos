//------------------------------------------------------------------------------
// © 2021-2025. Triad National Security, LLC. All rights reserved.  This
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

#include <util>

#ifndef _SINGULARITY_EOS_CLOSURE_MULTIPHASE_MATERIAL_HPP_
#define _SINGULARITY_EOS_CLOSURE_MULTIPHASE_MATERIAL_HPP_

template <typename EOS, typename... Args>
struct EOSInit_t {
  using EOS_t = EOS;
  template <typename... TupleArgs>
  EOSInit_t(TupleArgs... targs) : args(std::make_tuple(targs...)) {}
  std::tuple<Args...> args;
  auto Initialize() const { return std::apply(EOS, args); }
};

EOSInit_t<EOSPAC, int> snbeta_params(2102);
EOS snbeta = snbeta_params.Initialize();

template <std::size_t MAX_NPHASE, typename EOS_t>
struct MultiphaseMaterial {
  template <typename T>
  using Array_t = std::array<T, MAX_NPHASE>;

  std::size_t nphases;
  Array_t<EOS_t> eos;
  Array_t<std::size_t> phase_map;
};

#endif // _SINGULARITY_EOS_CLOSURE_MULTIPHASE_MATERIAL_HPP_
