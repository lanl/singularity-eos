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

#ifndef _SINGULARITY_EOS_EOS_EOS_BUILDER_HPP_
#define _SINGULARITY_EOS_EOS_EOS_BUILDER_HPP_

#include <utility>

#include <ports-of-call/portability.hpp>

#include <singularity-eos/base/variadic_utils.hpp>
#include <singularity-eos/eos/eos_variant.hpp>

namespace singularity {
namespace EOSBuilder {

// Just import the variadic utils namespace int the EOSBuilder
// namespace.
using namespace variadic_utils;

template <template <class...> typename Mod, typename U, typename... Ts>
using IsModifiable =
    conjunction<contains<U, Ts...>::value, contains<Mod<U>, Ts...>::value>;
template <template <class...> typename Mod, typename U, typename... Ts>
constexpr bool is_modifiable(const type_list<Ts...> &var) {
  return IsModifiable<Mod, U, Ts...>::value;
}
template <template <class...> typename Mod, typename U, typename... Ts>
constexpr bool is_modifiable(const Variant<Ts...> &var) {
  return IsModifiable<Mod, U, Ts...>::value;
}
template <template <class...> typename Mod, typename U, typename... Ts>
constexpr bool is_modifiable(const U &u, const Variant<Ts...> &var) {
  return IsModifiable<Mod, U, Ts...>::value;
}

// Modifies the eos contained in the var object with the modifier Mod,
// assuming such modification is possible, i.e., Mod<T> is in the
// variant, for the underlying type T. If this modification is not
// possible, returns the unmodified EOS. Args are the additional
// arguments to the modifier's constructor.
// Intended usage:
//
// eos = Modify<Modifier>(eos, args);
//
// For example:
//
// EOS eos = IdealGas(gm1, Cv);
// if (shifted) {
//   eos = Modify<ShifteEOS>(eos, shift);
// }
template <template <class> typename Mod, typename... Ts, typename... Args>
Variant<Ts...> Modify(const Variant<Ts...> &var, Args &&...args) {
  return var.template Modify<Mod>(std::forward<Args>(args)...);
}

} // namespace EOSBuilder
} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_EOS_BUILDER_HPP_
