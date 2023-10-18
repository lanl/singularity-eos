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
constexpr bool is_modifiable(const Variant<Ts...> &var) {
  return IsModifiable<Mod, U, Ts...>::value;
}
template <template <class...> typename Mod, typename U, typename... Ts>
constexpr bool is_modifiable(const U &u, const Variant<Ts...> &var) {
  return IsModifiable<Mod, U, Ts...>::value;
}

// Recursive functions needed for the Modify function JMM: Note this
// machinery would be a LOT easier in C++17, with constexpr if.  To
// mimic constexpr if, I use tag dispatch here. I need two separate
// functions to make it work, however.
namespace detail {

// Prototypes
template <template <class...> typename Mod, typename U, typename... Rest, typename... Ts,
          typename... Args>
Variant<Ts...> ModifyHelper(std::true_type, const Variant<Ts...> &var,
                            const type_list<U, Rest...> &tl, Args &&...args);
template <template <class...> typename Mod, typename U, typename... Rest, typename... Ts,
          typename... Args>
Variant<Ts...> ModifyHelper(std::false_type, const Variant<Ts...> &var,
                            const type_list<U, Rest...> &tl, Args &&...args);

// Base case
template <template <class...> typename Mod, typename... Ts, typename... Args>
constexpr auto ModifyDispatcher(const Variant<Ts...> &var, const type_list<> &tl,
                                Args &&...args) {
  return var;
}

// The ModifyDispatcher function specifies which tag in tag dispatch
// to use. It's main purpose is to check whether or not U and Mod<U>
// are in the type list for the variant. If yes, it dispatches to the
// ModifyHelper<Mod>(yes,...) method, otherwise, the
// ModifyHelper<Mod>(no,...)  method.

template <template <class...> typename Mod, typename U, typename... Rest, typename... Ts,
          typename... Args>
constexpr auto ModifyDispatcher(const Variant<Ts...> &var,
                                const type_list<U, Rest...> &tl, Args &&...args) {
  constexpr bool type_in_variant = IsModifiable<Mod, U, Ts...>::value;
  return ModifyHelper<Mod>(bool_constant<type_in_variant>(), var, tl,
                           std::forward<Args>(args)...);
}

// The ModifyHelper function recursively walks through all types in the
// typelist tl. For each type in tl, it checks if the type is (a) in
// the variant typelist, and (b) currently held by the variant
// object. If both (a) and (b) are true, it returns a new variant
// containing Mod<u>. Otherwise, it moves to the next type.
template <template <class...> typename Mod, typename U, typename... Rest, typename... Ts,
          typename... Args>
Variant<Ts...> ModifyHelper(std::true_type, const Variant<Ts...> &var,
                            const type_list<U, Rest...> &tl, Args &&...args) {
  bool contains_u = var.template IsType<U>();
  if (contains_u) {
    return Variant<Ts...>(Mod<U>(var.template get<U>(), std::forward<Args>(args)...));
  } else {
    type_list<Rest...> tl_new{};
    return ModifyDispatcher<Mod>(var, tl_new, std::forward<Args>(args)...);
  }
}
template <template <class...> typename Mod, typename U, typename... Rest, typename... Ts,
          typename... Args>
Variant<Ts...> ModifyHelper(std::false_type, const Variant<Ts...> &var,
                            const type_list<U, Rest...> &tl, Args &&...args) {
  type_list<Rest...> tl_new{};
  return ModifyDispatcher<Mod>(var, tl_new, std::forward<Args>(args)...);
}
template <template <class...> typename Mod, typename U, typename... Rest, typename... Ts,
          typename... Args>
Variant<Ts...> ModifyHelper(const Variant<Ts...> &var, const type_list<U, Rest...> &tl,
                            Args &&...args) {
  constexpr bool is_modifiable = IsModifiable<Mod, U, Ts...>::value;
  return ModifyHelper<Mod>(bool_constant<is_modifiable>(), var, tl,
                           std::forward<Args>(args)...);
}
} // namespace detail

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
// Note that this is potentially a pretty expensive operation,
// as it must walk through the entire variant, both at compile, and run-time.
// However, it is also extremely flexible, enabling user-selected variants
// with minimal boiler plate.
template <template <class...> typename Mod, typename... Ts, typename... Args>
Variant<Ts...> Modify(const Variant<Ts...> &var, Args &&...args) {
  type_list<Ts...> tl;
  return detail::ModifyHelper<Mod>(var, tl, std::forward<Args>(args)...);
}

} // namespace EOSBuilder
} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_EOS_BUILDER_HPP_
