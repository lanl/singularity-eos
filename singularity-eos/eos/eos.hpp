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

#ifndef _SINGULARITY_EOS_EOS_EOS_HPP_
#define _SINGULARITY_EOS_EOS_EOS_HPP_

#include "stdio.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <utility>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/eos/eos_base.hpp>
#include <singularity-eos/eos/eos_variant.hpp>

// Base stuff
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/eos_error.hpp>
#include <singularity-eos/base/variadic_utils.hpp>

// EOS models
#include <singularity-eos/eos/eos_davis.hpp>
#include <singularity-eos/eos/eos_eospac.hpp>
#include <singularity-eos/eos/eos_gruneisen.hpp>
#include <singularity-eos/eos/eos_ideal.hpp>
#include <singularity-eos/eos/eos_jwl.hpp>
#include <singularity-eos/eos/eos_spiner.hpp>
#include <singularity-eos/eos/eos_stellar_collapse.hpp>

// Modifiers
#include <singularity-eos/eos/modifiers/eos_unitsystem.hpp>
#include <singularity-eos/eos/modifiers/ramps_eos.hpp>
#include <singularity-eos/eos/modifiers/relativistic_eos.hpp>
#include <singularity-eos/eos/modifiers/scaled_eos.hpp>
#include <singularity-eos/eos/modifiers/shifted_eos.hpp>

namespace singularity {

// recreate variadic list
template <typename... Ts>
using tl = singularity::detail::type_list<Ts...>;

template <template <typename> typename... Ts>
using al = singularity::detail::adapt_list<Ts...>;

// transform variadic list: applies modifiers to eos's
using singularity::detail::transform_variadic_list;

// all eos's
static constexpr const auto full_eos_list =
    tl<IdealGas, Gruneisen, JWL, DavisReactants, DavisProducts
#ifdef SPINER_USE_HDF
       ,
       SpinerEOSDependsRhoT, SpinerEOSDependsRhoSie, StellarCollapse
#endif // SPINER_USE_HDF
#ifdef SINGULARITY_USE_EOSPAC
       ,
       EOSPAC
#endif // SINGULARITY_USE_EOSPAC
       >{};
// eos's that get relativistic and unit system modifiers
static constexpr const auto partial_eos_list =
    tl<IdealGas
#ifdef SPINER_USE_HDF
       ,
       SpinerEOSDependsRhoT, SpinerEOSDependsRhoSie, StellarCollapse
#endif // SPINER_USE_HDF
       >{};
// modifiers that get applied to all eos's
static constexpr const auto apply_to_all = al<ScaledEOS, ShiftedEOS>{};
// modifiers thet get applied to a subset of eos's
static constexpr const auto apply_to_partial = al<UnitSystem, RelativisticEOS>{};
// variadic list of eos's with shifted or scaled modifiers
static constexpr const auto shifted =
    transform_variadic_list(full_eos_list, al<ShiftedEOS>{});
static constexpr const auto scaled =
    transform_variadic_list(full_eos_list, al<ScaledEOS>{});
// variadic list of Scaled<Shifted<T>>'s
static constexpr const auto scaled_of_shifted =
    transform_variadic_list(shifted, al<ScaledEOS>{});
// variadic list of UnitSystem<T>'s
static constexpr const auto unit_system =
    transform_variadic_list(partial_eos_list, al<UnitSystem>{});
// variadic list of Relativistic<T>'s
static constexpr const auto relativistic =
    transform_variadic_list(partial_eos_list, al<RelativisticEOS>{});
// relativistic and unit system modifiers
static constexpr const auto unit_or_rel =
    transform_variadic_list(partial_eos_list, apply_to_partial);
// create combined list
static constexpr const auto combined_list_1 = singularity::detail::concat(
    full_eos_list, shifted, scaled, scaled_of_shifted, unit_or_rel);
// make a ramped eos of everything
static constexpr const auto ramped_all =
    transform_variadic_list(combined_list_1, al<BilinearRampEOS>{});
// final combined list
static constexpr const auto combined_list = singularity::detail::concat(
    full_eos_list, shifted, scaled, scaled_of_shifted, unit_or_rel, ramped_all);
// a function that returns a Variant from a typelist
template <typename... Ts>
struct tl_to_Variant_struct {
  using vt = Variant<Ts...>;
};

template <typename... Ts>
constexpr auto tl_to_Variant(tl<Ts...>) {
  return tl_to_Variant_struct<Ts...>{};
}

// create the alias
using EOS = typename decltype(tl_to_Variant(combined_list))::vt;
} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_EOS_HPP_
