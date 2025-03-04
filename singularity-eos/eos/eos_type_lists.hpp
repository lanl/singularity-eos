//------------------------------------------------------------------------------
// © 2021-2024. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_EOS_EOS_TYPE_LISTS_HPP_
#define _SINGULARITY_EOS_EOS_EOS_TYPE_LISTS_HPP_

#include <cassert>
#include <cmath>
#include <cstdio>
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
#include <singularity-eos/eos/eos_models.hpp>

namespace singularity {

// recreate variadic list
template <typename... Ts>
using tl = singularity::variadic_utils::type_list<Ts...>;

template <template <typename> class... Ts>
using al = singularity::variadic_utils::adapt_list<Ts...>;

// transform variadic list: applies modifiers to eos's
using singularity::variadic_utils::transform_variadic_list;

// all eos's
static constexpr const auto full_eos_list =
    tl<IdealGas, Gruneisen, Vinet, MGUsup, PowerMG, JWL, DavisReactants, DavisProducts
#ifdef SINGULARITY_USE_V_AND_V_EOS
       ,
       SAP_Polynomial, NobleAbel, CarnahanStarling, StiffGas
#endif // SINGULARITY_USE_V_AND_V_EOS
#ifdef SINGULARITY_USE_SPINER_WITH_HDF5
#ifdef SINGULARITY_USE_HELMHOLTZ
       ,
       Helmholtz
#endif // SINGULARITY_USE_HELMHOLTZ
       ,
       SpinerEOSDependsRhoT, SpinerEOSDependsRhoSie
#endif // SINGULARITY_USE_SPINER_WITH_HDF5
#ifdef SINGULARITY_USE_EOSPAC
       ,
       EOSPAC
#endif // SINGULARITY_USE_EOSPAC
       >{};
// modifiers that get applied to all eos's
static constexpr const auto apply_to_all = al<ScaledEOS, ShiftedEOS>{};
// variadic list of eos's with shifted or scaled modifiers
static constexpr const auto shifted =
    transform_variadic_list(full_eos_list, al<ShiftedEOS>{});
static constexpr const auto scaled_1 =
    transform_variadic_list(full_eos_list, al<ScaledEOS>{});
// variadic list of Scaled<Shifted<T>>'s
static constexpr const auto scaled_of_shifted =
    transform_variadic_list(shifted, al<ScaledEOS>{});
// combined list of all scaled EOS
static constexpr const auto scaled =
    singularity::variadic_utils::concat(scaled_1, scaled_of_shifted);
// create combined list
static constexpr const auto combined_list_1 =
    singularity::variadic_utils::concat(full_eos_list, shifted, scaled);
// make a ramped eos of everything
static constexpr const auto ramped_all =
    transform_variadic_list(combined_list_1, al<BilinearRampEOS>{});
// final combined list
static constexpr const auto combined_list =
    singularity::variadic_utils::concat(combined_list_1, ramped_all);

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_EOS_TYPE_LISTS_HPP_
