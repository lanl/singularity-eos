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

#include <singularity-eos/base/variadic_utils.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_builder.hpp>
#include <singularity-eos/eos/eos_variant.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

using singularity::Gruneisen;
using singularity::IdealGas;
using singularity::ScaledEOS;
using singularity::ShiftedEOS;
using Var_t = singularity::Variant<IdealGas, ShiftedEOS<IdealGas>>;

namespace variadic_utils = singularity::variadic_utils;
namespace EOSBuilder = singularity::EOSBuilder;

SCENARIO("Conjunction and disjunction work properly", "[VariadicUtils]") {
  constexpr bool true_and_false = variadic_utils::conjunction<true, false>::value;
  constexpr bool true_or_false = variadic_utils::disjunction<true, false>::value;
  REQUIRE(true_and_false == false);
  REQUIRE(true_or_false == true);
}

SCENARIO("Contains works as intended", "[VariadicUtils]") {
  REQUIRE(variadic_utils::contains_v<IdealGas, IdealGas, ShiftedEOS<IdealGas>>());
  REQUIRE(!variadic_utils::contains_v<Gruneisen, IdealGas, ShiftedEOS<IdealGas>>());
}

SCENARIO("EOS Variant can check its own modifiability",
         "[VariadicUtils][EOSBuilder][EOSBase][EOSVariant]") {
  constexpr Real Cv = 2.0;
  constexpr Real gm1 = 0.5;
  Var_t v = IdealGas(gm1, Cv);
  REQUIRE(v.ModifiedInVariant<ShiftedEOS>());
  REQUIRE(!v.ModifiedInVariant<ScaledEOS>());
}

SCENARIO("IsModifiable works as intended", "[VariadicUtils][EOSBuilder]") {
  using namespace EOSBuilder;
  const Var_t v{};
  REQUIRE(is_modifiable<ShiftedEOS, IdealGas>(v));
  REQUIRE(!is_modifiable<ShiftedEOS, Gruneisen>(v));
  REQUIRE(!is_modifiable<ScaledEOS, IdealGas>(v));
  REQUIRE(!is_modifiable<ScaledEOS, Gruneisen>(v));
}

SCENARIO("Finding a unique set of types in a type list") {
  WHEN("A type list contains duplicate types") {
    using U0 = variadic_utils::unique_types_list<int, float, int, const float, int &>;
    THEN("The unique type list removes duplicates and preserves order") {
      static_assert(
          std::is_same<U0,
                       variadic_utils::type_list<int, float, const float, int &>>::value,
          "order-preserving unique");
    }
  }
  WHEN("A type list contains types that are duplicates when `const` and reference "
       "qualifiers are removed") {
    using U0 =
        variadic_utils::unique_decayed_types_list<int, int &, const int, const float>;
    static_assert(std::is_same<U0, variadic_utils::type_list<int, float>>::value,
                  "order-preserving unique");

    static_assert(std::is_same<variadic_utils::front_t<U0>, int>::value,
                  "first type in unique list");
  }
}
