//------------------------------------------------------------------------------
// © 2025. Triad National Security, LLC. All rights reserved.  This
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

#include <array>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/indexable_types.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <test/eos_unit_test_helpers.hpp>

using namespace singularity::IndexerUtils;
using namespace singularity::IndexableTypes;

using Lambda_t = VariadicIndexer<MeanIonizationState, ElectronFraction, LogDensity>;

class ManualLambda_t {
 public:
  constexpr static size_t length = 3;
  ManualLambda_t() = default;
  const Real &operator[](const std::size_t idx) const { return data_[idx]; }
  Real &operator[](const std::size_t idx) { return data_[idx]; }
  Real &operator[](const MeanIonizationState &zbar) { return data_[0]; }
  Real &operator[](const ElectronFraction &Ye) { return data_[1]; }
  Real &operator[](const LogDensity &lRho) { return data_[2]; }

 private:
  std::array<Real, length> data_;
};

class NewManualLambda_t : public ManualLambda_t {
 public:
  // Enable recognition that this is type-indexable
  constexpr static bool is_type_indexable = true;
};

SCENARIO("IndexableTypes and VariadicIndexer", "[IndexableTypes][VariadicIndexer]") {
  GIVEN("A variadic indexer, filled with indices 0, 1, 2") {
    Lambda_t lambda;
    for (std::size_t i = 0; i < Lambda_t::size(); ++i) {
      lambda[i] = static_cast<Real>(i);
    }
    WHEN("We access by index") {
      THEN("We get the expected index") {
        for (std::size_t i = 0; i < Lambda_t::size(); ++i) {
          Real val = lambda[i];
          REQUIRE(val == static_cast<Real>(i));
        }
      }
    }
    WHEN("We access by name") {
      Real Zbar = lambda[MeanIonizationState()];
      Real Ye = lambda[ElectronFraction()];
      Real lRho = lambda[LogDensity()];
      THEN("We get the expected index") {
        REQUIRE(Zbar == static_cast<Real>(0));
        REQUIRE(Ye == static_cast<Real>(1));
        REQUIRE(lRho == static_cast<Real>(2));
      }
    }
    WHEN("We use the SafeMustGet functionality") {
      // Request a type that exists, but no integer index
      const Real Zbar = SafeMustGet<MeanIonizationState>(lambda);
      // Request a type that exists, but the wrong integer index
      Real Zbar_1 = SafeMustGet<MeanIonizationState>(lambda, 2);
      THEN("We get the correct values") {
        REQUIRE(Zbar == static_cast<Real>(0));
        REQUIRE_THAT(Zbar, Catch::Matchers::WithinRel(Zbar_1, 1.0e-14));
      }
      // This is probably a dumb check since Zbar_1 wasn't a Real* or a Real&
      AND_THEN("We don't modify the lambda values by changing the local values") {
        Zbar_1 = 5.3;
        REQUIRE_THAT(Zbar, Catch::Matchers::WithinRel(
                               SafeMustGet<MeanIonizationState>(lambda), 1.0e-14));
      }
    }
    WHEN("We use the SafeMustSet functionality") {
      // Request a type that exists, but no index
      Real Zbar = 5.567;
      SafeMustSet<MeanIonizationState>(lambda, Zbar);
      // Request a type that exists, but the wrong integer index
      Real lRho = 1.102;
      SafeMustSet<LogDensity>(lambda, 0, lRho);
      THEN("We get the correct values") {
        REQUIRE_THAT(Zbar, Catch::Matchers::WithinRel(
                               SafeMustGet<MeanIonizationState>(lambda), 1.0e-14));
        REQUIRE_THAT(
            lRho, Catch::Matchers::WithinRel(SafeMustGet<LogDensity>(lambda), 1.0e-14));
      }
    }
    WHEN("We use the SafeGet functionality") {
      constexpr Real unmodified = -1.0;
      Real destination = unmodified;
      WHEN("We request a type that exists") {
        const bool modified = SafeGet<MeanIonizationState>(lambda, 2, destination);
        THEN("The destination value will be modified") {
          CHECK(modified);
          REQUIRE(destination == lambda[MeanIonizationState{}]);
        }
      }
      WHEN("We request a type that doesn't exist") {
        const bool modified = SafeGet<TableStatus>(lambda, 2, destination);
        THEN("The destination value will remain UNmodified") {
          CHECK(!modified);
          REQUIRE(destination == unmodified);
        }
      }
      WHEN("A normal array-like lambda is used") {
        std::array<Real, Lambda_t::size()> lambda_arr{1, 2, 3};
        constexpr size_t my_index = 2;
        const bool modified = SafeGet<LogDensity>(lambda_arr, my_index, destination);
        THEN("The destination value will reflect the index from the array") {
          CHECK(modified);
          REQUIRE(destination == lambda_arr[my_index]);
        }
      }
    }
    WHEN("We use the SafeSet functionality") {
      constexpr Real new_value = -1.0;
      WHEN("We want to set a value for a type index that exists") {
        const bool modified = SafeSet<MeanIonizationState>(lambda, 2, new_value);
        THEN("The lambda index was modified") {
          CHECK(modified);
          REQUIRE(lambda[MeanIonizationState{}] == new_value);
        }
      }
      WHEN("We want to set a value for a type index that doesn't exist") {
        const bool modified = SafeSet<TableStatus>(lambda, 2, new_value);
        Lambda_t old_lambda;
        for (std::size_t i = 0; i < Lambda_t::size(); ++i) {
          old_lambda[i] = lambda[i];
        }
        THEN("None of the lambda values were modified") {
          CHECK(!modified);
          for (std::size_t i = 0; i < Lambda_t::size(); ++i) {
            INFO("i: " << i);
            CHECK(lambda[i] == old_lambda[i]);
          }
        }
      }
      WHEN("A normal array-like lambda is used") {
        std::array<Real, Lambda_t::size()> lambda_arr{4, 5, 6};
        constexpr size_t my_index = 1;
        const bool modified = SafeSet<LogDensity>(lambda_arr, my_index, new_value);
        THEN("The lambda value at the appropriate index has been modified") {
          CHECK(modified);
          REQUIRE(lambda_arr[my_index] == new_value);
        }
      }
    }
  }
  GIVEN("A normal array that does not support IndexableTypes") {
    constexpr size_t num_lambda = 4;
    std::array<Real, num_lambda> lambda{};
    for (size_t i = 0; i < num_lambda; i++) {
      lambda[i] = static_cast<Real>(i);
    }
    WHEN("We use the SafeMustGet functionality") {
      THEN("The type-based index is ignored and only the integer index is used") {
        for (size_t i = 0; i < num_lambda; i++) {
          INFO("i: " << i);
          const Real val = SafeMustGet<LogTemperature>(lambda, i);
          CHECK_THAT(lambda[i], Catch::Matchers::WithinRel(val, 1.0e-14));
        }
      }
    }
    WHEN("We use the SafeMustSet functionality") {
      THEN("The type-based index is ignored and only the integer index is used") {
        for (size_t i = 0; i < num_lambda; i++) {
          INFO("i: " << i);
          const Real val = i * i;
          SafeMustSet<LogDensity>(lambda, i, val);
          CHECK_THAT(lambda[i], Catch::Matchers::WithinRel(val, 1.0e-14));
        }
      }
    }
  }
}

SCENARIO("IndexableTypes and ManualLambda", "[IndexableTypes]") {
  GIVEN("A manually written indexer, filled with indices 0, 1, 2") {
    ManualLambda_t lambda;
    for (std::size_t i = 0; i < lambda.length; ++i) {
      lambda[i] = static_cast<Real>(i);
    }
    WHEN("We use the SafeGet functionality") {
      constexpr Real unmodified = -1.0;
      Real destination = unmodified;
      WHEN("We request a type that doesn't exist in the manual indexer") {
        constexpr size_t my_index = 1;
        const bool modified = SafeGet<TableStatus>(lambda, my_index, destination);
        THEN("The destination WILL be modified since the manual indexer doesn't have the "
             " `is_type_indexable` data member") {
          CHECK(modified);
          REQUIRE(destination == lambda[my_index]);
        }
      }
      WHEN("We define a new manual indexer that has the `is_type_indexable` "
           "data member") {
        NewManualLambda_t lambda_new;
        for (std::size_t i = 0; i < lambda.length; ++i) {
          lambda_new[i] = lambda[i];
        }
        WHEN("We request a type that doesn't exist in the manual indexer") {
          destination = unmodified;
          constexpr size_t my_index = 1;
          const bool modified = SafeGet<TableStatus>(lambda_new, my_index, destination);
          THEN("The destination will NOT be modified") {
            CHECK(!modified);
            REQUIRE(destination == unmodified);
          }
        }
      }
    }
    WHEN("We use the SafeSet functionality") {
      constexpr Real new_value = -1.0;
      WHEN("We want to set a value for a type index that doesn't exist") {
        constexpr size_t my_index = 1;
        const bool modified = SafeSet<TableStatus>(lambda, my_index, new_value);
        THEN("The lambda value WILL be modified since the manual indexer doesn't have "
             "the `is_type_indexable` data member") {
          CHECK(modified);
          REQUIRE(lambda[my_index] == new_value);
        }
      }
      WHEN("We define a new manual indexer that has the `is_type_indexable` data "
           "member") {
        NewManualLambda_t lambda_new;
        for (std::size_t i = 0; i < lambda.length; ++i) {
          lambda_new[i] = lambda[i];
        }
        WHEN("We request a type that doesn't exist in the manual indexer") {
          constexpr size_t my_index = 1;
          const bool modified = SafeSet<TableStatus>(lambda_new, my_index, new_value);
          THEN("The lambda will NOT be modified") {
            CHECK(!modified);
            for (std::size_t i = 0; i < lambda.length; ++i) {
              INFO("i: " << i);
              CHECK(lambda_new[i] == lambda[i]);
            }
          }
        }
      }
    }
  }
}
