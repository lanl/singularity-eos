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

using namespace singularity::IndexerUtils;
using namespace singularity::IndexableTypes;

using Lambda_t = VariadicIndexer<MeanIonizationState, ElectronFraction, LogDensity>;

class ManualLambda_t {
 public:
  ManualLambda_t() = default;
  Real &operator[](const std::size_t idx) { return data_[idx]; }
  Real &operator[](const MeanIonizationState &zbar) { return data_[0]; }
  Real &operator[](const ElectronFraction &Ye) { return data_[1]; }
  Real &operator[](const LogDensity &lRho) { return data_[2]; }

 private:
  std::array<Real, 3> data_;
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
    WHEN("We use the Get functionality") {
      // Request a type that exists, but an incorrect index
      Real Zbar = Get<MeanIonizationState>(lambda, 2);
      // Request a type that doesn't exist but an index that does
      Real lRho = Get<LogTemperature>(lambda, 2);
      THEN("We get the correct values") {
        REQUIRE(Zbar == static_cast<Real>(0));
        REQUIRE(lRho == static_cast<Real>(2));
      }
    }
  }
}

SCENARIO("IndexableTypes and ManualLambda", "[IndexableTypes]") {
  GIVEN("A manually written indexer, filled with indices 0, 1, 2") {
    ManualLambda_t lambda;
    for (std::size_t i = 0; i < 3; ++i) {
      lambda[i] = static_cast<Real>(i);
    }
    WHEN("We use the Get functionality") {
      // Request a type that exists, but an incorrect index
      Real Zbar = Get<MeanIonizationState>(lambda, 2);
      // Request a type that doesn't exist but an index that does
      Real lRho = Get<LogTemperature>(lambda, 2);
      THEN("We get the correct values") {
        REQUIRE(Zbar == static_cast<Real>(0));
        REQUIRE(lRho == static_cast<Real>(2));
      }
    }
  }
}
