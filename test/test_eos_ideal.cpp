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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/indexable_types.hpp>
#include <singularity-eos/eos/eos.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <test/eos_unit_test_helpers.hpp>

using singularity::IdealElectrons;
using singularity::IdealGas;
using singularity::MeanAtomicProperties;
using EOS = singularity::Variant<IdealGas, IdealElectrons>;

using singularity::IndexableTypes::MeanIonizationState;
using Lambda_t = singularity::IndexerUtils::VariadicIndexer<MeanIonizationState>;

SCENARIO("Ideal gas entropy", "[IdealGas][Entropy][GibbsFreeEnergy]") {
  GIVEN("Parameters for an ideal gas with entropy reference states") {
    // Create ideal gas EOS ojbect
    constexpr Real Cv = 5.0;
    constexpr Real gm1 = 0.4;
    constexpr Real EntropyT0 = 100;
    constexpr Real EntropyRho0 = 1e-03;
    EOS host_eos = IdealGas(gm1, Cv, EntropyT0, EntropyRho0);
    THEN("The entropy at the reference state should be zero") {
      auto entropy = host_eos.EntropyFromDensityTemperature(EntropyRho0, EntropyT0);
      INFO("Entropy should be zero but it is " << entropy);
      CHECK(isClose(entropy, 0.0, 1.e-14));
    }
    GIVEN("A state at the reference temperature and a density whose cube root is the "
          "reference density") {
      constexpr Real T = EntropyT0;
      constexpr Real rho = 0.1; // rho**3 = EntropyRho0
      THEN("The entropy should be 2. / 3. * gm1 * Cv * log(EntropyRho0)") {
        const Real entropy_true = 2. / 3. * gm1 * Cv * log(EntropyRho0);
        auto entropy = host_eos.EntropyFromDensityTemperature(rho, T);
        INFO("Entropy: " << entropy << "  True entropy: " << entropy_true);
        CHECK(isClose(entropy, entropy_true, 1e-12));

        AND_THEN("The free energy agrees") {
          const Real sie = host_eos.InternalEnergyFromDensityTemperature(rho, T);
          const Real P = host_eos.PressureFromDensityTemperature(rho, T);
          const Real G_true = sie + (P / rho) - T * entropy;
          const Real GT = host_eos.GibbsFreeEnergyFromDensityTemperature(rho, T);
          CHECK(isClose(GT, G_true, 1e-12));
          const Real Gsie = host_eos.GibbsFreeEnergyFromDensityInternalEnergy(rho, sie);
          CHECK(isClose(Gsie, G_true, 1e-12));
        }
      }
    }
    GIVEN("A state at the reference density and a temperature whose square is the "
          "reference temperature") {
      constexpr Real T = 10; // T**2 = EntropyT0
      constexpr Real rho = EntropyRho0;
      THEN("The entropy should be -1. / 2. * Cv * log(EntropyT0)") {
        const Real entropy_true = -1. / 2. * Cv * log(EntropyT0);
        auto entropy = host_eos.EntropyFromDensityTemperature(rho, T);
        INFO("Entropy: " << entropy << "  True entropy: " << entropy_true);
        CHECK(isClose(entropy, entropy_true, 1e-12));
      }
    }
  }
}

SCENARIO("Ideal gas mean atomic properties",
         "[IdealGas][MeanAtomicMass][MeanAtomicNumber]") {
  constexpr Real Cv = 5.0;
  constexpr Real gm1 = 0.4;
  constexpr Real Abar = 4.0; // Helium
  constexpr Real Zbar = 2.0;
  const MeanAtomicProperties azbar(Abar, Zbar);
  GIVEN("An ideal gas initialized with mean atomic poroperties") {
    EOS host_eos = IdealGas(gm1, Cv, azbar);
    WHEN("We evaluate it on host") {
      Real Ab_eval = host_eos.MeanAtomicMass();
      Real Zb_eval = host_eos.MeanAtomicNumber();
      THEN("We get the right answer") {
        REQUIRE(isClose(Ab_eval, Abar, 1e-12));
        REQUIRE(isClose(Zb_eval, Zbar, 1e-12));
      }
    }
    WHEN("We evaluate it on device, using a loop") {
      constexpr int N = 100;
      auto device_eos = host_eos.GetOnDevice();
      int nwrong = 0;
      portableReduce(
          "Check mean atomic number", 0, N,
          PORTABLE_LAMBDA(const int i, int &nw) {
            Real rho = i;
            Real T = 100.0 * i;
            Real Ab_eval = device_eos.MeanAtomicMassFromDensityTemperature(rho, T);
            Real Zb_eval = device_eos.MeanAtomicNumberFromDensityTemperature(rho, T);
            bool needs_zbar = device_eos.NeedsLambda(MeanIonizationState());
            nw += !(isClose(Ab_eval, Abar, 1e-12)) + !(isClose(Zb_eval, Zbar, 1e-12));
            nw += needs_zbar;
          },
          nwrong);
      REQUIRE(nwrong == 0);
      device_eos.Finalize();
    }
    host_eos.Finalize();
  }
}

SCENARIO("Ideal gas density energy from prssure temperature",
         "[IdealGas][DensityEnergyFromPressureTemperature]") {
  constexpr Real Cv = 5.0;
  constexpr Real gm1 = 0.4;
  GIVEN("An ideal gas") {
    EOS host_eos = IdealGas(gm1, Cv);
    auto device_eos = host_eos.GetOnDevice();
    WHEN("We compute density and energy from pressure and temperature") {
      constexpr int N = 100;
      int nwrong = 0;
      portableReduce(
          "Check density energy from pressure temperature", 1, N,
          PORTABLE_LAMBDA(const int i, int &nw) {
            Real rho = i;
            Real T = 100.0 * i;
            nw += !CheckRhoSieFromPT(device_eos, rho, T);
          },
          nwrong);
      THEN("There are no errors") { REQUIRE(nwrong == 0); }
    }
  }
}

// A non-standard pattern where we do a reduction
class CheckPofRE {
 public:
  CheckPofRE(Real *P, Real *rho, Real *sie, int N) : P_(P), rho_(rho), sie_(sie), N_(N) {}
  template <typename T>
  void operator()(const T &eos) {
    Real *P = P_;
    Real *rho = rho_;
    Real *sie = sie_;
    portableReduce(
        "MyCheckPofRE", 0, N_,
        PORTABLE_LAMBDA(const int i, int &nw) {
          nw += !(isClose(P[i],
                          eos.PressureFromDensityInternalEnergy(rho[i], sie[i], nullptr),
                          1e-15));
        },
        nwrong);
  }

  // must be initialized to zero because portableReduce is a simple for loop on host
  int nwrong = 0;

 private:
  Real *P_;
  Real *rho_;
  Real *sie_;
  int N_;
};
SCENARIO("Ideal gas vector Evaluate call", "[IdealGas][Evaluate]") {
  GIVEN("An ideal gas, and some device memory") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.5;
    constexpr int N = 100;
    constexpr Real rhomin = 1;
    constexpr Real rhomax = 11;
    constexpr Real drho = (rhomax - rhomin) / (N - 1);
    constexpr Real siemin = 1;
    constexpr Real siemax = 11;
    constexpr Real dsie = (siemax - siemin) / (N - 1);

    EOS eos_host = IdealGas(gm1, Cv);
    auto eos_device = eos_host.GetOnDevice();

    Real *P = (Real *)PORTABLE_MALLOC(N * sizeof(Real));
    Real *rho = (Real *)PORTABLE_MALLOC(N * sizeof(Real));
    Real *sie = (Real *)PORTABLE_MALLOC(N * sizeof(Real));

    WHEN("We set P, rho, sie via the scalar API") {
      portableFor(
          "Set rho and sie", 0, N, PORTABLE_LAMBDA(const int i) {
            rho[i] = rhomin + drho * i; // just some diagonal in the 2d rho/sie space
            sie[i] = siemin + dsie * i;
            P[i] = eos_device.PressureFromDensityInternalEnergy(rho[i], sie[i]);
          });
      THEN("The vector Evaluate API can be used to compare") {
        CheckPofRE my_op(P, rho, sie, N);
        eos_device.EvaluateHost(my_op);
        REQUIRE(my_op.nwrong == 0);
      }
    }

    eos_device.Finalize();
    PORTABLE_FREE(P);
    PORTABLE_FREE(rho);
    PORTABLE_FREE(sie);
  }
}

struct Dummy {};
SCENARIO("Ideal gas serialization", "[IdealGas][Serialization]") {
  GIVEN("An ideal gas object on host and a variant on host") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.5;
    IdealGas eos_bare(gm1, Cv);
    EOS eos_variant = IdealGas(gm1, Cv);

    THEN("They both report zero dynamic memory size") {
      REQUIRE(eos_bare.DynamicMemorySizeInBytes() == 0);
      REQUIRE(eos_variant.DynamicMemorySizeInBytes() == 0);
    }

    THEN("They both report sizes larger than a trivial struct, such that the eos variant "
         "size >= eos_bare size") {
      REQUIRE(eos_bare.SerializedSizeInBytes() > sizeof(Dummy));
      REQUIRE(eos_variant.SerializedSizeInBytes() > sizeof(Dummy));
      REQUIRE(eos_variant.SerializedSizeInBytes() >= eos_bare.SerializedSizeInBytes());
    }

    WHEN("We serialize each") {
      auto [size_bare, data_bare] = eos_bare.Serialize();
      auto [size_var, data_var] = eos_variant.Serialize();

      THEN("The reported sizes are what we expect") {
        REQUIRE(size_bare == eos_bare.SerializedSizeInBytes());
        REQUIRE(size_var == eos_variant.SerializedSizeInBytes());
      }

      THEN("We can de-serialize new objects from them") {
        IdealGas new_bare;
        new_bare.DeSerialize(data_bare);

        EOS new_variant;
        new_variant.DeSerialize(data_var);

        AND_THEN("The bare eos has the right Cv and Gruneisen params") {
          REQUIRE(new_bare.SpecificHeatFromDensityTemperature(1.0, 1.0) == Cv);
          REQUIRE(new_bare.GruneisenParamFromDensityTemperature(1.0, 1.0) == gm1);
        }

        AND_THEN("The variant has the right type") {
          REQUIRE(new_variant.IsType<IdealGas>());
        }

        AND_THEN("The bare eos has the right Cv and Gruneisen params") {
          REQUIRE(new_variant.SpecificHeatFromDensityTemperature(1.0, 1.0) == Cv);
          REQUIRE(new_variant.GruneisenParamFromDensityTemperature(1.0, 1.0) == gm1);
        }
      }

      // cleanup
      free(data_bare);
      free(data_var);
    }

    // cleanup
    eos_bare.Finalize();
    eos_variant.Finalize();
  }
}

SCENARIO("Ideal electron gas", "[IdealGas][IdealEelctrons]") {
  GIVEN("An ideal electron gas from partially ionized iron") {
    constexpr Real Abar = 26;
    constexpr Real Zbar = 55.8;
    constexpr Real rho = 1;
    constexpr Real T = 4000;

    MeanAtomicProperties AZbar(Abar, Zbar);
    EOS eos = IdealElectrons(AZbar);

    THEN("The gruneisen coefficient is for 3 DOF") {
      Real gm1 = eos.GruneisenParamFromDensityTemperature(rho, T);
      Real gamma = gm1 + 1;
      REQUIRE(isClose(gamma, 5. / 3., 1e-12));
    }

    WHEN("We evaluate the specific heat for different partial ionizations") {
      Real lambda[1] = {1};
      const Real cv1 = eos.SpecificHeatFromDensityTemperature(rho, T, lambda);
      int nwrong = 0;
      constexpr int N = 55;
      portableReduce(
          "Check Cv vs Z", 2, N,
          PORTABLE_LAMBDA(const int i, int &nw) {
            Lambda_t ll;
            ll[MeanIonizationState()] = static_cast<Real>(i);
            Real Cv = eos.SpecificHeatFromDensityTemperature(rho, T, ll);
            if (!isClose(Cv, i * cv1, 1e-12)) nw += 1;
            bool needs_zbar = eos.NeedsLambda(MeanIonizationState());
            nw += !needs_zbar;
          },
          nwrong);
      THEN("The specific heat should scale linearly with the ionization state") {
        REQUIRE(nwrong == 0);
      }
    }

    WHEN("We compute density and energy from pressure and temperature") {
      constexpr int N = 100;
      int nwrong = 0;
      portableReduce(
          "Check density energy from pressure temperature", 1, N,
          PORTABLE_LAMBDA(const int i, int &nw) {
            Real ll[1] = {static_cast<Real>(i)};
            nw += !CheckRhoSieFromPT(eos, rho, T, ll);
          },
          nwrong);
      THEN("There are no errors") { REQUIRE(nwrong == 0); }
    }
  }
}
