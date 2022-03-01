//------------------------------------------------------------------------------
// © 2021. Triad National Security, LLC. All rights reserved.  This
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
#include <array>
#include <iostream> // debug

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <root-finding-1d/root_finding.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_builder.hpp>

#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"

using singularity::EOS;
using singularity::IdealGas;
using singularity::ScaledEOS;
using singularity::ShiftedEOS;

#ifdef SPINER_USE_HDF
using singularity::SpinerEOSDependsRhoSie;
using singularity::SpinerEOSDependsRhoT;
#endif

#ifdef SINGULARITY_USE_EOSPAC
using singularity::EOSPAC;
#endif

namespace EOSBuilder = singularity::EOSBuilder;
namespace thermalqs = singularity::thermalqs;

const std::string eosName = "../materials.sp5";
const std::string airName = "air";
const std::string steelName = "stainless steel 347";

#ifdef SPINER_USE_HDF
#ifdef SINGULARITY_TEST_SESAME
constexpr int steelID = 4272;
constexpr int airID = 5030;
constexpr int DTID = 5267;
constexpr int gID = 2700;
constexpr Real ev2k = 1.160451812e4;
#endif // SINGULARITY_TEST_SESAME
#endif // SPINER_USE_HDF

constexpr Real EPS = 5e-2; // within a few percent
PORTABLE_INLINE_FUNCTION bool isClose(Real a, Real b) {
  return fabs(b - a) / (fabs(a + b) + 1e-20) <= EPS;
}

PORTABLE_INLINE_FUNCTION Real myAtan(Real x, Real shift, Real scale, Real offset) {
  return scale * atan(x - shift) + offset;
}

/*
Notes for creating tests within the Catch2 framework:

Within each "SCENARIO", the "GIVEN", "WHEN", "THEN" , etc. blocks all are
aliases for "SECTION" and the Catch2 documentation states that all "SECTION"s
must have unique names. Even if the "SECTION" is a subsection of another
"SECTION", it still needs to have a completley unique string associated with it.
If this doesn't happen, the resulting namespace conflict can cause undefined
behavior. This can include running bilions of assertions and unreproducible test
results.

I (JHP) haven't tested whether this also happens across SCENARIOS, but Jonah
thinks that is probably okay.
*/

SCENARIO("Rudimentary test of the root finder", "[RootFinding1D]") {

  GIVEN("A root counts object") {
    using namespace RootFinding1D;
    RootCounts counts;

    THEN("A root can be found for shift = 1, scale = 2, offset = 0.5") {
      int ntimes = 100;
      Real guess = 0;
      Real root;
      Status status;
      Real shift = 1;
      Real scale = 2;
      Real offset = 0.5;

      auto f = PORTABLE_LAMBDA(const Real x) { return myAtan(x, shift, scale, offset); };
      Status *statusesp = (Status *)PORTABLE_MALLOC(ntimes * sizeof(Status));
      Real *rootsp = (Real *)PORTABLE_MALLOC(ntimes * sizeof(Real));
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<Status *, Kokkos::MemoryTraits<Kokkos::Unmanaged>> statuses(statusesp,
                                                                               ntimes);
      Kokkos::View<Real *, Kokkos::MemoryTraits<Kokkos::Unmanaged>> roots(rootsp, ntimes);
#else
      PortableMDArray<Status> statuses;
      PortableMDArray<Real> roots;
      statuses.NewPortableMDArray(statusesp, ntimes);
      roots.NewPortableMDArray(rootsp, ntimes);
#endif
      portableFor(
          "find roots", 0, ntimes, PORTABLE_LAMBDA(const int i) {
            RootCounts per_thread_counts;
            statuses(i) =
                findRoot(f, 0, guess, -1, 3, 1e-10, 1e-10, roots(i), per_thread_counts);
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<Status> s_copy(statuses, 0);
      Kokkos::View<Real> r_copy(roots, 0);
      Kokkos::deep_copy(root, r_copy);
      Kokkos::deep_copy(status, s_copy);
#else
      status = statuses(ntimes - 1);
      root = roots(ntimes - 1);
#endif

      PORTABLE_FREE(statusesp);
      PORTABLE_FREE(rootsp);
      REQUIRE(status == Status::SUCCESS);
      REQUIRE(isClose(root, 0.744658));
      REQUIRE(100. * counts[counts.more()] / counts.total() <= 10);
    }
  }
}

SCENARIO("EOS Builder and Modifiers", "[EOSBuilder],[Modifiers][IdealGas]") {

  GIVEN("Parameters for a shifted and scaled ideal gas") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.5;
    constexpr Real scale = 2.0;
    constexpr Real shift = 0.1;
    constexpr Real rho = 2.0;
    constexpr Real sie = 0.5;
    WHEN("We construct a shifted, scaled IdealGas by hand") {
      IdealGas a = IdealGas(gm1, Cv);
      ScaledEOS<IdealGas> b = ScaledEOS<IdealGas>(std::move(a), scale);
      EOS eos = ShiftedEOS<ScaledEOS<IdealGas>>(std::move(b), shift);
      THEN("The shift and scale parameters pass through correctly") {

        REQUIRE(eos.PressureFromDensityInternalEnergy(rho, sie) == 0.4);
      }
    }
    WHEN("We use the EOSBuilder") {
      EOSBuilder::EOSType type = EOSBuilder::EOSType::IdealGas;
      EOSBuilder::modifiers_t modifiers;
      EOSBuilder::params_t base_params, shifted_params, scaled_params;
      base_params["Cv"].emplace<Real>(Cv);
      base_params["gm1"].emplace<Real>(gm1);
      shifted_params["shift"].emplace<Real>(shift);
      scaled_params["scale"].emplace<Real>(scale);
      modifiers[EOSBuilder::EOSModifier::Shifted] = shifted_params;
      modifiers[EOSBuilder::EOSModifier::Scaled] = scaled_params;
      EOS eos = EOSBuilder::buildEOS(type, base_params, modifiers);
      THEN("The shift and scale parameters pass through correctly") {
        REQUIRE(eos.PressureFromDensityInternalEnergy(rho, sie) == 0.4);
      }
    }
  }
}

SCENARIO("Relativistic EOS", "[EOSBuilder][RelativisticEOS][IdealGas]") {
  GIVEN("Parameters for an ideal gas") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.5;
    WHEN("We construct a relativistic IdealGas with EOSBuilder") {
      EOSBuilder::EOSType type = EOSBuilder::EOSType::IdealGas;
      EOSBuilder::modifiers_t modifiers;
      EOSBuilder::params_t base_params, relativity_params;
      base_params["Cv"].emplace<Real>(Cv);
      base_params["gm1"].emplace<Real>(gm1);
      relativity_params["cl"].emplace<Real>(1.0);
      modifiers[EOSBuilder::EOSModifier::Relativistic] = relativity_params;
      EOS eos = EOSBuilder::buildEOS(type, base_params, modifiers);
      THEN("The EOS has finite sound speeds") {
        constexpr Real rho = 1e3;
        constexpr Real sie = 1e3;
        Real bmod = eos.BulkModulusFromDensityInternalEnergy(rho, sie);
        Real cs2 = bmod / rho;
        REQUIRE(cs2 < 1);
      }
    }
  }
}

SCENARIO("EOS Unit System", "[EOSBuilder][UnitSystem][IdealGas]") {
  GIVEN("Parameters for an ideal gas") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.5;
    EOSBuilder::EOSType type = EOSBuilder::EOSType::IdealGas;
    EOSBuilder::modifiers_t modifiers;
    EOSBuilder::params_t base_params, units_params;
    base_params["Cv"].emplace<Real>(Cv);
    base_params["gm1"].emplace<Real>(gm1);
    GIVEN("Units with a thermal unit system") {
      constexpr Real rho_unit = 1e1;
      constexpr Real sie_unit = 1e-1;
      constexpr Real temp_unit = 123;
      WHEN("We construct an IdealGas with EOSBuilder") {
        units_params["rho_unit"].emplace<Real>(rho_unit);
        units_params["sie_unit"].emplace<Real>(sie_unit);
        units_params["temp_unit"].emplace<Real>(temp_unit);
        modifiers[EOSBuilder::EOSModifier::UnitSystem] = units_params;
        EOS eos = EOSBuilder::buildEOS(type, base_params, modifiers);
        THEN("Units cancel out for an ideal gas") {
          Real rho = 1e3;
          Real sie = 1e3;
          Real P = eos.PressureFromDensityInternalEnergy(rho, sie);
          Real Ptrue = gm1 * rho * sie;
          REQUIRE(std::abs(P - Ptrue) / Ptrue < 1e-3);
        }
      }
    }
    GIVEN("Units with length and time units") {
      constexpr Real time_unit = 456;
      constexpr Real length_unit = 1e2;
      constexpr Real mass_unit = 1e6;
      constexpr Real temp_unit = 789;
      WHEN("We construct an IdealGas with EOSBuilder") {
        units_params["use_length_time"].emplace<bool>(true);
        units_params["time_unit"].emplace<Real>(time_unit);
        units_params["length_unit"].emplace<Real>(length_unit);
        units_params["mass_unit"].emplace<Real>(mass_unit);
        units_params["temp_unit"].emplace<Real>(temp_unit);
        modifiers[EOSBuilder::EOSModifier::UnitSystem] = units_params;
        EOS eos = EOSBuilder::buildEOS(type, base_params, modifiers);
        THEN("Units cancel out for an ideal gas") {
          Real rho = 1e3;
          Real sie = 1e3;
          Real P = eos.PressureFromDensityInternalEnergy(rho, sie);
          Real Ptrue = gm1 * rho * sie;
          REQUIRE(std::abs(P - Ptrue) / Ptrue < 1e-3);
        }
      }
    }
  }
}

SCENARIO("Vector EOS", "[VectorEOS][IdealGas]") {

  GIVEN("Parameters for an ideal gas") {
    // Create ideal gas EOS ojbect
    constexpr Real Cv = 5.0;
    constexpr Real gm1 = 0.4;
    EOS eos = IdealGas(gm1, Cv);

    GIVEN("Energies and densities") {
      // Input arrays and pointers
      constexpr int num = 3;
      constexpr std::array<Real, num> density {1.0, 2.0, 5.0};
      constexpr std::array<Real, num> energy {5.0, 10.0, 15.0};
      std::array<Real*, num> lambdas;

      auto plambdas = lambdas.data();
      auto pdensity = density.data();
      auto penergy = energy.data();

      // Gold standard values
      constexpr std::array<Real, num> pressure_true {2.0, 8.0, 30.0};
      constexpr std::array<Real, num> temperature_true {1., 2., 3.};
      constexpr std::array<Real, num> bulkmodulus_true {2.8, 11.2, 42.};

      // Output arrays and pointers
      std::array<Real, num> temperature;
      std::array<Real, num> pressure;
      std::array<Real, num> heatcapacity;
      std::array<Real, num> bulkmodulus;
      std::array<Real, num> gruneisen;
      auto ptemperature = temperature.data();
      auto ppressure = pressure.data();
      auto pheatcapacity = heatcapacity.data();
      auto pbulkmodulus = bulkmodulus.data();
      auto pgruneisen = gruneisen.data();

      WHEN("A T(rho, e) lookup is performed") {
        eos.TemperatureFromDensityInternalEnergy(pdensity, penergy,
                                                 ptemperature, num, lambdas);
        THEN("The returned T(rho, e) should be equal to the true "
             "temperature") {
          for (int i; i < num; i++) {
            INFO("i: " << i << " Density: " << density[i] << " Energy: "
                 << energy[i]);
            REQUIRE(temperature[i] == Approx(temperature_true[i]));
          }
        }
      }

      WHEN("A P(rho, e) lookup is performed") {
        eos.PressureFromDensityInternalEnergy(pdensity, penergy, ppressure,
                                              num, lambdas);
        THEN("The returned P(rho, e) should be equal to the true pressure") {
          for (int i; i < num; i++) {
            INFO("i: " << i << " Density: " << density[i] << " Energy: "
                 << energy[i]);
            REQUIRE(pressure[i] == Approx(pressure_true[i]));
          }
        }
      }

      WHEN("A C_v(rho, e) lookup is performed") {
        eos.SpecificHeatFromDensityInternalEnergy(pdensity, penergy,
                                                  pheatcapacity, num, lambdas);
        THEN("The returned C_v(rho, e) should be constant") {
          for (int i; i < num; i++) {
            INFO("i: " << i << " Density: " << density[i] << " Energy: "
                 << energy[i]);
            REQUIRE(heatcapacity[i] == Approx(Cv));
          }
        }
      }

      WHEN("A B_S(rho, e) lookup is performed") {
        eos.BulkModulusFromDensityInternalEnergy(pdensity, penergy,
                                                 pbulkmodulus, num, lambdas);
        THEN("The returned B_S(rho, e) should be equal to the true bulk "
              "modulus") {
          for (int i; i < num; i++) {
            INFO("i: " << i << " Density: " << density[i] << " Energy: "
                 << energy[i]);
            REQUIRE(bulkmodulus[i] == Approx(bulkmodulus_true[i]));
          }
        }
      }

      WHEN("A Gamma(rho, e) lookup is performed") {
        eos.GruneisenParamFromDensityInternalEnergy(pdensity, penergy,
                                                    pgruneisen, num, lambdas);
        THEN("The returned Gamma(rho, e) should be constant") {
          for (int i; i < num; i++) {
            INFO("i: " << i << " Density: " << density[i] << " Energy: "
                 << energy[i]);
            REQUIRE(gruneisen[i] == Approx(gm1));
          }
        }
      }
    }
    GIVEN("Densities and temperatures") {
      // Input arrays and pointers
      constexpr int num = 3;
      constexpr std::array<Real, num> density {1.0, 2.0, 5.0};
      constexpr std::array<Real, num> temperature {50., 100., 150.};
      std::array<Real*, num> lambdas;
      auto plambdas = lambdas.data();
      auto pdensity = density.data();
      auto ptemperature = temperature.data();

      // Gold standard values
      constexpr std::array<Real, num> energy_true {250., 500., 750.};
      constexpr std::array<Real, num> pressure_true {100., 400., 1500.};
      constexpr std::array<Real, num> bulkmodulus_true {140., 560., 2100.};

      // Output arrays and pointers
      std::array<Real, num> energy;
      std::array<Real, num> pressure;
      std::array<Real, num> heatcapacity;
      std::array<Real, num> bulkmodulus;
      std::array<Real, num> gruneisen;
      auto penergy = energy.data();
      auto ppressure = pressure.data();
      auto pheatcapacity = heatcapacity.data();
      auto pbulkmodulus = bulkmodulus.data();
      auto pgruneisen = gruneisen.data();

      WHEN("A e(rho, T) lookup is performed") {
        eos.InternalEnergyFromDensityTemperature(pdensity, ptemperature,
                                                 penergy, num, lambdas);
        THEN("The returned e(rho, T) should be equal to the true energy") {
          for (int i; i < num; i++) {
            INFO("i: " << i << " Density: " << density[i] << " Temperature: "
                 << temperature[i]);
            REQUIRE(energy[i] == Approx(energy_true[i]));
          }
        }
      }

      WHEN("A P(rho, T) lookup is performed") {
        eos.PressureFromDensityTemperature(pdensity, ptemperature,
                                           ppressure, num, lambdas);
        THEN("The returned P(rho, T) should be equal to the true pressure") {
          for (int i; i < num; i++) {
            INFO("i: " << i << " Density: " << density[i] << " Temperature: "
                 << temperature[i]);
            REQUIRE(pressure[i] == Approx(pressure_true[i]));
          }
        }
      }

      WHEN("A C_v(rho, T) lookup is performed") {
        eos.SpecificHeatFromDensityTemperature(pdensity, ptemperature,
                                               pheatcapacity, num, lambdas);
        THEN("The returned C_v(rho, T) should be constant") {
          for (int i; i < num; i++) {
            INFO("i: " << i << " Density: " << density[i] << " Temperature: "
                 << temperature[i]);
            REQUIRE(heatcapacity[i] == Approx(Cv));
          }
        }
      }

      WHEN("A B_S(rho, T) lookup is performed") {
        eos.BulkModulusFromDensityTemperature(pdensity, ptemperature,
                                              pbulkmodulus, num, lambdas);
        THEN("The returned B_S(rho, T) should be equal to the true bulk "
              "modulus") {
          for (int i; i < num; i++) {
            INFO("i: " << i << " Density: " << density[i] << " Temperature: "
                 << temperature[i]);
            REQUIRE(bulkmodulus[i] == Approx(bulkmodulus_true[i]));
          }
        }
      }

      WHEN("A Gamma(rho, T) lookup is performed") {
        eos.GruneisenParamFromDensityTemperature(pdensity, ptemperature,
                                                 pgruneisen, num, lambdas);
        THEN("The returned Gamma(rho, T) should be constant") {
          for (int i; i < num; i++) {
            INFO("i: " << i << " Density: " << density[i] << " Temperature: "
                 << temperature[i]);
            REQUIRE(gruneisen[i] == Approx(gm1));
          }
        }
      }
    }
  }
}

#ifdef SPINER_USE_HDF
#ifdef SINGULARITY_TEST_SESAME
#ifdef SINGULARITY_USE_EOSPAC
SCENARIO("SpinerEOS depends on Rho and T", "[SpinerEOS],[DependsRhoT][EOSPAC]") {

  GIVEN("SpinerEOS and EOSPAC EOS for steel can be initialized with matid") {
    EOS steelEOS_host_polymorphic = SpinerEOSDependsRhoT(eosName, steelID);
    EOS steelEOS = steelEOS_host_polymorphic.GetOnDevice();
    auto steelEOS_host = steelEOS_host_polymorphic.get<SpinerEOSDependsRhoT>();

    EOS eospac = EOSPAC(steelID);

    THEN("The correct metadata is read in") {
      REQUIRE(steelEOS_host.matid() == steelID);
      REQUIRE(steelEOS_host.filename() == eosName);

      AND_THEN("We can get a reference density and temperature") {
        Real rho, T, sie, P, cv, bmod, dpde, dvdt;
        Real rho_pac, T_pac, sie_pac, P_pac, cv_pac, bmod_pac, dpde_pac, dvdt_pac;
        steelEOS_host.ValuesAtReferenceState(rho, T, sie, P, cv, bmod, dpde, dvdt);
        eospac.ValuesAtReferenceState(rho_pac, T_pac, sie_pac, P_pac, cv_pac, bmod_pac,
                                      dpde_pac, dvdt_pac);
        REQUIRE(isClose(rho, rho_pac));
        REQUIRE(isClose(T, T_pac));
      }

      // TODO: this needs to be a much more rigorous test
      AND_THEN("Quantities can be read from density and temperature") {
        const Real sie_pac = eospac.InternalEnergyFromDensityTemperature(1e0, 1e6);

        int nw_ie{0};
#ifdef PORTABILITY_STRATEGY_KOKKOS
        using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
        Kokkos::View<int, atomic_view> n_wrong_ie("wrong_ie");
#else
        PortableMDArray<int> n_wrong_ie(&nw_ie, 1);
#endif
        portableFor(
            "calc ie's", 0, 100, PORTABLE_LAMBDA(const int &i) {
              const double ie = steelEOS.InternalEnergyFromDensityTemperature(1e0, 1e6);
              if (!isClose(ie, sie_pac)) n_wrong_ie() += 1;
            });
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::deep_copy(nw_ie, n_wrong_ie);
#endif
        REQUIRE(nw_ie == 0);
      }

      AND_THEN("rho(P,T) correct for P=1atm, T=freezing") {
        Real T = 273;  // Kelvin
        Real P = 1e6;  // barye
        Real rho, sie; // output
        Real rho_pac, sie_pac;
        std::vector<Real> lambda(steelEOS_host_polymorphic.nlambda());
        steelEOS_host_polymorphic.DensityEnergyFromPressureTemperature(
            P, T, lambda.data(), rho, sie);
        eospac.DensityEnergyFromPressureTemperature(P, T, nullptr, rho_pac, sie_pac);
        REQUIRE(isClose(rho, rho_pac));
      }
    }
    // Failing to call finalize leads to a memory leak,
    // but otherwise behaviour is as expected.
    // It's possible to this automatically clean up with
    // some form of reference counting. If this is a priority,
    // we can re-examine.
    steelEOS_host_polymorphic.Finalize(); // host and device must be
                                          // finalized separately.
    steelEOS.Finalize();                  // cleans up memory on device.
  }

  GIVEN("SpinerEOS and EOSPAC for air can be initialized with matid") {
    SpinerEOSDependsRhoT airEOS_host = SpinerEOSDependsRhoT(eosName, airID);
    EOS airEOS = airEOS_host.GetOnDevice();
    EOS eospac = EOSPAC(airID);
    constexpr Real rho = 1e0;
    constexpr Real sie = 2.5e-4;
    THEN("We can get a reference state") {
      Real rho, T, sie, P, cv, bmod, dpde, dvdt;
      Real rho_pac, T_pac, sie_pac, P_pac, cv_pac, bmod_pac, dpde_pac, dvdt_pac;
      airEOS_host.ValuesAtReferenceState(rho, T, sie, P, cv, bmod, dpde, dvdt);
      eospac.ValuesAtReferenceState(rho_pac, T_pac, sie_pac, P_pac, cv_pac, bmod_pac,
                                    dpde_pac, dvdt_pac);
      REQUIRE(isClose(rho, rho_pac));
      REQUIRE(isClose(T, T_pac));
    }
    THEN("Quantities of rho and sie look sane on both host and device") {
      Real gm1_host = airEOS_host.GruneisenParamFromDensityInternalEnergy(rho, sie);
      Real T_host = airEOS_host.TemperatureFromDensityInternalEnergy(rho, sie);
      Real cv_host = airEOS_host.SpecificHeatFromDensityInternalEnergy(rho, sie);

      int nw_gm1{0}, nw_cv{0};
#ifdef PORTABILITY_STRATEGY_KOKKOS
      using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
      Kokkos::View<int, atomic_view> n_wrong_gm1("wrong_gm1");
      Kokkos::View<int, atomic_view> n_wrong_cv("wrong_cv");
#else
      PortableMDArray<int> n_wrong_gm1(&nw_gm1, 1);
      PortableMDArray<int> n_wrong_cv(&nw_cv, 1);
#endif
      portableFor(
          "calc gm1 and cv", 0, 100, PORTABLE_LAMBDA(const int &i) {
            const double gm1 = airEOS.GruneisenParamFromDensityInternalEnergy(rho, sie);
            const double cv = airEOS.SpecificHeatFromDensityInternalEnergy(rho, sie);
            if (!isClose(gm1, gm1_host)) n_wrong_gm1() += 1;
            if (!isClose(cv, cv_host)) n_wrong_cv() += 1;
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::deep_copy(nw_gm1, n_wrong_gm1);
      Kokkos::deep_copy(nw_cv, n_wrong_cv);
#endif
      REQUIRE(nw_gm1 == 0);
      REQUIRE(nw_cv == 0);
    }
    AND_THEN("P(rho, sie) correct for extrapolation regime") {
      Real rho = 1;
      Real sie = 2.43e16;
      Real P_pac = eospac.PressureFromDensityInternalEnergy(rho, sie);
      Real P_spi = airEOS.PressureFromDensityInternalEnergy(rho, sie);
      REQUIRE(isClose(P_pac, P_spi));
    }
    airEOS_host.Finalize();
    airEOS.Finalize();
  }

  GIVEN("EOS initialized with matid") {
    SpinerEOSDependsRhoT eos_spiner = SpinerEOSDependsRhoT(eosName, DTID);
    EOS eos_eospac = EOSPAC(DTID);
    Real P = 1e8;
    Real rho = 1.28e-3;
    THEN("Inversion for T(rho,P) works on host") {
      Real T, sie, cv, bmod;
      Real T_pac, sie_pac, cv_pac, bmod_pac;
      const unsigned long output =
          (thermalqs::temperature | thermalqs::specific_internal_energy |
           thermalqs::specific_heat | thermalqs::bulk_modulus);
      eos_spiner.FillEos(rho, T, sie, P, cv, bmod, output);
      eos_eospac.FillEos(rho, T_pac, sie_pac, P, cv_pac, bmod_pac, output);
      REQUIRE(isClose(T, T_pac));
      REQUIRE(isClose(sie, sie_pac));
      REQUIRE(isClose(cv, cv_pac));
    }
    eos_spiner.Finalize();
  }

  GIVEN("EOS initialized with matid") {
    SpinerEOSDependsRhoT eos_spiner = SpinerEOSDependsRhoT(eosName, gID);
    EOS eos_eospac = EOSPAC(gID);
    THEN("PT lookup works on the host") {
      Real P = 1e6;          // cgs
      Real T = 0.025 / ev2k; // K
      Real rho, sie;
      Real rho_pac, sie_pac;
      std::vector<Real> lambda(eos_spiner.nlambda());
      eos_spiner.DensityEnergyFromPressureTemperature(P, T, lambda.data(), rho, sie);
      eos_eospac.DensityEnergyFromPressureTemperature(P, T, lambda.data(), rho_pac,
                                                      sie_pac);
      REQUIRE(isClose(rho, rho_pac));
    }
    eos_spiner.Finalize();
  }
}

// Disabling these tests for now as the DependsRhoSie code is not well-maintained
SCENARIO("SpinerEOS depends on rho and sie", "[SpinerEOS],[DependsRhoSie]") {

  GIVEN("SpinerEOSes for steel can be initialised with matid") {
    SpinerEOSDependsRhoSie steelEOS_host(eosName, steelID);
    EOS steelEOS = steelEOS_host.GetOnDevice();
    THEN("The correct metadata is read in") {
      REQUIRE(steelEOS_host.matid() == steelID);
      REQUIRE(steelEOS_host.filename() == eosName);

      int nw_ie2{0}, nw_te2{0};
#ifdef PORTABILITY_STRATEGY_KOKKOS
      using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
      Kokkos::View<int, atomic_view> n_wrong_ie("wrong_ie");
      Kokkos::View<int, atomic_view> n_wrong_te("wrong_te");
#else
      PortableMDArray<int> n_wrong_ie(&nw_ie2, 1);
      PortableMDArray<int> n_wrong_te(&nw_te2, 1);
#endif
      portableFor(
          "calc ie's", 0, 100, PORTABLE_LAMBDA(const int &i) {
            const Real ie{steelEOS.InternalEnergyFromDensityTemperature(1e0, 1e6)};
            const Real te{steelEOS.TemperatureFromDensityInternalEnergy(1e0, 1e12)};
            if (!isClose(ie, 4.96415e13)) n_wrong_ie() += 1;
            if (!isClose(te, 75594.6)) n_wrong_te() += 1;
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::deep_copy(nw_ie2, n_wrong_ie);
      Kokkos::deep_copy(nw_te2, n_wrong_te);
#endif
      REQUIRE(nw_ie2 == 0);
      REQUIRE(nw_te2 == 0);
    }
    // this can be removed with with reference counting or other tricks
    steelEOS_host.Finalize(); // cleans up host memory
    steelEOS.Finalize();      // cleans up device memory
  }
}

SCENARIO("EOS Builder and SpinerEOS",
         "[SpinerEOS],[EOSBuilder],[GetOnDevice],[Finalize]") {
  GIVEN("Parameters for shift and scale") {
    constexpr Real shift = 0.0;
    constexpr Real scale = 1.0;
    WHEN("We construct a SpinerEOS with EOSBuilder") {
      EOSBuilder::EOSType type = EOSBuilder::EOSType::SpinerEOSDependsRhoT;
      EOSBuilder::modifiers_t modifiers;
      EOSBuilder::params_t base_params, shifted_params, scaled_params;
      base_params["filename"].emplace<std::string>(eosName);
      base_params["matid"].emplace<int>(steelID);
      shifted_params["shift"].emplace<Real>(shift);
      scaled_params["scale"].emplace<Real>(scale);
      modifiers[EOSBuilder::EOSModifier::Shifted] = shifted_params;
      modifiers[EOSBuilder::EOSModifier::Scaled] = scaled_params;
      EOS eosHost = EOSBuilder::buildEOS(type, base_params, modifiers);
      THEN("The EOS is consistent.") {
        REQUIRE(
            isClose(4.96416e13, eosHost.InternalEnergyFromDensityTemperature(1e0, 1e6)));
      }
      WHEN("We get an EOS on device") {
        EOS eosDevice = eosHost.GetOnDevice();
        THEN("EOS calls match raw access") {
          int nw_bm{0};
#ifdef PORTABILITY_STRATEGY_KOKKOS
          using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
          Kokkos::View<int, atomic_view> n_wrong_bm("wrong_bm");
#else
          PortableMDArray<int> n_wrong_bm(&nw_bm, 1);
#endif
          portableFor(
              "calc ie's", 0, 100, PORTABLE_LAMBDA(const int &i) {
                const Real bm{eosDevice.BulkModulusFromDensityTemperature(1e0, 1e6)};
                if (!isClose(bm, 2.55268e13)) n_wrong_bm() += 1;
              });
#ifdef PORTABILITY_STRATEGY_KOKKOS
          Kokkos::deep_copy(nw_bm, n_wrong_bm);
#endif
          REQUIRE(nw_bm == 0);
        }
        eosDevice.Finalize();
      }
      eosHost.Finalize();
    }
  }
}
#endif // SINGULARITY_USE_EOSPAC
#endif // SINGULARITY_TEST_SESAME

#ifdef SINGULARITY_TEST_STELLAR_COLLAPSE
SCENARIO("Stellar Collapse EOS", "[StellarCollapse][EOSBuilder]") {
  using singularity::IdealGas;
  using singularity::StellarCollapse;
  const std::string savename = "stellar_collapse_ideal_2.sp5";
  GIVEN("A stellar collapse EOS") {
    const std::string filename = "../stellar_collapse_ideal.h5";
    THEN("We can load the file") { // don't bother filtering bmod here.
      StellarCollapse sc(filename, false, false);
      AND_THEN("Some properties we expect for ideal gas hold") {
        Real lambda[2];
        Real rho, t, sie, p, cv, b, dpde, dvdt;
        sc.ValuesAtReferenceState(rho, t, sie, p, cv, b, dpde, dvdt, lambda);
        Real yemin = sc.YeMin();
        Real yemax = sc.YeMax();
        int N = 123;
        Real dY = (yemax - yemin) / (N + 1);
        for (int i = 0; i < N; ++i) {
          Real Ye = yemin + i * dY;
          lambda[0] = Ye;
          REQUIRE(isClose(sie, sc.InternalEnergyFromDensityTemperature(rho, t, lambda)));
        }
        Real rhomin = sc.rhoMin();
        Real rhomax = sc.rhoMax();
        Real drho = (rhomax - rhomin) / (N + 1);
        for (int i = 0; i < N; ++i) {
          Real rho = rhomin + i * drho;
          REQUIRE(isClose(sie, sc.InternalEnergyFromDensityTemperature(rho, t, lambda)));
        }
      }
      GIVEN("An Ideal Gas equation of state") {
        constexpr Real gamma = 1.4;
        constexpr Real mp = 1.67262171e-24;
        constexpr Real Cv = 1. / (mp * (gamma - 1)); // assumes cgs
        IdealGas ig(gamma - 1, Cv);
        auto ig_d = ig.GetOnDevice();
        THEN("The tabulated gamma Stellar Collapse and the gamma agree roughly") {
          Real yemin = sc.YeMin();
          Real yemax = sc.YeMax();
          Real tmin = sc.TMin();
          Real tmax = sc.TMax();
          Real ltmin = std::log10(tmin);
          Real ltmax = std::log10(tmax);
          Real lrhomin = sc.lRhoMin();
          Real lrhomax = sc.lRhoMax();
          auto sc_d = sc.GetOnDevice();

          int nwrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
          using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
          Kokkos::View<int, atomic_view> nwrong_d("wrong");
#else
          PortableMDArray<int> nwrong_d(&nwrong_h, 1);
#endif

          const int N = 123;
          const Real dY = (yemax - yemin) / (N + 1);
          const Real dlT = (ltmax - ltmin) / (N + 1);
          const Real dlR = (lrhomax - lrhomin) / (N + 1);
          portableFor(
              "fill eos", 0, N, 0, N, 0, N,
              PORTABLE_LAMBDA(const int &k, const int &j, const int &i) {
                Real lambda[2];
                Real Ye = yemin + k * dY;
                Real lT = ltmin + j * dlT;
                Real lR = lrhomin + i * dlR;
                Real T = std::pow(10., lT);
                Real R = std::pow(10., lR);
                Real e1, e2, p1, p2, cv1, cv2, b1, b2;
                unsigned long output = (singularity::thermalqs::pressure |
                                        singularity::thermalqs::specific_internal_energy |
                                        singularity::thermalqs::specific_heat |
                                        singularity::thermalqs::bulk_modulus);
                lambda[0] = Ye;

                sc_d.FillEos(R, T, e1, p1, cv1, b1, output, lambda);
                ig_d.FillEos(R, T, e2, p2, cv2, b2, output, lambda);
                if (!isClose(e1, e2)) {
                  nwrong_d() += 1;
                }
                if (!isClose(p1, p2)) {
                  nwrong_d() += 1;
                }
                if (!isClose(cv1, cv2)) {
                  nwrong_d() += 1;
                }
                if (!isClose(b1, b2)) {
                  nwrong_d() += 1;
                }
              });
#ifdef PORTABILITY_STRATEGY_KOKKOS
          Kokkos::deep_copy(nwrong_h, nwrong_d);
#endif
          REQUIRE(nwrong_h == 0);
          sc_d.Finalize();
        }
        ig_d.Finalize();
        ig.Finalize();
      }
      AND_THEN("We can save the file to SP5") {
        sc.Save(savename);
        AND_THEN("We can load the sp5 file") {
          StellarCollapse sc2(savename, true);
          AND_THEN("The two stellar collapse EOS's agree") {

            Real yemin = sc.YeMin();
            Real yemax = sc.YeMax();
            Real tmin = sc.TMin();
            Real tmax = sc.TMax();
            Real ltmin = std::log10(tmin);
            Real ltmax = std::log10(tmax);
            Real lrhomin = sc.lRhoMin();
            Real lrhomax = sc.lRhoMax();
            REQUIRE(yemin == sc2.YeMin());
            REQUIRE(yemax == sc2.YeMax());
            REQUIRE(sc.lTMin() == sc2.lTMin());
            REQUIRE(sc.lTMax() == sc2.lTMax());
            REQUIRE(lrhomin == sc2.lRhoMin());
            REQUIRE(lrhomax == sc2.lRhoMax());

            auto sc1_d = sc.GetOnDevice();
            auto sc2_d = sc2.GetOnDevice();

            int nwrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
            using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
            Kokkos::View<int, atomic_view> nwrong_d("wrong");
#else
            PortableMDArray<int> nwrong_d(&nwrong_h, 1);
#endif

            const int N = 123;
            const Real dY = (yemax - yemin) / (N + 1);
            const Real dlT = (ltmax - ltmin) / (N + 1);
            const Real dlR = (lrhomax - lrhomin) / (N + 1);
            portableFor(
                "fill eos", 0, N, 0, N, 0, N,
                PORTABLE_LAMBDA(const int &k, const int &j, const int &i) {
                  Real lambda[2];
                  Real Ye = yemin + k * dY;
                  Real lT = ltmin + j * dlT;
                  Real lR = lrhomin + i * dlR;
                  Real T = std::pow(10., lT);
                  Real R = std::pow(10., lR);
                  Real e1, e2, p1, p2, cv1, cv2, b1, b2;
                  unsigned long output =
                      (singularity::thermalqs::pressure |
                       singularity::thermalqs::specific_internal_energy |
                       singularity::thermalqs::specific_heat |
                       singularity::thermalqs::bulk_modulus);
                  lambda[0] = Ye;

                  sc1_d.FillEos(R, T, e1, p1, cv1, b1, output, lambda);
                  sc2_d.FillEos(R, T, e2, p2, cv2, b2, output, lambda);
                  if (!isClose(e1, e2)) nwrong_d() += 1;
                  if (!isClose(p1, p2)) nwrong_d() += 1;
                  if (!isClose(cv1, cv2)) nwrong_d() += 1;
                  if (!isClose(b1, b2)) nwrong_d() += 1;
                });
#ifdef PORTABILITY_STRATEGY_KOKKOS
            Kokkos::deep_copy(nwrong_h, nwrong_d);
#endif
            REQUIRE(nwrong_h == 0);

            sc1_d.Finalize();
            sc2_d.Finalize();
          }
          sc2.Finalize();
        }
      }
      sc.Finalize();
    }
  }
}
#endif // SINGULARITY_TEST_STELLAR_COLLAPSE
#endif // USE_HDF5

int main(int argc, char *argv[]) {

#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
#endif
  int result;
  { result = Catch::Session().run(argc, argv); }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif
  return result;
}
