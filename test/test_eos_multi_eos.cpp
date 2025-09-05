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

#include <array>
#include <cmath>
#include <tuple>
#include <type_traits>
#include <utility>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#endif

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>

#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/indexable_types.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/base/variadic_utils.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_multi_eos.hpp>
#include <test/eos_unit_test_helpers.hpp>

using singularity::DavisProducts;
using singularity::DavisReactants;
using singularity::make_MultiEOS;
using singularity::MassFracAverageFunctor;
using singularity::ShiftedEOS;
using singularity::Variant;
using singularity::VolumeFracHarmonicAverageFunctor;
using singularity::IndexableTypes::MassFraction;
using singularity::IndexerUtils::VariadicIndexer;

using Catch::Matchers::ContainsSubstring;

SCENARIO("Test the MultiEOS object with reactants and products EOS",
         "[MultiEOS][DavisReactants][DavisProducts]") {
  GIVEN("A pair of EOS for reactants and products and a resulting MultiEOS "
        "object") {
    // Unit conversions
    [[maybe_unused]] constexpr Real cm = 1.;
    [[maybe_unused]] constexpr Real us = 1e-06;
    [[maybe_unused]] constexpr Real Mbcc_per_g = 1e12;
    constexpr Real GPa = 1.0e10;
    constexpr Real sqrtGPa = 1.0e5; // sqrt isn't constexpr...
    constexpr Real MJ_per_kg = 1.0e10;

    // Davis Reactants EOS
    constexpr Real rho0_DP = 1.890;   // g/cm^3
    constexpr Real e0_DP = 0.;        // erg / g
    constexpr Real P0_DP = 0.;        // microbar
    constexpr Real T0_DP = 297;       // K
    constexpr Real A = 1.8 * sqrtGPa; // sqrt(microbar)
    constexpr Real B = 4.6;
    constexpr Real C = 0.34;
    constexpr Real G0 = 0.56;
    constexpr Real Z = 0.0;
    constexpr Real alpha = 0.4265;
    constexpr Real Cv_DP = 0.001074 * MJ_per_kg; // erg / g / K
    auto davis_r_eos =
        DavisReactants(rho0_DP, e0_DP, P0_DP, T0_DP, A, B, C, G0, Z, alpha, Cv_DP);

    // Davis Products EOS
    constexpr Real a = 0.798311;
    constexpr Real b = 0.58;
    constexpr Real k = 1.35;
    constexpr Real n = 2.66182;
    constexpr Real vc = 0.75419;              // cm^3 / g
    constexpr Real pc = 3.2 * GPa;            // microbar
    constexpr Real Cv = 0.001072 * MJ_per_kg; // erg / g / K
    constexpr Real Q = 4.115 * MJ_per_kg;
    auto davis_p_eos = DavisProducts(a, b, k, n, vc, pc, Cv).Modify<ShiftedEOS>(-Q);

    // Create the multiEOS object
    constexpr size_t num_eos = 2;
    using LambdaT = VariadicIndexer<MassFraction<0>, MassFraction<1>>;
    auto multi_eos = make_MultiEOS(davis_r_eos, davis_p_eos);

    // Lookup result tolerances
    constexpr Real lookup_tol = 1.0e-12;
    constexpr Real deriv_tol = lookup_tol * 1e3;

    WHEN("A mass fraction cutoff less than 0 is specified") {
      [[maybe_unused]] constexpr Real mf_cutoff = -0.01;
      THEN("Constructing the MultiEOS object should throw an exception") {
        REQUIRE_MAYBE_THROWS(make_MultiEOS(mf_cutoff, davis_r_eos, davis_p_eos));
      }
    }

    WHEN("A mass fraction cutoff equal to 1 is specified") {
      [[maybe_unused]] constexpr Real mf_cutoff = 1.0;
      THEN("Constructing the MultiEOS object should throw an exception") {
        REQUIRE_MAYBE_THROWS(make_MultiEOS(mf_cutoff, davis_r_eos, davis_p_eos));
      }
    }

    WHEN("A mass fraction cutoff greater than 1 is specified") {
      [[maybe_unused]] constexpr Real mf_cutoff = 2.0;
      THEN("Constructing the MultiEOS object should throw an exception") {
        REQUIRE_MAYBE_THROWS(make_MultiEOS(mf_cutoff, davis_r_eos, davis_p_eos));
      }
    }

    WHEN("A set of mass fractions is provided in a lambda") {
      std::array<Real, num_eos> set_mass_fracs{};
      set_mass_fracs.fill(1.0 / num_eos);
      LambdaT lambda{set_mass_fracs};

      THEN("The MultiEOS object can translate these from the lambda to an array") {
        std::array<Real, num_eos> mass_fracs{};
        multi_eos.PopulateMassFracArray(mass_fracs, lambda);
        for (size_t m = 0; m < num_eos; m++) {
          INFO("Mass fraction test: mass_fracs[" << m << "] = " << mass_fracs[m]);
          CHECK_THAT(set_mass_fracs[m],
                     Catch::Matchers::WithinRel(mass_fracs[m], 1.0e-14));
        }
      }
    }

    WHEN("A mass fraction cutoff is specified in the constructor") {
      constexpr Real mf_cutoff = 0.001;
      auto multi_eos_largeMF = make_MultiEOS(mf_cutoff, davis_r_eos, davis_p_eos);

      THEN("The 'SetUpPTE' member function will correctly determine which materials are "
           "participating in PTE") {
        constexpr auto small_mf = mf_cutoff / 10.;
        constexpr auto large_mf = (1.0 - small_mf) / (num_eos - 1);

        // Somewhat realistic but irrelvant values
        constexpr Real density_tot = 1.5;
        constexpr Real sie_tot = 1.0e10;

        // Arrays
        std::array<Real, num_eos> mass_fracs{};
        mass_fracs.fill(large_mf);
        std::array<Real, num_eos> vol_fracs{};   // zero initialized
        std::array<Real, num_eos> density_mat{}; // zero initialized
        std::array<Real, num_eos> sie_mat{};     // zero initialized
        std::array<size_t, num_eos> pte_mats{};  // zero initialized
        size_t npte = 0;
        Real vfrac_tot = 0;
        auto eos_arr = multi_eos_largeMF.CreateEOSArray();

        // Select first index to have a small mass fraction so that it won't
        // participate
        mass_fracs[0] = small_mf;

        multi_eos_largeMF.SetUpPTE(mass_fracs, vol_fracs, density_mat, sie_mat, pte_mats,
                                   npte, vfrac_tot, density_tot, sie_tot, eos_arr);

        INFO("mass_fracs[0] = " << mass_fracs[0]);
        INFO("Small Mass Fraction = " << small_mf);
        INFO("Mass Fraction Cutoff = " << mf_cutoff);
        CHECK(npte == num_eos - 1);
        CHECK(pte_mats[0] != 0); // First PTE material should not be index 0
      }
    }

    WHEN("The EosType function is called") {
      auto davis_r_eostype = davis_r_eos.EosType();
      auto davis_p_eostype = davis_r_eos.EosType();
      auto multi_eostype = multi_eos.EosType();

      THEN("The MultiEOS EosType should contain the component EOS names") {
        CHECK_THAT(multi_eostype, ContainsSubstring(davis_r_eostype));
        CHECK_THAT(multi_eostype, ContainsSubstring(davis_p_eostype));
        CHECK_THAT(multi_eostype, ContainsSubstring("MultiEOS"));
      }
    }

    WHEN("The EosPyType function is called") {
      auto davis_r_eospytype = davis_r_eos.EosPyType();
      auto davis_p_eospytype = davis_r_eos.EosPyType();
      auto multi_eospytype = multi_eos.EosPyType();

      THEN("The MultiEOS EosType should contain the component EOS names") {
        CHECK_THAT(multi_eospytype, ContainsSubstring(davis_r_eospytype));
        CHECK_THAT(multi_eospytype, ContainsSubstring(davis_p_eospytype));
        CHECK_THAT(multi_eospytype, ContainsSubstring("MultiEOS"));
      }
    }

    THEN("The MultiEOS object can be placed in an EOS Variant") {
      using EOS =
          Variant<DavisReactants, ShiftedEOS<DavisProducts>,
                  decltype(make_MultiEOS(DavisReactants{}, ShiftedEOS<DavisProducts>{}))>;
      EOS multi_eos_in_variant = make_MultiEOS(davis_r_eos, davis_p_eos);

      AND_WHEN("A pressure-temperature lookup is performed") {
        // Populate lambda with mass fractions (only)
        std::array<Real, num_eos> set_mass_fracs{};
        set_mass_fracs.fill(1.0 / num_eos);
        LambdaT lambda{set_mass_fracs};

        // A high pressure and temperature
        constexpr Real P = 8.0e10;
        constexpr Real T = 8000;
        Real rho;
        Real sie;
        multi_eos_in_variant.DensityEnergyFromPressureTemperature(P, T, lambda, rho, sie);
        INFO("Original lookup P: " << P << "   T: " << T);
        INFO("MultiEOS density: " << rho << "   sie: " << sie);
        CHECK(rho > 0.);
        CHECK(std::isnormal(rho));
        CHECK(std::isnormal(sie));

        THEN("The result is consistent with doing the same for the individual EOS") {

          // Create an array of EOS
          const auto eos_arr = multi_eos.CreateEOSArray();

          // Create an array of mass fractions for adding up the individual values
          std::array<Real, num_eos> mass_fracs{};
          multi_eos.PopulateMassFracArray(mass_fracs, lambda);

          // Sum the individual EOS contributions weighted by mass fractions
          Real spvol_bulk = 0.;
          Real sie_bulk = 0.;
          for (size_t m = 0; m < num_eos; m++) {
            Real rho_m;
            Real sie_m;
            eos_arr[m].DensityEnergyFromPressureTemperature(P, T, nullptr, rho_m, sie_m);
            INFO("EOS[" << m << "]:\n" << eos_arr[m].EosType());
            CHECK(rho_m > 0.);
            CHECK(sie_m != 0.);
            spvol_bulk += mass_fracs[m] / rho_m;
            sie_bulk += mass_fracs[m] * sie_m;
          }

          // Quick check that we can invert the specific volume
          REQUIRE(spvol_bulk > 0.);
          Real rho_bulk = 1.0 / spvol_bulk;

          // Make sure the weighted sums add up to the result from the MultiEOS
          CHECK_THAT(rho_bulk, Catch::Matchers::WithinRel(rho, lookup_tol));
          CHECK_THAT(sie_bulk, Catch::Matchers::WithinRel(sie, lookup_tol));
        }

        AND_THEN("Density-energy lookups using the MultiEOS object all yield results "
                 "consistent with the original conditions") {
          INFO("Lookup -- Density: " << rho << "  Energy: " << sie);
          auto T_RE =
              multi_eos_in_variant.TemperatureFromDensityInternalEnergy(rho, sie, lambda);
          CHECK_THAT(T_RE, Catch::Matchers::WithinRel(T, lookup_tol));
          auto P_RE =
              multi_eos_in_variant.PressureFromDensityInternalEnergy(rho, sie, lambda);
          CHECK_THAT(P_RE, Catch::Matchers::WithinRel(P, lookup_tol));
        }

        AND_THEN("Density-temperature lookups using the MultiEOS object all yield "
                 "results consistent with the original conditions") {
          INFO("Lookup -- Density: " << rho << "  Temperature: " << T);
          auto P_RT = multi_eos_in_variant.PressureFromDensityTemperature(rho, T, lambda);
          CHECK_THAT(P_RT, Catch::Matchers::WithinRel(P, lookup_tol));
          auto E_RT =
              multi_eos_in_variant.InternalEnergyFromDensityTemperature(rho, T, lambda);
          CHECK_THAT(E_RT, Catch::Matchers::WithinRel(sie, lookup_tol));
        }

        AND_WHEN("FillEos() is called") {
          using namespace singularity;
          Real rho_FillEos;
          Real temp_FillEos;
          Real sie_FillEos;
          Real pres_FillEos;
          Real cv_FillEos;
          Real bmod_FillEos;

          THEN("P-T lookups are consistent") {
            const unsigned long input = thermalqs::pressure | thermalqs::temperature;
            const unsigned long output = ~input;
            temp_FillEos = T;
            pres_FillEos = P;
            multi_eos_in_variant.FillEos(rho_FillEos, temp_FillEos, sie_FillEos,
                                         pres_FillEos, cv_FillEos, bmod_FillEos, output,
                                         lambda);
            CHECK_THAT(rho_FillEos, Catch::Matchers::WithinRel(rho, lookup_tol));
            CHECK_THAT(sie_FillEos, Catch::Matchers::WithinRel(sie, lookup_tol));
          }

          THEN("rho-T lookups are consistent") {
            const unsigned long input = thermalqs::density | thermalqs::temperature;
            const unsigned long output = ~input;
            rho_FillEos = rho;
            temp_FillEos = T;
            multi_eos_in_variant.FillEos(rho_FillEos, temp_FillEos, sie_FillEos,
                                         pres_FillEos, cv_FillEos, bmod_FillEos, output,
                                         lambda);
            CHECK_THAT(pres_FillEos, Catch::Matchers::WithinRel(P, lookup_tol));
            CHECK_THAT(sie_FillEos, Catch::Matchers::WithinRel(sie, lookup_tol));
          }

          THEN("rho-sie lookups are consistent") {
            const unsigned long input =
                thermalqs::density | thermalqs::specific_internal_energy;
            const unsigned long output = ~input;
            rho_FillEos = rho;
            sie_FillEos = sie;
            multi_eos_in_variant.FillEos(rho_FillEos, temp_FillEos, sie_FillEos,
                                         pres_FillEos, cv_FillEos, bmod_FillEos, output,
                                         lambda);
            CHECK_THAT(pres_FillEos, Catch::Matchers::WithinRel(P, lookup_tol));
            CHECK_THAT(temp_FillEos, Catch::Matchers::WithinRel(T, lookup_tol));
          }

          THEN("rho-P lookups are consistent") {
            const unsigned long input = thermalqs::density | thermalqs::pressure;
            const unsigned long output = ~input;
            rho_FillEos = rho;
            pres_FillEos = P;
            multi_eos_in_variant.FillEos(rho_FillEos, temp_FillEos, sie_FillEos,
                                         pres_FillEos, cv_FillEos, bmod_FillEos, output,
                                         lambda);
            CHECK_THAT(temp_FillEos, Catch::Matchers::WithinRel(T, lookup_tol));
            CHECK_THAT(sie_FillEos, Catch::Matchers::WithinRel(sie, lookup_tol));
          }
        }
      }
    }

    WHEN("Volume-fraction harmonic averages are used for the bulk modulus and Gruneisen "
         "averaging") {
      using BAvgF = VolumeFracHarmonicAverageFunctor;
      using GAvgF = VolumeFracHarmonicAverageFunctor;
      auto models = std::make_tuple(davis_r_eos, davis_p_eos);
      auto eos = make_MultiEOS<BAvgF, GAvgF>(davis_r_eos, davis_p_eos);

      // Note that this constructor is not allowed
      // auto eos = make_MultiEOS<BAvgF, GAvgF>(models);

      AND_WHEN("A pressure-temperature lookup is performed") {
        using namespace singularity;

        // Populate lambda with mass fractions (only)
        std::array<Real, num_eos> set_mass_fracs{};
        set_mass_fracs.fill(1.0 / num_eos);
        LambdaT lambda{set_mass_fracs};

        // A high pressure and temperature
        constexpr Real P = 1e10;
        constexpr Real T = 5000;

        // Output values
        Real rho;
        Real sie;
        std::array<Real, num_eos> density_mat{};
        std::array<Real, num_eos> sie_mat{};

        // Get the material values
        eos.GetStatesFromPressureTemperature(rho, sie, P, T, density_mat, sie_mat,
                                             lambda);

        // Don't continue test if this function isn't already working
        REQUIRE(rho > 0);

        // The lookup functors for averaging
        auto bmod_func = [](auto const &model, Real const rho, Real const sie,
                            Real const temp, auto const lambda) {
          if (model.PreferredInput() == (thermalqs::density | thermalqs::temperature)) {
            return model.BulkModulusFromDensityTemperature(rho, temp);
          } else {
            return model.BulkModulusFromDensityInternalEnergy(rho, sie);
          }
        };

        auto gruneisen_func = [](auto const &model, Real const rho, Real const sie,
                                 Real const temp, auto const lambda) {
          if (model.PreferredInput() == (thermalqs::density | thermalqs::temperature)) {
            return model.GruneisenParamFromDensityTemperature(rho, temp);
          } else {
            return model.GruneisenParamFromDensityInternalEnergy(rho, sie);
          }
        };

        auto cv_func = [](auto const &model, Real const rho, Real const sie,
                          Real const temp, auto const lambda) {
          if (model.PreferredInput() == (thermalqs::density | thermalqs::temperature)) {
            return model.SpecificHeatFromDensityTemperature(rho, temp);
          } else {
            return model.SpecificHeatFromDensityInternalEnergy(rho, sie);
          }
        };

        // auto entropy_func = [](auto const &model, Real const rho, Real const sie,
        //                        Real const temp, auto const lambda) {
        //   if (model.PreferredInput() == (thermalqs::density | thermalqs::temperature))
        //   {
        //     return model.EntropyFromDensityTemperature(rho, temp);
        //   } else {
        //     return model.EntropyFromDensityInternalEnergy(rho, sie);
        //   }
        // };

        // auto min_sie_func = [](auto const &model, Real const rho, Real const sie,
        //                        Real const temp, auto const lambda) {
        //   return model.MinInternalEnergyFromDensity(rho);
        // };

        auto zbar_func = [](auto const &model, Real const rho, Real const sie,
                            Real const temp, auto const lambda) {
          if (model.PreferredInput() == (thermalqs::density | thermalqs::temperature)) {
            return model.SpecificHeatFromDensityTemperature(rho, temp);
          } else {
            return model.SpecificHeatFromDensityInternalEnergy(rho, sie);
          }
        };

        auto b_avg_f = BAvgF{};
        Real const bmod_avg = b_avg_f(bmod_func, models, density_mat, sie_mat, rho, T,
                                      lambda, std::make_index_sequence<num_eos>{});

        auto g_avg_f = GAvgF{};
        Real const gruneisen_avg =
            g_avg_f(gruneisen_func, models, density_mat, sie_mat, rho, T, lambda,
                    std::make_index_sequence<num_eos>{});

        auto cv_avg_f = MassFracAverageFunctor{};
        Real const cv_avg = cv_avg_f(cv_func, models, density_mat, sie_mat, rho, T,
                                     lambda, std::make_index_sequence<num_eos>{});

        // TODO: Enable these with Davis reactants/products to make sure the answer is
        // correct. At the moment, these have been run, but produce an error since the
        // member functions aren't enabled

        // auto entropy_avg_f = MassFracAverageFunctor{};
        // Real const entropy_avg = entropy_avg_f(entropy_func, models, density_mat,
        // sie_mat,
        //                                        rho, T, lambda,
        //                                        std::make_index_sequence<num_eos>{});

        // auto min_sie_avg_f = MassFracAverageFunctor{};
        // Real const min_sie_avg = min_sie_avg_f(min_sie_func, models, density_mat,
        // sie_mat,
        //                                        rho, T, lambda,
        //                                        std::make_index_sequence<num_eos>{});

        THEN("The MultiEOS object will return a bulk modulus consistent with the "
             "averaging functor and density-energy lookups") {

          Real bmod_RE = eos.BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
          CHECK_THAT(bmod_RE, Catch::Matchers::WithinRel(bmod_avg, deriv_tol));

          AND_WHEN("FillEos is used") {
            Real rho_FillEos;
            Real temp_FillEos;
            Real sie_FillEos;
            Real pres_FillEos;
            Real cv_FillEos;
            Real bmod_FillEos;

            THEN("The density-energy result is still consistent") {
              const unsigned long input =
                  thermalqs::density | thermalqs::specific_internal_energy;
              const unsigned long output = ~input;
              rho_FillEos = rho;
              sie_FillEos = sie;
              eos.FillEos(rho_FillEos, temp_FillEos, sie_FillEos, pres_FillEos,
                          cv_FillEos, bmod_FillEos, output, lambda);
              CHECK_THAT(bmod_FillEos, Catch::Matchers::WithinRel(bmod_avg, deriv_tol));
            }
          }
        }

        THEN("The MultiEOS object will return a bulk modulus consistent with the "
             "averaging functor and density-temperature lookups") {
          Real bmod_RT = eos.BulkModulusFromDensityTemperature(rho, T, lambda);
          CHECK_THAT(bmod_RT, Catch::Matchers::WithinRel(bmod_avg, deriv_tol));

          AND_WHEN("FillEos is used") {
            Real rho_FillEos;
            Real temp_FillEos;
            Real sie_FillEos;
            Real pres_FillEos;
            Real cv_FillEos;
            Real bmod_FillEos;

            THEN("The density-temperature result is still consistent") {
              const unsigned long input = thermalqs::density | thermalqs::temperature;
              const unsigned long output = ~input;
              rho_FillEos = rho;
              temp_FillEos = T;
              eos.FillEos(rho_FillEos, temp_FillEos, sie_FillEos, pres_FillEos,
                          cv_FillEos, bmod_FillEos, output, lambda);
              CHECK_THAT(bmod_FillEos, Catch::Matchers::WithinRel(bmod_avg, deriv_tol));
            }
          }
        }

        THEN("The MultiEOS object will return a Gruneisen parameter consistent with the "
             "averaging functor and a density-energy lookup") {

          Real gruneisen_RE =
              eos.GruneisenParamFromDensityInternalEnergy(rho, sie, lambda);
          CHECK_THAT(gruneisen_RE, Catch::Matchers::WithinRel(gruneisen_avg, deriv_tol));
        }

        THEN("The MultiEOS object will return a Gruneisen parameter consistent with the "
             "averaging functor and a density-temperature lookup") {
          Real gruneisen_RT = eos.GruneisenParamFromDensityTemperature(rho, T, lambda);
          CHECK_THAT(gruneisen_RT, Catch::Matchers::WithinRel(gruneisen_avg, deriv_tol));
        }

        THEN("The MultiEOS object will return a mass-fraction-averaged specific heat "
             "capacity from density-energy lookups") {

          Real cv_RE = eos.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
          CHECK_THAT(cv_RE, Catch::Matchers::WithinRel(cv_avg, deriv_tol));

          AND_WHEN("FillEos is used") {
            Real rho_FillEos;
            Real temp_FillEos;
            Real sie_FillEos;
            Real pres_FillEos;
            Real cv_FillEos;
            Real bmod_FillEos;

            THEN("The density-energy result is still consistent") {
              const unsigned long input =
                  thermalqs::density | thermalqs::specific_internal_energy;
              const unsigned long output = ~input;
              rho_FillEos = rho;
              sie_FillEos = sie;
              eos.FillEos(rho_FillEos, temp_FillEos, sie_FillEos, pres_FillEos,
                          cv_FillEos, bmod_FillEos, output, lambda);
              CHECK_THAT(cv_FillEos, Catch::Matchers::WithinRel(cv_avg, deriv_tol));
            }
          }
        }

        THEN("The MultiEOS object will return a mass-fraction-averaged specific heat "
             "capacity from density-temperature lookups") {
          Real cv_RT = eos.SpecificHeatFromDensityTemperature(rho, T, lambda);
          CHECK_THAT(cv_RT, Catch::Matchers::WithinRel(cv_avg, deriv_tol));

          AND_WHEN("FillEos is used") {
            Real rho_FillEos;
            Real temp_FillEos;
            Real sie_FillEos;
            Real pres_FillEos;
            Real cv_FillEos;
            Real bmod_FillEos;

            THEN("The density-temperature result is still consistent") {
              const unsigned long input = thermalqs::density | thermalqs::temperature;
              const unsigned long output = ~input;
              rho_FillEos = rho;
              temp_FillEos = T;
              eos.FillEos(rho_FillEos, temp_FillEos, sie_FillEos, pres_FillEos,
                          cv_FillEos, bmod_FillEos, output, lambda);
              CHECK_THAT(cv_FillEos, Catch::Matchers::WithinRel(cv_avg, deriv_tol));
            }
          }
        }

        // THEN("The MultiEOS object will return a mass-fraction-averaged entropy from "
        //      "density-energy lookups") {

        //   Real entropy_RE = eos.EntropyFromDensityInternalEnergy(rho, sie, lambda);
        //   CHECK_THAT(entropy_RE, Catch::Matchers::WithinRel(entropy_avg, lookup_tol));
        // }

        // THEN("The MultiEOS object will return a mass-fraction-averaged entropy from "
        //      "density-temperature lookups") {
        //   Real entropy_RT = eos.EntropyFromDensityTemperature(rho, T, lambda);
        //   CHECK_THAT(entropy_RT, Catch::Matchers::WithinRel(entropy_avg, lookup_tol));
        // }

        // THEN("The MultiEOS object will return a mass-fraction-averaged minimum specific
        // "
        //      "internal energy from a density lookup") {

        //   Real min_sie_RE = eos.MinInternalEnergyFromDensity(rho, lambda);
        //   CHECK_THAT(min_sie_RE, Catch::Matchers::WithinRel(min_sie_avg, lookup_tol));
        // }
      }
    }

    WHEN("The nlambda() member function is called") {
      constexpr auto nlambda = multi_eos.nlambda();
      // The Davis EOS don't have lambda requirements
      STATIC_REQUIRE(nlambda == num_eos);
    }

    WHEN("The PreferredInput() member function is called") {
      using namespace singularity;
      constexpr auto pref_in = multi_eos.PreferredInput();
      // The preferred input is _always_ P-T. See comments in code member function
      STATIC_REQUIRE(pref_in == (thermalqs::pressure | thermalqs::temperature));
    }

    WHEN("The scratch_size() member function is called") {
      const auto scratch_sz = multi_eos.scratch_size("Any Method", 100);
      // Scratch is internally allocated
      REQUIRE(scratch_sz == 0);
    }

    WHEN("The PrintParams() member function is called") {
      multi_eos.PrintParams();
      THEN("Nothing goes wrong"){};
    }

    WHEN("The MinimumDensity() member function is called") {
      const auto min_rho = multi_eos.MinimumDensity();
      // Davis EOS minimum density is zero
      REQUIRE(min_rho == 0);
    }

    WHEN("The MinimumDensity() member function is called") {
      const Real eos_val = multi_eos.MinimumDensity();
      auto eos_arr = multi_eos.CreateEOSArray();
      Real calc = -1e300;
      for (size_t m = 0; m < num_eos; m++) {
        calc = std::max(calc, eos_arr[m].MinimumDensity());
      }
      REQUIRE_THAT(eos_val, Catch::Matchers::WithinRel(calc, lookup_tol));
    }

    WHEN("The MinimumTemperature() member function is called") {
      const auto eos_val = multi_eos.MinimumTemperature();
      auto eos_arr = multi_eos.CreateEOSArray();
      Real calc = -1e300;
      for (size_t m = 0; m < num_eos; m++) {
        calc = std::max(calc, eos_arr[m].MinimumTemperature());
      }
      REQUIRE_THAT(eos_val, Catch::Matchers::WithinRel(calc, lookup_tol));
    }

    WHEN("The MaximumDensity() member function is called") {
      const auto eos_val = multi_eos.MaximumDensity();
      auto eos_arr = multi_eos.CreateEOSArray();
      Real calc = 1e300;
      for (size_t m = 0; m < num_eos; m++) {
        calc = std::min(calc, eos_arr[m].MaximumDensity());
      }
      REQUIRE_THAT(eos_val, Catch::Matchers::WithinRel(calc, lookup_tol));
    }

    WHEN("The MinimumPressure() member function is called") {
      const auto eos_val = multi_eos.MinimumPressure();
      auto eos_arr = multi_eos.CreateEOSArray();
      Real calc = -1e300;
      for (size_t m = 0; m < num_eos; m++) {
        calc = std::max(calc, eos_arr[m].MinimumPressure());
      }
      REQUIRE_THAT(eos_val, Catch::Matchers::WithinRel(calc, lookup_tol));
    }

    WHEN("The MaximumPressureAtTemperature() member function is called") {
      constexpr Real T = 2000;
      const auto eos_val = multi_eos.MaximumPressureAtTemperature(T);
      auto eos_arr = multi_eos.CreateEOSArray();
      Real calc = 1e300;
      for (size_t m = 0; m < num_eos; m++) {
        calc = std::min(calc, eos_arr[m].MaximumPressureAtTemperature(T));
      }
      REQUIRE_THAT(eos_val, Catch::Matchers::WithinRel(calc, lookup_tol));
    }

    WHEN("The IsModified() member function is called") {
      const auto modified = multi_eos.IsModified();
      // By definition, EOS is not modified
      STATIC_REQUIRE(!modified);
    }
    WHEN("The UnmodifyOnce() member function is called") {
      using namespace ::singularity::variadic_utils;
      // By definition, EOS is not modified, and thus we should get the same
      // type back
      STATIC_REQUIRE(std::is_same<remove_cvref_t<decltype(multi_eos.UnmodifyOnce())>,
                                  remove_cvref_t<decltype(multi_eos)>>::value);
    }
  }
}

// TODO: Write a new test case using EOS with dynamic memory and probably also
// probably include 3T properties
