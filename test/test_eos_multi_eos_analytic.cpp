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

SCENARIO("Test the MultiEOS object with reactants and products EOS",
         "[MultiEOS][DavisReactants][DavisProducts]") {
  using Catch::Matchers::ContainsSubstring;
  using namespace singularity;
  using namespace singularity::IndexableTypes;
  using singularity::IndexerUtils::VariadicIndexer;

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
    constexpr Real rho0 = 1.890;   // g/cm^3
    constexpr Real e0 = 0.;        // erg / g
    constexpr Real P0 = 0.;        // microbar
    constexpr Real T0 = 297;       // K
    constexpr Real A = 1.8 * sqrtGPa; // sqrt(microbar)
    constexpr Real B = 4.6;
    constexpr Real C = 0.34;
    constexpr Real G0 = 0.56;
    constexpr Real Z = 0.0;
    constexpr Real alpha = 0.4265;
    constexpr Real Cv_r = 0.001074 * MJ_per_kg; // erg / g / K
    auto davis_r_eos =
        DavisReactants(rho0, e0, P0, T0, A, B, C, G0, Z, alpha, Cv_r);

    // Davis Products EOS
    constexpr Real a = 0.798311;
    constexpr Real b = 0.58;
    constexpr Real k = 1.35;
    constexpr Real n = 2.66182;
    constexpr Real vc = 0.75419;              // cm^3 / g
    constexpr Real pc = 3.2 * GPa;            // microbar
    constexpr Real Cv_p = 0.001072 * MJ_per_kg; // erg / g / K
    constexpr Real Q = 4.115 * MJ_per_kg;
    auto davis_p_eos = DavisProducts(a, b, k, n, vc, pc, Cv_p).Modify<ShiftedEOS>(-Q);

    // Create the multiEOS object
    constexpr size_t num_eos = 2;
    using LambdaT = VariadicIndexer<MassFraction<0>, MassFraction<1>>;
    auto multi_eos = MultiEOS(davis_r_eos, davis_p_eos);

    // Lookup result tolerances
    constexpr Real lookup_tol = 1.0e-12;

    // Tolerance for derivatives is fairly loose since we are comparing P-T,
    // rho-T, and rho-sie derivatives
    constexpr Real deriv_tol = 5.0e-04;

    // There is a lot of noise in numerical finite differences that can be
    // amplified when using thermodynamic identities. This tolerance is fairly
    // loose, which motivates direct access to P-T derivatives in the future.
    // Those can be mass-averaged and then thermo identities can be used with
    // those to compute mixture quantities. Really this is the most accurate way
    // to compute mixture derivatives for PTE
    constexpr Real thermo_identity_tol = 0.9; // Essentially order of magntiude

    WHEN("A mass fraction cutoff less than 0 is specified") {
      [[maybe_unused]] constexpr Real mf_cutoff = -0.01;
      THEN("Constructing the MultiEOS object should throw an exception") {
        REQUIRE_MAYBE_THROWS(MultiEOS(mf_cutoff, davis_r_eos, davis_p_eos));
      }
    }

    WHEN("A mass fraction cutoff equal to 1 is specified") {
      [[maybe_unused]] constexpr Real mf_cutoff = 1.0;
      THEN("Constructing the MultiEOS object should throw an exception") {
        REQUIRE_MAYBE_THROWS(MultiEOS(mf_cutoff, davis_r_eos, davis_p_eos));
      }
    }

    WHEN("A mass fraction cutoff greater than 1 is specified") {
      [[maybe_unused]] constexpr Real mf_cutoff = 2.0;
      THEN("Constructing the MultiEOS object should throw an exception") {
        REQUIRE_MAYBE_THROWS(MultiEOS(mf_cutoff, davis_r_eos, davis_p_eos));
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
      auto multi_eos_largeMF = MultiEOS(mf_cutoff, davis_r_eos, davis_p_eos);

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
                  decltype(MultiEOS(DavisReactants{}, ShiftedEOS<DavisProducts>{}))>;
      EOS multi_eos_in_variant = MultiEOS(davis_r_eos, davis_p_eos);

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
          Real rho_FillEos;
          Real temp_FillEos;
          Real sie_FillEos;
          Real pres_FillEos;
          Real cv_FillEos;
          Real bmod_FillEos;

          THEN("P-T lookups are consistent") {
            constexpr unsigned long input = thermalqs::pressure | thermalqs::temperature;
            constexpr unsigned long output = ~input;
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

    WHEN("A pressure-temperature lookup is performed") {
      // Populate lambda with mass fractions (only)
      std::array<Real, num_eos> set_mass_fracs{};
      set_mass_fracs.fill(1.0 / num_eos);
      LambdaT lambda{set_mass_fracs};

      // Get the P-T state at a high P and T. Use the material densities and
      // energies so that we can ensure we're perturbing from the same PTE state
      constexpr Real P = 1e10;
      constexpr Real T = 5000;
      Real rho_scr;
      Real sie_scr;
      std::array<Real, num_eos> density_mat{};
      std::array<Real, num_eos> sie_mat{};
      multi_eos.GetStatesFromPressureTemperature(rho_scr, sie_scr, P, T, density_mat,
                                                 sie_mat, lambda);
      // Make const to make sure we're not screwing anything up
      const Real rho = rho_scr;
      const Real sie = sie_scr;

      // Don't continue test if this function isn't already working
      REQUIRE(rho > 0);

      INFO("rho: " << rho);
      INFO("sie: " << sie);

      AND_WHEN("The PTE state is perturbed in density-temperature space") {
        Real dedT_R;
        Real dedR_T;
        Real dPdT_R;
        Real dPdR_T;
        multi_eos.PerturbPTEStateRT(rho, sie, P, T, density_mat, sie_mat, lambda, dedT_R,
                                    dedR_T, dPdT_R, dPdR_T);
        INFO("dedT_R =  " << dedT_R);
        INFO("dedR_T =  " << dedR_T);
        INFO("dPdT_R =  " << dPdT_R);
        INFO("dPdR_T =  " << dPdR_T);
        const Real B_T = rho * dPdR_T;
        const Real B_S = B_T + robust::ratio(T * dPdT_R * dPdT_R, rho * dedT_R);
        const Real Cv = dedT_R;
        const Real dpde = robust::ratio(dPdT_R, dedT_R);

        REQUIRE(B_S > B_T);

        AND_WHEN("The PTE state is perturbed in pressure-temperature space") {
          Real dedT_P;
          Real dedP_T;
          Real dRdT_P;
          Real dRdP_T;
          multi_eos.PerturbPTEStatePT(rho, sie, P, T, density_mat, sie_mat, lambda,
                                      dedT_P, dedP_T, dRdT_P, dRdP_T);

          INFO("dedT_P = " << dedT_P);
          INFO("dedP_T = " << dedP_T);
          INFO("dRdT_P = " << dRdT_P);
          INFO("dRdP_T = " << dRdP_T);

          const Real Cp = dedT_P;
          const Real alpha = -dRdT_P / rho;
          const Real K_T = dRdP_T / rho;

          THEN("The thermodynamic inequalities should hold") {
            INFO("K_T = " << K_T << "  1.0 / B_S = " << 1.0 / B_S);
            INFO("Cp = " << Cp << "  Cv = " << Cv);
            CHECK(K_T >= 1.0 / B_S);
            CHECK(Cp >= Cv);
          }

          THEN("Inverses should appropriate be related") {
            CHECK_THAT(dPdR_T, Catch::Matchers::WithinRel(1.0 / dRdP_T, deriv_tol));
            CHECK_THAT(K_T, Catch::Matchers::WithinRel(1.0 / B_T, deriv_tol));
          }

          THEN("The cyclic relations should hold") {
            CHECK_THAT(dRdP_T * dPdT_R / dRdT_P,
                       Catch::Matchers::WithinRel(-1.0, deriv_tol));
            CHECK_THAT(dPdR_T * dRdT_P / dPdT_R,
                       Catch::Matchers::WithinRel(-1.0, deriv_tol));
          }

          INFO("Cp = dedT_P = " << Cp);
          INFO("alpha = -dRdT_P / rho = " << alpha);
          INFO("K_T = dRdP_T / rho = " << K_T);

          THEN("The isentropic compressiblity should obey the thermodynamic identity") {
            INFO("K_T - T * alpha * alpha / rho / Cp = " << K_T - T * alpha * alpha /
                                                                      rho / Cp);
            CHECK_THAT(1.0 / B_S,
                       Catch::Matchers::WithinRel(K_T - T * alpha * alpha / rho / Cp,
                                                  thermo_identity_tol));
          }

          THEN("The difference in heat capacities should be correct") {
            INFO("Cp - Cv = " << Cp - Cv);
            INFO("T * alpha * alpha / K_T / rho = " << T * alpha * alpha / K_T / rho);
            CHECK_THAT(Cp - Cv, Catch::Matchers::WithinRel(T * alpha * alpha / K_T / rho,
                                                           thermo_identity_tol));
          }
          THEN("The ratio of heat capacities and moduli should be equal") {
            INFO("B_S / B_T = " << B_S / B_T);
            INFO("Cp / Cv = " << Cp / Cv);
            CHECK_THAT(B_S / B_T,
                       Catch::Matchers::WithinRel(Cp / Cv, thermo_identity_tol));
          }
        }

        Real rho_FillEos;
        Real temp_FillEos;
        Real sie_FillEos;
        Real pres_FillEos;
        Real cv_FillEos;
        Real bmod_FillEos;

        AND_WHEN("Density-temperature lookups are used") {
          const unsigned long input = thermalqs::density | thermalqs::temperature;
          const unsigned long output = ~input;
          rho_FillEos = rho;
          temp_FillEos = T;
          multi_eos.FillEos(rho_FillEos, temp_FillEos, sie_FillEos, pres_FillEos,
                            cv_FillEos, bmod_FillEos, output, lambda);

          THEN("The bulk modulus is consistent with the thermodynamic formulas using "
               "finite differeces") {
            const Real bmod_calc =
                multi_eos.BulkModulusFromDensityTemperature(rho, T, lambda);

            CHECK_THAT(bmod_calc, Catch::Matchers::WithinRel(B_S, deriv_tol));

            AND_THEN("The result from FillEos is also consistent") {
              CHECK_THAT(bmod_FillEos, Catch::Matchers::WithinRel(B_S, deriv_tol));
            }
          }

          THEN("The dpde value is consistent with the thermodynamic formulas using "
               "finite differeces") {
            const Real Gruneisen =
                multi_eos.GruneisenParamFromDensityTemperature(rho, T, lambda);
            const Real dpde_calc = rho * Gruneisen;

            CHECK_THAT(dpde_calc, Catch::Matchers::WithinRel(dpde, deriv_tol));
          }

          THEN("The heat capacity is consistent with the thermodynamic formulas using "
               "finite differeces") {
            const Real Cv_calc =
                multi_eos.SpecificHeatFromDensityTemperature(rho, T, lambda);

            CHECK_THAT(Cv_calc, Catch::Matchers::WithinRel(Cv, deriv_tol));

            AND_THEN("The result from FillEos is also consistent") {
              CHECK_THAT(cv_FillEos, Catch::Matchers::WithinRel(Cv, deriv_tol));
            }
          }
        }

        AND_WHEN("Density-energy lookups are used") {
          const unsigned long input =
              thermalqs::density | thermalqs::specific_internal_energy;
          const unsigned long output = ~input;
          rho_FillEos = rho;
          sie_FillEos = sie;
          multi_eos.FillEos(rho_FillEos, temp_FillEos, sie_FillEos, pres_FillEos,
                            cv_FillEos, bmod_FillEos, output, lambda);

          THEN("The bulk modulus is consistent with the thermodynamic formulas using "
               "finite differeces") {
            const Real bmod_calc =
                multi_eos.BulkModulusFromDensityInternalEnergy(rho, sie, lambda);

            CHECK_THAT(bmod_calc, Catch::Matchers::WithinRel(B_S, deriv_tol));

            AND_THEN("The result from FillEos is also consistent") {
              CHECK_THAT(bmod_FillEos, Catch::Matchers::WithinRel(B_S, deriv_tol));
            }
          }

          THEN("The dpde value is consistent with the thermodynamic formulas using "
               "finite differeces") {
            const Real Gruneisen =
                multi_eos.GruneisenParamFromDensityInternalEnergy(rho, sie, lambda);
            const Real dpde_calc = rho * Gruneisen;

            CHECK_THAT(dpde_calc, Catch::Matchers::WithinRel(dpde, deriv_tol));
          }

          THEN("The heat capacity is consistent with the thermodynamic formulas using "
               "finite differeces") {
            const Real Cv_calc =
                multi_eos.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);

            CHECK_THAT(Cv_calc, Catch::Matchers::WithinRel(Cv, deriv_tol));

            AND_THEN("The result from FillEos is also consistent") {
              CHECK_THAT(cv_FillEos, Catch::Matchers::WithinRel(Cv, deriv_tol));
            }
          }
        }
      }
    }

    // TODO: Check mass fraction averaging for minimum energy and entropy. Those
    // aren't enabled for Davis right now

    WHEN("The nlambda() member function is called") {
      constexpr auto nlambda = multi_eos.nlambda();
      // The Davis EOS don't have lambda requirements
      STATIC_REQUIRE(nlambda == num_eos);
    }

    WHEN("The PreferredInput() member function is called") {
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
