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

#ifdef SINGULARITY_TEST_SESAME
#ifdef SINGULARITY_USE_SPINER_WITH_HDF5

#include <array>
#include <cmath>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>

#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/indexable_types.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/base/variadic_utils.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_multi_eos.hpp>
#include <singularity-eos/eos/eos_spiner_rho_temp.hpp>
#include <test/eos_unit_test_helpers.hpp>

SCENARIO("Test the MultiEOS object with a binary alloy", "[MultiEOS][SpinerEOS]") {
  using namespace singularity;
  using singularity::IndexerUtils::VariadicIndexer;

  GIVEN("A pair of metal EOS for the binary alloy") {
    constexpr int num_eos = 2;
    constexpr int Al_matid = 3720;
    constexpr int Cu_matid = 3337;
    const std::string eos_file = "../materials.sp5";
    auto Cu_eos = SpinerEOSDependsRhoT{eos_file, Cu_matid};
    auto Al_eos = SpinerEOSDependsRhoT{eos_file, Al_matid};
    using LambdaT =
        VariadicIndexer<IndexableTypes::LogDensity, IndexableTypes::LogTemperature,
                        IndexableTypes::MassFraction<0>, IndexableTypes::MassFraction<1>>;
    constexpr size_t lambda_mf_offset = 2;

    // Lookup result tolerances
    constexpr Real lookup_tol = 1.0e-11;
    constexpr Real deriv_tol = lookup_tol * 1e3;
    constexpr Real dx_factor = 1.0e-07;
    constexpr Real finite_diff_tol = 1.0e-02;

    WHEN("The two EOS are combined in a MultiEOS object") {
      auto alloy = MultiEOS(Cu_eos, Al_eos);

      // Set equal mass fractions
      std::array<Real, num_eos> set_mass_fracs{};
      set_mass_fracs.fill(1.0 / num_eos);
      LambdaT lambda{};
      for (size_t i = 0; i < num_eos; i++) {
        lambda[lambda_mf_offset + i] = set_mass_fracs[i];
      }
      const LambdaT lambda_init{lambda};

      AND_WHEN("The alloy reference state values are returned") {
        Real rho_ref;
        Real temp_ref;
        Real sie_ref;
        Real pres_ref;
        Real cv_ref;
        Real bmod_ref;
        Real dpde_ref;
        Real dvdt_ref;
        alloy.ValuesAtReferenceState(rho_ref, temp_ref, sie_ref, pres_ref, cv_ref,
                                     bmod_ref, dpde_ref, dvdt_ref, lambda);

        // Sanity check to make sure the EOS are not overwriting mass fraction
        // values
        for (size_t i = 0; i < num_eos; i++) {
          INFO("i: " << i);
          REQUIRE_THAT(
              lambda[i + lambda_mf_offset],
              Catch::Matchers::WithinRel(lambda_init[i + lambda_mf_offset], lookup_tol));
        }

        INFO("pres_ref = " << pres_ref << "   temp_ref = " << temp_ref);
        INFO("rho_ref  = " << rho_ref << "  sie_ref = " << sie_ref);

        THEN("They agree with the states at the reference pressure and temperature") {
          std::array<Real, num_eos> density_mat;
          std::array<Real, num_eos> sie_mat;
          Real rho_PT;
          Real sie_PT;
          alloy.GetStatesFromPressureTemperature(rho_PT, sie_PT, pres_ref, temp_ref,
                                                 density_mat, sie_mat, lambda);

          // Sanity check to make sure the EOS are not overwriting mass fraction
          // values
          for (size_t i = 0; i < num_eos; i++) {
            INFO("i: " << i);
            REQUIRE_THAT(lambda[i + lambda_mf_offset],
                         Catch::Matchers::WithinRel(lambda_init[i + lambda_mf_offset],
                                                    lookup_tol));
          }

          CHECK_THAT(rho_ref, Catch::Matchers::WithinRel(rho_PT, lookup_tol));
          CHECK_THAT(sie_ref, Catch::Matchers::WithinRel(sie_PT, lookup_tol));

          AND_THEN("The lookups on the individual EOS are consistent") {
            const auto eos_arr = alloy.CreateEOSArray();
            for (size_t i = 0; i < num_eos; i++) {
              INFO("EOS[" << i << "]: " << eos_arr[i].EosType());
              Real rho_mat_PT;
              Real sie_mat_PT;
              eos_arr[i].DensityEnergyFromPressureTemperature(pres_ref, temp_ref, lambda,
                                                              rho_mat_PT, sie_mat_PT);

              // Sanity check to make sure the EOS are not overwriting mass fraction
              // values
              for (size_t j = 0; j < num_eos; j++) {
                INFO("j: " << j);
                REQUIRE_THAT(lambda[j + lambda_mf_offset],
                             Catch::Matchers::WithinRel(lambda_init[j + lambda_mf_offset],
                                                        lookup_tol));
              }

              CHECK_THAT(rho_mat_PT,
                         Catch::Matchers::WithinRel(density_mat[i], lookup_tol));
              CHECK_THAT(sie_mat_PT, Catch::Matchers::WithinRel(sie_mat[i], lookup_tol));
            }
          }
        }

        THEN("They agree with values from FillEos at the reference pressure and "
             "temperature") {
          constexpr unsigned long input = thermalqs::pressure | thermalqs::temperature;
          constexpr unsigned long output = ~input;
          Real pres_FillEos = pres_ref;
          Real temp_FillEos = temp_ref;
          Real rho_FillEos;
          Real sie_FillEos;
          Real cv_FillEos;
          Real bmod_FillEos;
          alloy.FillEos(rho_FillEos, temp_FillEos, sie_FillEos, pres_FillEos, cv_FillEos,
                        bmod_FillEos, output, lambda);

          // Sanity check to make sure the EOS are not overwriting mass fraction
          // values
          for (size_t i = 0; i < num_eos; i++) {
            INFO("i: " << i);
            REQUIRE_THAT(lambda[i + lambda_mf_offset],
                         Catch::Matchers::WithinRel(lambda_init[i + lambda_mf_offset],
                                                    lookup_tol));
          }

          CHECK_THAT(rho_ref, Catch::Matchers::WithinRel(rho_FillEos, lookup_tol));
          CHECK_THAT(sie_ref, Catch::Matchers::WithinRel(sie_FillEos, lookup_tol));
          CHECK_THAT(cv_ref, Catch::Matchers::WithinRel(cv_FillEos, deriv_tol));
          CHECK_THAT(bmod_ref, Catch::Matchers::WithinRel(bmod_FillEos, deriv_tol));
        }

        THEN("dpde agrees with the Gruneisen from density and temperature") {
          const Real Gruneisen =
              alloy.GruneisenParamFromDensityTemperature(rho_ref, temp_ref, lambda);
          const Real dpde_calc = rho_ref * Gruneisen;

          // We need a looser tolerance here because the Gruneisen parameter by itself
          // comes from density-energy lookups, which require a P-T lookup and root-find.
          const Real gru_tol = 1.0e-03;
          CHECK_THAT(dpde_ref, Catch::Matchers::WithinRel(dpde_calc, gru_tol));
        }

        AND_WHEN("A density-temperature finite difference approximation is used") {

          std::array<Real, num_eos> density_mat{};
          std::array<Real, num_eos> sie_mat{};
          SolverStatus status;

          const Real dT = temp_ref * dx_factor;
          const Real dR = rho_ref * dx_factor;
          const Real T_pert = temp_ref + dT;
          const Real R_pert = rho_ref + dR;
          INFO("T_pert = " << T_pert << "  R_pert = " << R_pert);
          INFO("  dT = " << dT);
          INFO("  dR = " << dR);

          // Tighter tolerance is required for PTE solvers to produce accurate
          // derivatives
          constexpr bool doing_derivs = true;

          // Recalculate the pressure so we can get some nice accurate
          // derivatives
          Real sie_ref_RT;
          Real pres_ref_RT;
          status = alloy.GetStatesFromDensityTemperature(rho_ref, sie_ref_RT, pres_ref_RT,
                                                         temp_ref, density_mat, sie_mat,
                                                         lambda, doing_derivs);

          // Convergence sanity check
          THEN("The density-temperature lookup at the reference state should converge") {
            INFO("SOLVER STATUS (reference state):");
            INFO("    slow_convergence_detected: " << status.slow_convergence_detected);
            INFO("    max_niter: " << status.max_niter);
            INFO("    max_line_niter: " << status.max_line_niter);
            INFO("    small_step_iters: " << status.small_step_iters);
            INFO("    residual: " << status.residual);
            REQUIRE(status.converged);
          }

          // Pressure sanity check
          THEN("Density-temperature lookups for the individual EOS should agree with the "
               "state returned by the reference density and temperature lookup") {
            const auto eos_arr = alloy.CreateEOSArray();
            for (size_t m = 0; m < num_eos; m++) {
              INFO("m: " << m);
              INFO("  rho[m] = " << density_mat[m]);
              REQUIRE_THAT(pres_ref_RT, Catch::Matchers::WithinRel(
                                            eos_arr[m].PressureFromDensityTemperature(
                                                density_mat[m], temp_ref, lambda),
                                            1.0e-9 // Slightly looser tolerance needed
                                            ));
            }
          }

          INFO("  P(rho_ref, T_ref) = " << pres_ref_RT);
          INFO("  e(rho_ref, T_ref) = " << sie_ref_RT);

          // Sanity check: let's make sure the pressures for the perturbed states are
          // reasonable. We'll also save them so we can look at them in the debugger if
          // needed
          std::array<Real, num_eos> pmat_higherT{};
          std::array<Real, num_eos> pmat_higherR{};
          const auto eos_arr = alloy.CreateEOSArray();
          for (size_t m = 0; m < num_eos; m++) {
            pmat_higherT[m] = eos_arr[m].PressureFromDensityTemperature(
                density_mat[m], temp_ref + dT, lambda);
            pmat_higherR[m] = eos_arr[m].PressureFromDensityTemperature(
                density_mat[m] + dR, temp_ref, lambda);
          }
          THEN("Density-temperature lookups for the individual EOS at the perturbed "
               "temperature should yield slightly higher pressures") {
            for (size_t m = 0; m < num_eos; m++) {
              INFO("  m: " << m);
              INFO("    rho[m] = " << density_mat[m]);
              INFO("    P(rho[m], T+dT) = " << pmat_higherT[m]);
              CHECK(pmat_higherT[m] > pres_ref_RT);
            }
          }
          THEN("Perturbing the densities of the individual EOS should yield slightly "
               "higher pressures") {
            for (size_t m = 0; m < num_eos; m++) {
              INFO("  m: " << m);
              INFO("    rho[m] = " << density_mat[m]);
              INFO("    P(rho[m]+dR, T) = " << pmat_higherR[m]);
              CHECK(pmat_higherR[m] > pres_ref_RT);
            }
          }

          // Perturb pressure in temperature
          Real sie_dT;
          Real pres_dT;
          status = alloy.GetStatesFromDensityTemperature(rho_ref, sie_dT, pres_dT, T_pert,
                                                         density_mat, sie_mat, lambda,
                                                         doing_derivs);
          THEN("The density-temperature lookup at a perturbed temperature should "
               "converge") {
            INFO("SOLVER STATUS (perturbed in temperature):");
            INFO("    slow_convergence_detected: " << status.slow_convergence_detected);
            INFO("    max_niter: " << status.max_niter);
            INFO("    max_line_niter: " << status.max_line_niter);
            INFO("    small_step_iters: " << status.small_step_iters);
            INFO("    residual: " << status.residual);
            REQUIRE(status.converged);
          }
          INFO("  P(rho_ref, T_pert) = " << pres_dT);
          INFO("  e(rho_ref, T_pert) = " << sie_dT);
          INFO("    rho[0] @ (rho_ref, T_pert) = " << density_mat[0]);
          INFO("    rho[1] @ (rho_ref, T_pert) = " << density_mat[1]);

          // Perturb pressure in density
          Real sie_dR;
          Real pres_dR;
          status = alloy.GetStatesFromDensityTemperature(R_pert, sie_dR, pres_dR,
                                                         temp_ref, density_mat, sie_mat,
                                                         lambda, doing_derivs);
          THEN("The density-temperature lookup at a perturbed density should converge") {
            INFO("SOLVER STATUS (perturbed in density):");
            INFO("    slow_convergence_detected: " << status.slow_convergence_detected);
            INFO("    max_niter: " << status.max_niter);
            INFO("    max_line_niter: " << status.max_line_niter);
            INFO("    small_step_iters: " << status.small_step_iters);
            INFO("    residual: " << status.residual);
            REQUIRE(status.converged);
          }

          // Sanity check to make sure the EOS are not overwriting mass fraction
          // values
          for (size_t i = 0; i < num_eos; i++) {
            INFO("i: " << i);
            REQUIRE_THAT(lambda[i + lambda_mf_offset],
                         Catch::Matchers::WithinRel(lambda_init[i + lambda_mf_offset],
                                                    lookup_tol));
          }

          THEN("The reference dpde and dvdt values should agree with the finite "
               "difference") {

            INFO("  P(R_pert, temp_ref) = " << pres_dR);
            INFO("  e(R_pert, temp_ref) = " << sie_dR);
            INFO("    rho[0] @ (R_pert, temp_ref) = " << density_mat[0]);
            INFO("    rho[1] @ (R_pert, temp_ref) = " << density_mat[1]);

            const Real dP_R = pres_dT - pres_ref_RT;
            const Real dP_T = pres_dR - pres_ref_RT;
            const Real de_R = sie_dT - sie_ref_RT;
            const Real de_T = sie_dR - sie_ref_RT;

            const Real dPdT_R = dP_R / dT;
            const Real dPdR_T = dP_T / dR;
            const Real dVdP_T = robust::ratio(-1.0, rho_ref * rho_ref * dPdR_T);
            const Real dVdT_FD = -dPdT_R * dVdP_T;

            INFO("  de(dT) @ const. R = " << de_R);
            INFO("  de(dR) @ const. T = " << de_T);
            INFO("  dP(dT) @ const. R = " << dP_R);
            INFO("  dP(dR) @ const. T = " << dP_T);

            INFO("dPdT_R = " << dPdT_R);
            INFO("dPdR_T = " << dPdR_T);
            INFO("dVdP_T = " << dVdP_T);
            CHECK_THAT(dvdt_ref, Catch::Matchers::WithinRel(dVdT_FD, finite_diff_tol));

            const Real dpde_FD = robust::ratio(dP_R, de_R);
            CHECK_THAT(dpde_ref, Catch::Matchers::WithinRel(dpde_FD, finite_diff_tol));
          }
        }
      }

      AND_WHEN("We use the MinimumDensity() function") {
        auto min_rho = alloy.MinimumDensity();
        THEN("The MultiEOS minimum is the maximum of the individual minimums") {
          const auto eos_arr = alloy.CreateEOSArray();
          Real min_value = std::numeric_limits<Real>::lowest();
          for (auto eos : eos_arr) {
            min_value = std::max(min_value, eos.MinimumDensity());
          }
          CHECK_THAT(min_value, Catch::Matchers::WithinRel(min_rho, 1.0e-14));
        }
      }
      AND_WHEN("We use the MinimumTemperature() function") {
        auto min_T = alloy.MinimumTemperature();
        THEN("The MultiEOS minimum is the maximum of the individual minimums") {
          const auto eos_arr = alloy.CreateEOSArray();
          Real min_value = std::numeric_limits<Real>::lowest();
          for (auto eos : eos_arr) {
            min_value = std::max(min_value, eos.MinimumTemperature());
          }
          CHECK_THAT(min_value, Catch::Matchers::WithinRel(min_T, 1.0e-14));
        }
      }
      AND_WHEN("We use the MaximumDensity() function") {
        auto max_rho = alloy.MaximumDensity();
        THEN("The MultiEOS maximum is the minimum of the individual maxmums") {
          const auto eos_arr = alloy.CreateEOSArray();
          Real max_value = std::numeric_limits<Real>::max();
          for (auto eos : eos_arr) {
            max_value = std::min(max_value, eos.MaximumDensity());
          }
          CHECK_THAT(max_value, Catch::Matchers::WithinRel(max_rho, 1.0e-14));
        }
      }
      AND_WHEN("We use the MinimumPressure() function") {
        auto min_P = alloy.MinimumPressure();
        THEN("The MultiEOS minimum is the maximum of the individual minimums") {
          const auto eos_arr = alloy.CreateEOSArray();
          Real min_value = std::numeric_limits<Real>::lowest();
          for (auto eos : eos_arr) {
            min_value = std::max(min_value, eos.MinimumPressure());
          }
          CHECK_THAT(min_value, Catch::Matchers::WithinRel(min_P, 1.0e-14));
        }
      }
      AND_WHEN("We use the MaximumPressureAtTemperature() function") {
        // Initialize some temperatures
        auto Ts = std::array<Real, 10>{};
        for (size_t i = 0; i < Ts.size(); i++) {
          Ts[i] = (i + 1) * 100;
        }

        THEN("The MultiEOS maximum is the minimum of the individual maximums") {
          const auto eos_arr = alloy.CreateEOSArray();
          for (const auto T : Ts) {
            const auto max_P = alloy.MaximumPressureAtTemperature(T);
            Real max_value = std::numeric_limits<Real>::max();
            for (const auto eos : eos_arr) {
              max_value = std::min(max_value, eos.MaximumPressureAtTemperature(T));
            }
            INFO("T = " << T);
            CHECK_THAT(max_value, Catch::Matchers::WithinRel(max_P, 1.0e-14));
          }
        }
      }

      MassFracAverageFunctor mass_avg_funct{};
      AND_WHEN("We use the MeanAtomicMassFromDensityTemperature() function") {
        std::array<Real, num_eos> mfracs{};
        mfracs.fill(1.0 / num_eos);
        LambdaT lambda_equal_mass_fracs{};
        constexpr size_t lambda_mf_offset = 2;
        for (size_t i = 0; i < num_eos; i++) {
          lambda_equal_mass_fracs[lambda_mf_offset + i] = mfracs[i];
        }
        // For a standard spiner EOS, there is no density or temperature dependence
        auto Ts = std::array<Real, 10>{};
        auto Ds = std::array<Real, 10>{};
        for (size_t i = 0; i < Ts.size(); i++) {
          Ts[i] = (i * i + 1) * 100;
        }
        for (size_t i = 0; i < Ds.size(); i++) {
          Ds[i] = (i * i * i * i + 1) * 0.01;
        }

        THEN("The returned value is the mass-fraction average of the individual values") {
          const auto eos_arr = alloy.CreateEOSArray();
          for (auto T : Ts) {
            for (auto rho : Ds) {
              Real mean_value = 0;
              for (size_t i = 0; i < num_eos; i++) {
                mean_value +=
                    set_mass_fracs[i] * eos_arr[i].MeanAtomicMassFromDensityTemperature(
                                            rho, T, lambda_equal_mass_fracs);
              }
              INFO("rho = " << rho << "  T = " << T);
              Real abar_RT = alloy.MeanAtomicMassFromDensityTemperature(
                  rho, T, lambda_equal_mass_fracs);
              CHECK_THAT(mean_value, Catch::Matchers::WithinRel(abar_RT, lookup_tol));
              // Note that this will only mass when the mass fractions are equal
              Real abar = alloy.MeanAtomicMass();
              CHECK_THAT(abar_RT, Catch::Matchers::WithinRel(abar, 1.0e-14));
            }
          }
        }
      }
      AND_WHEN("We use the MeanAtomicNumberFromDensityTemperature() function") {
        std::array<Real, num_eos> mfracs{};
        mfracs.fill(1.0 / num_eos);
        LambdaT lambda_equal_mass_fracs{};
        constexpr size_t lambda_mf_offset = 2;
        for (size_t i = 0; i < num_eos; i++) {
          lambda_equal_mass_fracs[lambda_mf_offset + i] = mfracs[i];
        }
        // For a standard spiner EOS, there is no density or temperature dependence
        auto Ts = std::array<Real, 10>{};
        auto Ds = std::array<Real, 10>{};
        for (size_t i = 0; i < Ts.size(); i++) {
          Ts[i] = (i * i + 1) * 100;
        }
        for (size_t i = 0; i < Ds.size(); i++) {
          Ds[i] = (i * i * i * i + 1) * 0.01;
        }

        THEN("The returned value is the mass-fraction average of the individual values") {
          const auto eos_arr = alloy.CreateEOSArray();
          for (auto T : Ts) {
            for (auto rho : Ds) {
              Real mean_value = 0;
              for (size_t i = 0; i < num_eos; i++) {
                mean_value +=
                    set_mass_fracs[i] * eos_arr[i].MeanAtomicNumberFromDensityTemperature(
                                            rho, T, lambda_equal_mass_fracs);
              }
              INFO("rho = " << rho << "  T = " << T);
              Real zbar_RT = alloy.MeanAtomicNumberFromDensityTemperature(
                  rho, T, lambda_equal_mass_fracs);
              CHECK_THAT(mean_value, Catch::Matchers::WithinRel(zbar_RT, lookup_tol));
              // Note that this will only mass when the mass fractions are equal
              Real zbar = alloy.MeanAtomicNumber();
              CHECK_THAT(zbar_RT, Catch::Matchers::WithinRel(zbar, 1.0e-14));
            }
          }
        }
      }

      AND_WHEN("We use the IsModified() function") {
        THEN("The EOS isn't modified") { REQUIRE(!alloy.IsModified()); }
        THEN("UnmodifyOnce returns the same EOS") {
          auto alloy_copy = alloy.UnmodifyOnce();
          REQUIRE(alloy_copy.EosType() == alloy.EosType());
        }
      }

      AND_WHEN("We want to perform lookups on device") {
        auto alloy_device = alloy.GetOnDevice();
        int nwrong = 0;
        portableReduce(
            "Test MultiEOS on device", 0, 1,
            PORTABLE_LAMBDA(int dummy_idx, int &nwrong) {
              // Initialize mass fractions
              std::array<Real, num_eos> mfracs{};
              mfracs.fill(1.0 / num_eos);
              LambdaT lambda_device{};
              for (size_t i = 0; i < num_eos; i++) {
                lambda_device[lambda_mf_offset + i] = set_mass_fracs[i];
              }

              const Real rho = 6.0;
              const Real T = 500;
              nwrong += !CheckRhoSieFromPT(alloy_device, rho, T, lambda_device);
            },
            nwrong);
        REQUIRE(nwrong == 0);
        alloy_device.Finalize();
      }
      alloy.Finalize();
    }
  }
}

SCENARIO("Test the MultiEOS object dynamic memory features",
         "[MultiEOS][SpinerEOS][EOSPAC]") {
  using namespace singularity;
  using singularity::IndexerUtils::VariadicIndexer;

  GIVEN("A pair of metal EOS for the binary alloy") {
    constexpr int num_eos = 2;
    constexpr int Al_matid = 3720;
    constexpr int Cu_matid = 3337;
    const std::string eos_file = "../materials.sp5";
    auto Cu_eos = SpinerEOSDependsRhoT{eos_file, Cu_matid};
    using CuEOS_t = decltype(Cu_eos);
    auto Al_eos = EOSPAC{Al_matid};
    using AlEOS_t = decltype(Al_eos);

    using LambdaT =
        VariadicIndexer<IndexableTypes::LogDensity, IndexableTypes::LogTemperature,
                        IndexableTypes::MassFraction<0>, IndexableTypes::MassFraction<1>>;

    WHEN("The two EOS are combined in either a bare MultiEOS object or one in a "
         "Variant") {
      auto alloy = MultiEOS(Cu_eos, Al_eos);
      using EOS = Variant<decltype(alloy)>;
      EOS alloy_in_variant = alloy;
      const auto Cu_serialized_size = Cu_eos.SerializedSizeInBytes();
      const auto Al_serialized_size = Al_eos.SerializedSizeInBytes();
      const auto Cu_dynamic_size = Cu_eos.DynamicMemorySizeInBytes();
      const auto Al_dynamic_size = Al_eos.DynamicMemorySizeInBytes();
      const auto Cu_shared_size = Cu_eos.SharedMemorySizeInBytes();
      const auto Al_shared_size = Al_eos.SharedMemorySizeInBytes();

      THEN("The MultiEOS serialized size should be greater than the sum of the "
           "individual models") {
        const auto serialized_size = alloy.SerializedSizeInBytes();
        INFO("MultiEOS serlialized size: " << serialized_size);
        INFO("Copper EOS serialized size: " << Cu_serialized_size);
        INFO("Aluminum EOS serialized size: " << Al_serialized_size);
        CHECK(serialized_size >= Al_serialized_size + Cu_serialized_size);
      }

      THEN("The MultiEOS dynamic size should just be the sum of the individual models") {
        const auto dynamic_size = alloy.DynamicMemorySizeInBytes();
        INFO("MultiEOS dynamic size: " << dynamic_size);
        INFO("Copper EOS dynamic size: " << Cu_dynamic_size);
        INFO("Aluminum EOS dynamic size: " << Al_dynamic_size);
        CHECK(dynamic_size == Al_dynamic_size + Cu_dynamic_size);
      }

      THEN("The MultiEOS shared size should just be the sum of the individual models") {
        const auto shared_size = alloy.SharedMemorySizeInBytes();
        INFO("MultiEOS shared size: " << shared_size);
        INFO("Copper EOS shared size: " << Cu_shared_size);
        INFO("Aluminum EOS shared size: " << Al_shared_size);
        CHECK(shared_size == Al_shared_size + Cu_shared_size);
      }

      // AND_WHEN("We serialize") {
      //   auto [bare_size, bare_data] = alloy.Serialize();
      //   auto [variant_size, variant_data] = alloy_in_variant.Serialize();

      //   THEN("The serialized sizes are correct") {
      //     REQUIRE(bare_size == alloy.SerializedSizeInBytes());
      //     REQUIRE(variant_size == alloy_in_variant.SerializedSizeInBytes());
      //     alloy.Finalize();
      //   }

      //   const std::size_t bare_dynamic_size = alloy.DynamicMemorySizeInBytes();
      //   THEN("The dynamic size is less than the total serialized size") {
      //     REQUIRE(bare_size > bare_dynamic_size);
      //     const std::size_t variant_dynamic_size =
      //         alloy_in_variant.DynamicMemorySizeInBytes();
      //     REQUIRE(bare_size > variant_dynamic_size);
      //     REQUIRE(bare_dynamic_size == variant_dynamic_size);
      //     alloy.Finalize();
      //   }

      //   THEN("We can deserialize into shared memory") {
      //     using singularity::SharedMemSettings;

      //     // Create location for shared memory
      //     std::size_t dynamic_size = alloy.DynamicMemorySizeInBytes();
      //     REQUIRE(dynamic_size <= bare_dynamic_size);
      //     char *shared_loc = (char *)malloc(dynamic_size);

      //     // Finalize step is for EOSPAC mostly
      //     alloy.Finalize();

      //     // Create empty objects and deserialize into them
      //     EOS alloy2 = MultiEOS<CuEOS_t, AlEOS_t>{};
      //     auto alloy3 = MultiEOS<CuEOS_t, AlEOS_t>{};
      //     std::size_t read_size_2 = alloy2.DeSerialize(
      //         variant_data, SharedMemSettings(shared_loc, /* copy into shared */
      //         true));
      //     REQUIRE(read_size_2 == bare_size);
      //     std::size_t read_size_3 = alloy3.DeSerialize(
      //         bare_data, SharedMemSettings(shared_loc, /* copy into shared */ false));
      //     REQUIRE(read_size_3 == bare_size);

      //     AND_THEN("EOS lookups work") {
      //       constexpr Real rho_trial = 5;
      //       constexpr Real sie_trial = 1e12;
      //       std::array<Real, num_eos> mfracs{};
      //       mfracs.fill(1.0 / num_eos);
      //       LambdaT lambda{};
      //       constexpr size_t lambda_mf_offset = 2;
      //       for (size_t i = 0; i < num_eos; i++) {
      //         lambda[lambda_mf_offset + i] = mfracs[i];
      //       }

      //       const Real P2 =
      //           alloy2.PressureFromDensityTemperature(rho_trial, sie_trial, lambda);
      //       const Real P3 =
      //           alloy3.PressureFromDensityTemperature(rho_trial, sie_trial, lambda);

      //       REQUIRE_THAT(P2, Catch::Matchers::WithinRel(P3, 1.0e-12));
      //     }

      //     alloy2.Finalize();
      //     alloy3.Finalize();
      //   }
      // }
    }
  }
}

#endif // SINGULARITY_USE_SPINER_WITH_HDF5
#endif // SINGULARITY_TEST_SESAME
