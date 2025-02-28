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
#include <vector>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/modifiers/floored_energy.hpp>
#include <test/eos_unit_test_helpers.hpp>

using singularity::DavisReactants;
using singularity::FlooredEnergy;
using singularity::Gruneisen;
using singularity::IdealGas;
using singularity::JWL;
using singularity::ShiftedEOS;
using singularity::robust::ratio;
#ifdef SINGULARITY_USE_SPINER_WITH_HDF5
using singularity::SpinerEOSDependsRhoT;
#endif

using EOS =
    singularity::Variant<FlooredEnergy<IdealGas>, FlooredEnergy<Gruneisen>,
                         FlooredEnergy<DavisReactants>, FlooredEnergy<ShiftedEOS<JWL>>
#ifdef SINGULARITY_USE_SPINER_WITH_HDF5
                         ,
                         FlooredEnergy<SpinerEOSDependsRhoT>
#endif
                         >;

// Helper function that returns the difference between a P(rho, T) and P(rho, e)
// lookup for a set of EOS. The density is at rho_ref * rho_factor while the
// temperature for the lookup is T_lookup. The energy for the lookup is given by
// e(rho, T_lookup) - e_offset
PORTABLE_INLINE_FUNCTION
auto diff_pressures(const int n_eos, EOS *v_EOS, const Real T_lookup,
                    const Real e_offset = 0., const Real rho_factor = 1.2) {
  // Create storage for relative diffs
  std::vector<Real> P_rdiffs(n_eos); // zero initialized
  const size_t bytes = n_eos * sizeof(Real);
  Real *v_P_rdiffs = (Real *)PORTABLE_MALLOC(bytes);
  portableCopyToDevice(v_P_rdiffs, P_rdiffs.data(), bytes);

  // Loop over EOS
  portableFor(
      "Positive temperature test", 0, n_eos, PORTABLE_LAMBDA(int i) {
        const auto this_eos = v_EOS[i];

        // Find the reference state (really we just want the density :shrug:)
        Real rho_ref;
        Real temp_ref;
        Real sie_ref;
        Real press_ref;
        Real cv_ref;
        Real bmod_ref;
        Real dpde_ref;
        Real dvdt_ref;
        this_eos.ValuesAtReferenceState(rho_ref, temp_ref, sie_ref, press_ref, cv_ref,
                                        bmod_ref, dpde_ref, dvdt_ref);

        const Real density_lookup = rho_factor * rho_ref;

        // Find energy at reference density and specified temperature
        const Real e_lookup =
            this_eos.InternalEnergyFromDensityTemperature(density_lookup, T_lookup) -
            e_offset;

        // Diff the P(rho, e) and P(rho, T) lookups
        const Real P_from_e =
            this_eos.PressureFromDensityInternalEnergy(density_lookup, e_lookup);
        const Real P_from_T =
            this_eos.PressureFromDensityTemperature(density_lookup, T_lookup);
        v_P_rdiffs[i] = ratio(P_from_e - P_from_T, (P_from_e + P_from_T) / 2.);
      });

  // Transfer to host
  portableCopyToHost(P_rdiffs.data(), v_P_rdiffs, bytes);

  // Free device memory
  PORTABLE_FREE(v_P_rdiffs);

  return P_rdiffs;
}

// Helper functor to get the name of the EOS
struct GetName {
  template <typename eosT>
  void operator()(const eosT &eos) {
    name = typeid(eos).name();
  }

  std::string name;
};

SCENARIO("Test the floored energy modifer for a suite of EOS",
         "[FlooredEnergy][IdealGas][Gruneisen][Davis][JWL]") {

  GIVEN("A set of EOS with the energy floor applied") {
    // Unit conversions
    constexpr Real cm = 1.;
    constexpr Real us = 1e-06;
    constexpr Real Mbcc_per_g = 1e12;
    constexpr Real GPa = 1.0e10;
    constexpr Real MJ_per_kg = 1.0e10;

    // Ideal gas air
    constexpr Real P0 = 1.0e6;       // 1 bar
    constexpr Real T0 = 296;         // K
    constexpr Real rho0_air = 1e-03; // g/cc
    constexpr Real Gruneisen_air = 0.4;
    constexpr Real Cv_air = P0 / rho0_air / (Gruneisen_air * T0);
    EOS air_eos = FlooredEnergy<IdealGas>(IdealGas(Gruneisen_air, Cv_air));

    // Davis Reactants EOS
    constexpr Real rho0_DP = 1.890;
    constexpr Real e0_DP = 0.;      // erg / g
    constexpr Real P0_DP = 0.;      // microbar
    constexpr Real T0_DP = 297;     // K
    constexpr Real A = 1.8 * 1.0e5; // 1.8 * sqrt(GPa) -> sqrt(microbar)
    constexpr Real B = 4.6;
    constexpr Real C = 0.34;
    constexpr Real G0 = 0.56;
    constexpr Real Z = 0.0;
    constexpr Real alpha = 0.4265;
    constexpr Real Cv_DP = 0.001074 * MJ_per_kg; // erg / g / K
    EOS davis_r_eos = FlooredEnergy<DavisReactants>(
        DavisReactants(rho0_DP, e0_DP, P0_DP, T0_DP, A, B, C, G0, Z, alpha, Cv_DP));

    // JWL EOS (Tarver & McGuire, 2002)
    constexpr Real A_JWL = 6.3207e4 * GPa;
    constexpr Real B_JWL = -4.472 * GPa;
    constexpr Real R1_JWL = 11.3;
    constexpr Real R2_JWL = 1.13;
    constexpr Real rho0_JWL = 1.895;
    constexpr Real w_JWL = 0.8938;
    constexpr Real Cv_JWL = 2.487e-3 / rho0_JWL * MJ_per_kg;
    constexpr Real E0_JWL = 0.246929 * MJ_per_kg;
    EOS jwl_eos = FlooredEnergy<ShiftedEOS<JWL>>(ShiftedEOS<JWL>(
        JWL(A_JWL, B_JWL, R1_JWL, R2_JWL, w_JWL, rho0_JWL, Cv_JWL), E0_JWL));

    // Gruneisen parameters for copper
    constexpr Real C0_G = 0.394 * cm / us;
    constexpr Real S1_G = 1.489;
    constexpr Real S2_G = 0.;
    constexpr Real S3_G = 0.;
    constexpr Real Gamma0_G = 2.02;
    constexpr Real b_G = 0.47;
    constexpr Real rho0_G = 8.93;
    constexpr Real Cv_G = 0.383e-05 * Mbcc_per_g;
    EOS gruneisen_eos = FlooredEnergy<Gruneisen>(
        Gruneisen(C0_G, S1_G, S2_G, S3_G, Gamma0_G, b_G, rho0_G, T0, P0, Cv_G));

    // Tabular EOS parameters (when used)
    constexpr int matid = 3337;
    const std::string eos_file = "../materials.sp5";
#ifdef SINGULARITY_TEST_SESAME
#ifdef SINGULARITY_USE_SPINER_WITH_HDF5
    EOS spiner_eos =
        FlooredEnergy<SpinerEOSDependsRhoT>(SpinerEOSDependsRhoT(eos_file, matid));
#endif
#endif

    // Put EOS in a vector and put EOS on device
    std::vector<EOS> eos_vec = {air_eos, davis_r_eos, jwl_eos, gruneisen_eos
// #ifdef SINGULARITY_TEST_SESAME
// #ifdef SINGULARITY_USE_SPINER_WITH_HDF5
//                                 ,
//                                 spiner_eos
// #endif
// #endif
    };

    const size_t n_eos = eos_vec.size();
    EOS *v_EOS = copy_eos_arr_to_device<decltype(eos_vec), EOS>(n_eos, eos_vec);

    WHEN("The energy is associated with a temperature and density above the reference") {

      constexpr Real T_lookup = 500;
      constexpr Real e_offset = 0.;    // No offset
      constexpr Real rho_factor = 1.2; // Slightly larger than reference density

      constexpr Real tol = 1.0e-14;

      THEN("P(rho, e) lookups should agree with P(rho, T) lookups when the energy is "
           "floored") {

        auto P_diffs = diff_pressures(n_eos, v_EOS, T_lookup, e_offset, rho_factor);

        for (size_t i = 0; i < n_eos; i++) {
          GetName evaluate_func{};
          eos_vec[i].EvaluateHost(evaluate_func);
          INFO("EOS " << evaluate_func.name << " index: " << i);
          CHECK(fabs(P_diffs[i]) < tol);
        }
      }
    }
    WHEN("The energy is below the zero Kelvin isotherm at the reference density") {

      constexpr Real T_lookup = 0;
      constexpr Real e_offset = 1.0e10; // Large offset
      constexpr Real rho_factor = 1.0;  // Exactly reference density

      constexpr Real tol = 1.0e-14;

      THEN("P(rho, e) lookups should agree with P(rho, T) lookups") {

        auto P_diffs = diff_pressures(n_eos, v_EOS, T_lookup, e_offset, rho_factor);

        for (size_t i = 0; i < n_eos; i++) {
          GetName evaluate_func{};
          eos_vec[i].EvaluateHost(evaluate_func);
          INFO("EOS " << evaluate_func.name << " index: " << i);
          CHECK(fabs(P_diffs[i]) < tol);
        }
      }
    }
  }
}
