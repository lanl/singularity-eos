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

#include <array>
#include <string>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <ports-of-call/portability.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <test/eos_unit_test_helpers.hpp>

using namespace singularity;

constexpr Real GPa = 1.0e10;
constexpr Real MJ_per_kg = 1.0e10;

#ifdef SINGULARITY_BUILD_CLOSURE
#ifdef SINGULARITY_TEST_SESAME
#ifdef SINGULARITY_USE_SPINER
#ifdef SPINER_USE_HDF

SCENARIO("Density-Temperature PTE Solver", "[PTE]") {

  GIVEN("Equations of state") {
    // Set up the three EOS
    constexpr int num_eos = 3;
    std::array<EOS, num_eos> eos_arr;
    // Reference state
    constexpr Real P0 = 1.0e6; // 1 bar
    constexpr Real T0 = 296;   // K
    // Ideal gas air
    constexpr Real rho0_air = 1e-03; // g/cc
    constexpr Real Gruneisen_air = 0.4;
    constexpr Real CV_air = P0 / rho0_air / (Gruneisen_air * T0);
    eos_arr[0] = IdealGas(Gruneisen_air, CV_air);
    // Spiner copper EOS
    constexpr int Cu_matid = 3337;
    const std::string eos_file = "../materials.sp5";
    eos_arr[1] = SpinerEOSDependsRhoT(eos_file, Cu_matid);
    // Davis Products EOS
    constexpr Real a = 0.798311;
    constexpr Real b = 0.58;
    constexpr Real k = 1.35;
    constexpr Real n = 2.66182;
    constexpr Real vc = 0.75419;              // cm^3 / g
    constexpr Real pc = 3.2 * GPa;            // microbar
    constexpr Real Cv = 0.001072 * MJ_per_kg; // erg / g / K
    constexpr Real Q = 4.115 * MJ_per_kg;
    eos_arr[2] = DavisProducts(a, b, k, n, vc, pc, Cv).Modify<ShiftedEOS>(-Q);

    // Move EOS array from host to device
    for (auto eos : eos_arr) {
      eos.GetOnDevice();
    }
    constexpr size_t EOS_bytes = num_eos * sizeof(EOS);
    EOS *v_EOS = (EOS *)PORTABLE_MALLOC(EOS_bytes);
    portableCopyToDevice(v_EOS, eos_arr.data(), EOS_bytes);

    GIVEN("A difficult state for the density-temperature PTE solver") {
      // Define the state
      constexpr Real spvol_bulk = 6.256037280402093e-01;
      constexpr Real sie_bulk = -1.441692060406520e+10;
      const std::array<const Real, num_eos> mass_frac = {
          1.10382442033331e-10, 0.124935312146569, 0.875064687743048};

      // Initial guess: all materials at the cell density
      const std::array<const Real, num_eos> vol_frac = mass_frac;

      // Calculate material densities (and corresponding energies) and the total
      // volume fraction
      std::array<Real, num_eos> densities;
      std::array<Real, num_eos> sies;
      Real vfrac_sum = 0;
      for (auto i = 0; i < num_eos; ++i) {
        densities[i] = mass_frac[i] / spvol_bulk / vol_frac[i];
        sies[i] = sie_bulk * mass_frac[i];
        vfrac_sum += vol_frac[i];
      }

      // Initialize pressure and temperature arrays to zero
      std::array<Real, num_eos> temperatures = {0., 0., 0.};
      std::array<Real, num_eos> pressures = {0., 0., 0.};

      // Copy values to device (when available)
      constexpr size_t bytes = num_eos * sizeof(Real);
      Real *v_densities = (Real *)PORTABLE_MALLOC(bytes);
      portableCopyToDevice(v_densities, densities.data(), bytes);
      Real *v_vol_frac = (Real *)PORTABLE_MALLOC(bytes);
      portableCopyToDevice(v_vol_frac, vol_frac.data(), bytes);
      Real *v_sies = (Real *)PORTABLE_MALLOC(bytes);
      portableCopyToDevice(v_sies, sies.data(), bytes);
      Real *v_temperatures = (Real *)PORTABLE_MALLOC(bytes);
      portableCopyToDevice(v_temperatures, temperatures.data(), bytes);
      Real *v_pressures = (Real *)PORTABLE_MALLOC(bytes);
      portableCopyToDevice(v_pressures, pressures.data(), bytes);

      THEN("The PTE solver should converge") {
        // Allocate scratch space for the PTE solver
        int pte_solver_scratch_size = PTESolverRhoTRequiredScratch(num_eos);
        double *scratch = (double *)PORTABLE_MALLOC(pte_solver_scratch_size);
        double **lambdas = (double **)PORTABLE_MALLOC(MAX_NUM_LAMBDAS * num_eos);
        PTESolverRhoT<decltype(v_EOS), Real *, Real **> method(
            num_eos, v_EOS, vfrac_sum, sie_bulk, v_densities, v_vol_frac, v_sies,
            v_temperatures, v_pressures, lambdas, scratch);

        // Solve the PTE system and ensure it converged
        bool pte_converged = PTESolver(method);
        CHECK(pte_converged);

        // Free temp memory
        PORTABLE_FREE(lambdas);
        PORTABLE_FREE(scratch);
      }
    }
    // Call Finalize on each EOS
    for (auto eos : eos_arr) {
      eos.Finalize();
    }
    // Free EOS memory
    PORTABLE_FREE(v_EOS);
  }
}
#endif // SPINER_USE_HDF
#endif // SINGULARITY_USE_SPINER
#endif // SINGULARITY_TEST_SESAME
#endif // SINGULARITY_BUILD_CLOSURE
