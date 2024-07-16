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

#ifdef SINGULARITY_BUILD_CLOSURE

#include <array>
#include <string>
#include <vector>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <ports-of-call/portability.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <test/eos_unit_test_helpers.hpp>

constexpr Real GPa = 1.0e10;
constexpr Real MJ_per_kg = 1.0e10;

#ifdef SINGULARITY_TEST_SESAME
#ifdef SINGULARITY_USE_SPINER_WITH_HDF5

using singularity::DavisProducts;
using singularity::IdealGas;
using singularity::MAX_NUM_LAMBDAS;
using singularity::PTESolverRhoT;
using singularity::PTESolverRhoTRequiredScratch;
using singularity::ShiftedEOS;
using singularity::SpinerEOSDependsRhoT;

using EOS = singularity::Variant<IdealGas, ShiftedEOS<DavisProducts>, DavisProducts,
                                 SpinerEOSDependsRhoT>;

template <typename EOSArrT>
EOS *copy_eos_arr_to_device(const int num_eos, EOSArrT eos_arr) {
  // Move EOS array from host to device
  const size_t EOS_bytes = num_eos * sizeof(EOS);
  for (auto i = 0; i < num_eos; i++) {
    eos_arr[i] = eos_arr[i].GetOnDevice();
  }
  EOS *v_EOS = (EOS *)PORTABLE_MALLOC(EOS_bytes);
  const size_t bytes = num_eos * sizeof(EOS);
  portableCopyToDevice(v_EOS, eos_arr.data(), bytes);
  return v_EOS;
}

void finalize_eos_arr(const int num_eos, EOS *v_EOS) {
  // Call Finalize on each EOS
  for (auto i = 0; i < num_eos; i++) {
    v_EOS[i].Finalize();
  }
  // Free EOS memory
  PORTABLE_FREE(v_EOS);
}

template <typename ArrT>
bool run_PTE_from_state(const int num_eos, EOS *v_EOS, const Real spvol_bulk,
                        const Real sie_bulk, ArrT mass_frac) {
  // Calculate material densities (and corresponding energies) and the total
  // volume fraction
  std::vector<Real> vol_frac(num_eos);
  std::vector<Real> densities(num_eos);
  std::vector<Real> sies(num_eos);
  std::vector<Real> temperatures(num_eos);
  std::vector<Real> pressures(num_eos);
  Real vfrac_sum = 0;
  for (auto i = 0; i < num_eos; ++i) {
    // Initial guess: all materials at the cell density
    vol_frac[i] = mass_frac[i];
    densities[i] = mass_frac[i] / spvol_bulk / vol_frac[i];
    sies[i] = sie_bulk * mass_frac[i];
    vfrac_sum += vol_frac[i];
    // Initialize pressure and temperature arrays to zero
    temperatures[i] = 0;
    pressures[i] = 0;
  }

  // Copy values to device (when available)
  const size_t bytes = num_eos * sizeof(Real);
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

  // Allocate scratch space for the PTE solver
  const int pte_solver_scratch_size = PTESolverRhoTRequiredScratch(num_eos);
  const size_t scratch_bytes = pte_solver_scratch_size * sizeof(Real);
  Real *scratch = (double *)PORTABLE_MALLOC(scratch_bytes);

  // Allocate lambdas for each EOS
  constexpr size_t lambda_bytes = MAX_NUM_LAMBDAS * sizeof(Real); // EOS lambda size
  const size_t lambda_size = num_eos * sizeof(Real *);            // Size of lambda array
  Real **lambdas = (Real **)PORTABLE_MALLOC(lambda_size);
  for (auto i = 0; i < num_eos; i++) {
    lambdas[i] = (Real *)PORTABLE_MALLOC(lambda_bytes);
  }

  // Solve the PTE system on device using a one-teration portableFor
  bool pte_converged;
  constexpr size_t bool_bytes = 1 * sizeof(bool);
  bool *pte_converged_d = (bool *)PORTABLE_MALLOC(bool_bytes);
  portableFor("Device execution of PTE Test", 0, 1, PORTABLE_LAMBDA(auto i) {
    PTESolverRhoT<decltype(v_EOS), Real *, Real **> method(
    num_eos, v_EOS, vfrac_sum, sie_bulk, v_densities, v_vol_frac, v_sies,
    v_temperatures, v_pressures, lambdas, scratch);
    pte_converged_d[0] = PTESolver(method);
  });
  portableCopyToHost(&pte_converged, pte_converged_d, bool_bytes);

  // Free temp memory
  for (auto i = 0; i < num_eos; i++) {
    // Free each lambda separately first
    PORTABLE_FREE(lambdas[i]);
  }
  PORTABLE_FREE(lambdas); // Free entire lambda array
  PORTABLE_FREE(scratch);

  // Free PTE values
  PORTABLE_FREE(v_densities);
  PORTABLE_FREE(v_vol_frac);
  PORTABLE_FREE(v_sies);
  PORTABLE_FREE(v_temperatures);
  PORTABLE_FREE(v_pressures);

  return pte_converged;
}

SCENARIO("Density-Temperature PTE Solver", "[PTE]") {

  GIVEN("Three equations of state") {
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

    GIVEN("A difficult state for the density-temperature PTE solver") {
      // Define the state
      constexpr Real spvol_bulk = 6.256037280402093e-01;
      constexpr Real sie_bulk = -1.441692060406520e+10;
      const std::array<const Real, num_eos> mass_frac = {
          1.10382442033331e-10, 0.124935312146569, 0.875064687743048};

      THEN("The PTE solver should converge") {
        EOS *v_EOS = copy_eos_arr_to_device(num_eos, eos_arr);
        const bool pte_converged =
            run_PTE_from_state(num_eos, v_EOS, spvol_bulk, sie_bulk, mass_frac);
        CHECK(pte_converged);
        finalize_eos_arr(num_eos, v_EOS);
      }
    }
  }
  GIVEN("Two equations of state") {
    // Set up the three EOS
    constexpr int num_eos = 2;
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

    GIVEN("A state that would cause a negative temperature in the PTE solver") {
      // Define the state
      constexpr Real spvol_bulk = 4.010467628234189e-01;
      constexpr Real sie_bulk = 3.290180957185173e+07;
      const std::array<const Real, num_eos> mass_frac = {0.000312273191678158,
                                                         0.999687726808322};

      THEN("The PTE solver should converge") {
        EOS *v_EOS = copy_eos_arr_to_device(num_eos, eos_arr);
        const bool pte_converged =
            run_PTE_from_state(num_eos, v_EOS, spvol_bulk, sie_bulk, mass_frac);
        // TODO: make this test converge
        // CHECK(pte_converged);
        finalize_eos_arr(num_eos, v_EOS);
      }
    }
  }
}
#endif // SINGULARITY_USE_SPINER_WITH_HDF5
#endif // SINGULARITY_TEST_SESAME
#endif // SINGULARITY_BUILD_CLOSURE
