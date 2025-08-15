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

#ifdef SINGULARITY_BUILD_CLOSURE

#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/indexable_types.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <test/eos_unit_test_helpers.hpp>

#ifdef SINGULARITY_TEST_SESAME
#ifdef SINGULARITY_USE_SPINER_WITH_HDF5

constexpr Real GPa = 1.0e10;
constexpr Real MJ_per_kg = 1.0e10;

using singularity::DavisProducts;
using singularity::DavisReactants;
using singularity::IdealGas;
using singularity::MAX_NUM_LAMBDAS;
using singularity::MixParams;
using singularity::PTESolverPT;
using singularity::PTESolverPTRequiredScratch;
using singularity::PTESolverRhoT;
using singularity::PTESolverRhoTRequiredScratch;
using singularity::ShiftedEOS;
using singularity::SpinerEOSDependsRhoT;
using singularity::mix_impl::CacheAccessor;
using singularity::robust::ratio;

using EOS = singularity::Variant<IdealGas, ShiftedEOS<DavisProducts>, DavisProducts,
                                 DavisReactants, SpinerEOSDependsRhoT>;

template <int NEOS, typename... Ts>
struct LambdaIndexer {
 public:
  using Lambda_t = singularity::IndexerUtils::VariadicPointerIndexer<Ts...>;
  LambdaIndexer(Real *data) : data_(data) {}
  PORTABLE_FORCEINLINE_FUNCTION
  auto operator[](const std::size_t m) { return Lambda_t(data_ + Lambda_t::size() * m); }
  static inline constexpr std::size_t size() { return NEOS * Lambda_t::size(); }

 private:
  Real *data_;
};
// Specialized to the EOSs we actually want
template <int NEOS>
using MyLambdaIndexer = LambdaIndexer<NEOS, singularity::IndexableTypes::LogDensity,
                                      singularity::IndexableTypes::LogTemperature,
                                      singularity::IndexableTypes::RootStatus,
                                      singularity::IndexableTypes::TableStatus>;

// TODO(JMM): Clean this up to account for these two differetn APIs
template <template <typename... Types> class PTESolver_t, typename Scratch_t,
          typename ArrT>
bool run_PTE_from_state(const int num_pte, EOS *v_EOS, const Real spvol_bulk,
                        const Real sie_bulk, const Scratch_t &&RequiredScratch,
                        ArrT mass_frac, Real &u_bulk_out) {
  // Calculate material densities (and corresponding energies) and the total
  // volume fraction
  std::vector<Real> vol_frac(num_pte);
  std::vector<Real> densities(num_pte);
  std::vector<Real> sies(num_pte);
  std::vector<Real> temperatures(num_pte);
  std::vector<Real> pressures(num_pte);
  Real vfrac_sum = 0;
  for (auto i = 0; i < num_pte; ++i) {
    // Initial guess: all materials at the cell density
    vol_frac[i] = mass_frac[i];
    densities[i] = mass_frac[i] / spvol_bulk / vol_frac[i];
    sies[i] = sie_bulk * mass_frac[i];
    vfrac_sum += vol_frac[i];
    // Initialize pressure and temperature arrays to zero
    temperatures[i] = 0;
    pressures[i] = 0;
  }
  for (auto i = 0; i < num_pte; ++i) {
    vol_frac[i] /= vfrac_sum;
    densities[i] *= vfrac_sum;
  }
  vfrac_sum = 1.0;

  // Copy values to device (when available)
  const size_t bytes = num_pte * sizeof(Real);
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
  const int pte_solver_scratch_size = RequiredScratch(num_pte);
  const size_t scratch_bytes = pte_solver_scratch_size * sizeof(Real);
  Real *scratch = (double *)PORTABLE_MALLOC(scratch_bytes);

  // Allocate lambdas for all EOS and use an accessor to index into it
  const size_t lambda_bytes = num_pte * MAX_NUM_LAMBDAS * sizeof(Real);
  Real *lambda_memory = (Real *)PORTABLE_MALLOC(lambda_bytes);

  // Solve the PTE system on device using a one-teration portableFor
  bool pte_converged;
  constexpr size_t bool_bytes = 1 * sizeof(bool);
  bool *pte_converged_d = (bool *)PORTABLE_MALLOC(bool_bytes);
  constexpr size_t real_bytes = 1 * sizeof(Real);
  Real *pte_u_out_d = (Real *)PORTABLE_MALLOC(real_bytes);
  portableFor(
      "Device execution of PTE Test", 0, 1, PORTABLE_LAMBDA(int i) {
        CacheAccessor lambdas(lambda_memory);
        PTESolver_t<decltype(v_EOS), Real *, decltype(lambdas)> method(
            num_pte, v_EOS, vfrac_sum, sie_bulk, v_densities, v_vol_frac, v_sies,
            v_temperatures, v_pressures, lambdas, scratch);
        auto status = PTESolver(method);
        pte_converged_d[0] = status.converged;
        pte_u_out_d[0] = v_densities[0] * v_vol_frac[0] * v_sies[0];
        for (int m = 1; m < num_pte; ++m) {
          pte_u_out_d[0] += v_densities[m] * v_vol_frac[m] * v_sies[m];
        }
      });
  portableCopyToHost(&pte_converged, pte_converged_d, bool_bytes);
  portableCopyToHost(&u_bulk_out, pte_u_out_d, real_bytes);

  // Free temp memory
  PORTABLE_FREE(lambda_memory); // Free entire lambda array
  PORTABLE_FREE(scratch);

  // Free PTE values
  PORTABLE_FREE(v_densities);
  PORTABLE_FREE(v_vol_frac);
  PORTABLE_FREE(v_sies);
  PORTABLE_FREE(v_temperatures);
  PORTABLE_FREE(v_pressures);

  return pte_converged;
}

// TODO(JMM): Clean this up to account for these two differetn APIs
inline bool RunPTE2Mat(EOS eos1, EOS eos2, Real rhobar1, Real rhobar2, Real alpha_guess1,
                       Real alpha_guess2, Real sietot, Real Tguess, Real alpha1_true,
                       Real alpha2_true, Real Ttrue, Real Ptrue) {
  constexpr int NT = 1;
  constexpr int NEOS = 2;

  EOS *eos = (EOS *)PORTABLE_MALLOC(NEOS * sizeof(EOS));
  Real *rho = (Real *)PORTABLE_MALLOC(NEOS * sizeof(Real *));
  Real *alpha = (Real *)PORTABLE_MALLOC(NEOS * sizeof(Real *));
  Real *sie = (Real *)PORTABLE_MALLOC(NEOS * sizeof(Real *));
  Real *temp = (Real *)PORTABLE_MALLOC(NEOS * sizeof(Real *));
  Real *press = (Real *)PORTABLE_MALLOC(NEOS * sizeof(Real *));
  Real *plambda = (Real *)PORTABLE_MALLOC(MyLambdaIndexer<NEOS>::size() * sizeof(Real));

  // PTE solvers require internal scratch space. However, the solver
  // doesn't manage memory. We must provide it ourselves.
  const std::size_t pte_scratch_size = singularity::PTESolverRhoTRequiredScratch(NEOS);
  Real *pscratch = (Real *)PORTABLE_MALLOC(pte_scratch_size * sizeof(Real));

  // The pte_params object contains a number of settings you can
  // modify for the PTE solver, such as tolerances.
  singularity::MixParams pte_params;
  pte_params.pte_rel_tolerance_p = 1e-12;
  pte_params.pte_rel_tolerance_e = 1e-12;
  pte_params.pte_abs_tolerance_p = 0;

  int nsuccess = 0;
  portableReduce(
      "Run a difficult 2 material state", 0, NT,
      PORTABLE_LAMBDA(const int i, int &ns) {
        Real vfrac_sum = alpha_guess1 + alpha_guess2;
        alpha[0] = alpha_guess1 / vfrac_sum;
        alpha[1] = alpha_guess2 / vfrac_sum;
        rho[0] = vfrac_sum * rhobar1 / alpha_guess1;
        rho[1] = vfrac_sum * rhobar2 / alpha_guess2;

        eos[0] = eos1;
        eos[1] = eos2;
        MyLambdaIndexer<NEOS> lambda(plambda);

        singularity::PTESolverRhoT<EOS *, Real *, MyLambdaIndexer<NEOS>> method(
            NEOS, eos, 1.0, sietot, rho, alpha, sie, temp, press, lambda, pscratch,
            Tguess, pte_params);
        // Run the solver
        auto status = singularity::PTESolver(method);

        bool success = true;
        if (!status.converged) {
          printf("Solver did not converge!\n");
          success = false;
        } else {
          if (!(isClose(alpha[0], alpha1_true) && isClose(alpha[1], alpha2_true))) {
            printf("Volume fractions dop not match! [%.14e %.14e], [%.14e %.14e]\n",
                   alpha[0], alpha[1], alpha1_true, alpha2_true);
            success = false;
          }
          if (!isClose(temp[0], Ttrue)) {
            printf("Temperatures do not match! %.14e %.14e\n", temp[0], Ttrue);
            success = false;
          }
          if (!isClose(press[0], Ptrue)) {
            printf("Pressures do not match! %.14e %.14e\n", press[0], Ptrue);
            success = false;
          }
        }
        ns += success;
      },
      nsuccess);

  PORTABLE_FREE(eos);
  PORTABLE_FREE(rho);
  PORTABLE_FREE(alpha);
  PORTABLE_FREE(sie);
  PORTABLE_FREE(temp);
  PORTABLE_FREE(press);
  PORTABLE_FREE(plambda);

  return (nsuccess == NT);
}

SCENARIO("Density- and Pressure-Temperature PTE Solvers", "[PTE]") {

  GIVEN("Four equations of state") {
    // Reference state
    constexpr Real P0 = 1.0e6; // 1 bar
    constexpr Real T0 = 296;   // K
    // Ideal gas air
    constexpr Real rho0_air = 1e-03; // g/cc
    constexpr Real Gruneisen_air = 0.4;
    constexpr Real Cv_air = P0 / rho0_air / (Gruneisen_air * T0);
    EOS air_eos = IdealGas(Gruneisen_air, Cv_air);
    // Spiner copper EOS
    constexpr int Cu_matid = 3337;
    const std::string eos_file = "../materials.sp5";
    EOS copper_eos = SpinerEOSDependsRhoT(eos_file, Cu_matid);
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
    EOS davis_r_eos =
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
    EOS davis_p_eos = DavisProducts(a, b, k, n, vc, pc, Cv).Modify<ShiftedEOS>(-Q);

    const MixParams params;
    GIVEN("A difficult three-material state") {
      // Define the state
      constexpr Real spvol_bulk = 6.256037280402093e-01;
      constexpr Real sie_bulk = -1.441692060406520e+10;
      constexpr int num_pte = 3;
      const std::array<const Real, num_pte> mass_frac = {
          1.10382442033331e-10, 0.124935312146569, 0.875064687743048};
      std::array<EOS, num_pte> eos_arr = {air_eos.GetOnDevice(), copper_eos.GetOnDevice(),
                                          davis_p_eos.GetOnDevice()};

      THEN("The PTE Rho T solver should converge") {
        EOS *v_EOS = copy_eos_arr_to_device<decltype(eos_arr), EOS>(num_pte, eos_arr);
        Real u_bulk_out = std::numeric_limits<Real>::max();
        const bool pte_converged = run_PTE_from_state<PTESolverRhoT>(
            num_pte, v_EOS, spvol_bulk, sie_bulk, PTESolverRhoTRequiredScratch, mass_frac,
            u_bulk_out);
        CHECK(pte_converged);
        AND_THEN("The solution satisfies the bulk internal energy constraint") {
          const Real u_bulk = ratio(sie_bulk, spvol_bulk);
          const Real u_scale = std::abs(u_bulk);
          const Real u_bulk_scale = ratio(u_bulk, u_scale);
          const Real residual = std::abs(u_bulk_scale - ratio(u_bulk_out, u_scale));
          CHECK(residual < params.pte_rel_tolerance_e);
        }
        // Free EOS copies on device
        PORTABLE_FREE(v_EOS);
      }
      THEN("The PTE P T solver should converge") {
        EOS *v_EOS = copy_eos_arr_to_device<decltype(eos_arr), EOS>(num_pte, eos_arr);
        Real u_bulk_out = std::numeric_limits<Real>::max();
        const bool pte_converged = run_PTE_from_state<PTESolverPT>(
            num_pte, v_EOS, spvol_bulk, sie_bulk, PTESolverPTRequiredScratch, mass_frac,
            u_bulk_out);
        CHECK(pte_converged);
        AND_THEN("The solution satisfies the bulk internal energy constraint") {
          const Real u_bulk = ratio(sie_bulk, spvol_bulk);
          const Real u_scale = std::abs(u_bulk);
          const Real u_bulk_scale = ratio(u_bulk, u_scale);
          const Real residual = std::abs(u_bulk_scale - ratio(u_bulk_out, u_scale));
          CHECK(residual < params.pte_rel_tolerance_e);
        }
        // Free EOS copies on device
        PORTABLE_FREE(v_EOS);
      }
      finalize_eos_arr(eos_arr);
    }
    GIVEN("A difficult two-material state") {
      // TODO: make this test converge
      // Define the state
      constexpr Real spvol_bulk = 4.010467628234189e-01;
      constexpr Real sie_bulk = 3.290180957185173e+07;
      constexpr int num_pte = 2;
      const std::array<const Real, num_pte> mass_frac = {0.000312273191678158,
                                                         0.999687726808322};
      std::array<EOS, num_pte> eos_arr = {air_eos.GetOnDevice(),
                                          copper_eos.GetOnDevice()};
      // TODO(JMM): This test does not converge. See Issue 390. Possibly due to Spiner?
      THEN("The PTE solver should converge") {
        EOS *v_EOS = copy_eos_arr_to_device<decltype(eos_arr), EOS>(num_pte, eos_arr);
        Real u_bulk_out = std::numeric_limits<Real>::max();
        const bool pte_converged = run_PTE_from_state<PTESolverRhoT>(
            num_pte, v_EOS, spvol_bulk, sie_bulk, PTESolverRhoTRequiredScratch, mass_frac,
            u_bulk_out);
        CHECK(pte_converged);
        const MixParams params;
        AND_THEN("The solution satisfies the bulk internal energy constraint") {
          // NOTE(@pdmullen): The following fails prior to PR401
          const Real u_bulk = ratio(sie_bulk, spvol_bulk);
          const Real u_scale = std::abs(u_bulk);
          const Real u_bulk_scale = ratio(u_bulk, u_scale);
          const Real residual = std::abs(u_bulk_scale - ratio(u_bulk_out, u_scale));
          CHECK(residual < params.pte_rel_tolerance_e);
        }
        // Free EOS copies on device
        PORTABLE_FREE(v_EOS);
      }
      // Deallocates device memory for each EOS (if applicable)
      finalize_eos_arr(eos_arr);
    }
    // Clean up original EOS objects before they go out of scope
    air_eos.Finalize();     // irrelevant because no data allocated
    copper_eos.Finalize();  // actually does something
    davis_r_eos.Finalize(); // irrelevant because no data allocated
    davis_p_eos.Finalize(); // irrelevant because no data allocated
  }

  GIVEN("A difficult two material mixed cell") {
    const std::string eos_file = "../materials.sp5";

    constexpr int Al_matid = 3720;
    EOS al_eos_h = SpinerEOSDependsRhoT(eos_file, Al_matid);
    EOS al_eos = al_eos_h.GetOnDevice();

    constexpr int foam_matid = 7592;
    EOS foam_eos_h = SpinerEOSDependsRhoT(eos_file, foam_matid);
    EOS foam_eos = foam_eos_h.GetOnDevice();

    constexpr Real rhobar1 = 4.84896778419434e-04;
    constexpr Real rhobar2 = 3.85334470246594e-02;
    constexpr Real alpha_guess1 = 1.87984334110313e-03;
    constexpr Real alpha_guess2 = 9.98120156658897e-01;
    constexpr Real sietot = 9.38106182865728e+11;
    constexpr Real Tguess = 1.49041098541734e+03;

    WHEN("We call PTE") {
      constexpr Real alpha1_true = 2.65277486969419e-03;
      constexpr Real alpha2_true = 9.97347225130306e-01;
      constexpr Real Ttrue = 2.10866067749579e+04;
      constexpr Real Ptrue = 1.44504093939007e+10;

      bool success =
          RunPTE2Mat(al_eos, foam_eos, rhobar1, rhobar2, alpha_guess1, alpha_guess2,
                     sietot, Tguess, alpha1_true, alpha2_true, Ttrue, Ptrue);

      THEN("The solver converges") { REQUIRE(success); }
    }

    al_eos_h.Finalize();
    al_eos.Finalize();
    foam_eos_h.Finalize();
    foam_eos.Finalize();
  }

  GIVEN("A cell containing a difficult mixture with a tabulated and ideal gas EOS") {
    const std::string eos_file = "../materials.sp5";
    constexpr Real gm1 = 0.666666666666667;
    constexpr Real Cv = 3.1e7;
    EOS He_eos_h = IdealGas(gm1, Cv);
    EOS He_eos = He_eos_h.GetOnDevice();

    constexpr int foam_matid = 7592;
    EOS foam_eos_h = SpinerEOSDependsRhoT(eos_file, foam_matid);
    EOS foam_eos = foam_eos_h.GetOnDevice();

    WHEN("We request density at minimum pressure") {
      Real rhopmin_100 = foam_eos_h.RhoPmin(100);
      THEN("It's zero, as appropriate for a Maxwell constructed EOS") {
        REQUIRE(rhopmin_100 == 0);
      }
    }

    constexpr Real rhobar1 = 1.59761356859602e-04;
    constexpr Real rhobar2 = 7.81928505957464e-09;
    constexpr Real alpha_guess1 = 9.99999776549350e-01;
    constexpr Real alpha_guess2 = 2.23450650410351e-07;
    constexpr Real sietot = 9.37635921815154e+09;
    constexpr Real Tguess = 293;

    WHEN("We call PTE") {
      constexpr Real alpha1_true = 9.99999989159325e-01;
      constexpr Real alpha2_true = 1.08406754944454e-08;
      constexpr Real Ttrue = 3.02477729062119e+02;
      constexpr Real Ptrue = 9.98701227398616e+05;

      bool success =
          RunPTE2Mat(He_eos, foam_eos, rhobar1, rhobar2, alpha_guess1, alpha_guess2,
                     sietot, Tguess, alpha1_true, alpha2_true, Ttrue, Ptrue);

      THEN("The solver converges") { REQUIRE(success); }

      AND_WHEN("We call PTE but with a wildly incorrect temperature guess, deep in the "
               "Maxwell constructed region") {
        constexpr Real alpha1_true = 9.99999989159325e-01;
        constexpr Real alpha2_true = 1.08406754944454e-08;
        constexpr Real Ttrue = 3.02477729062119e+02;
        constexpr Real Ptrue = 9.98701227398616e+05;

        bool success =
            RunPTE2Mat(He_eos, foam_eos, rhobar1, rhobar2, alpha_guess1, alpha_guess2,
                       sietot, 5, alpha1_true, alpha2_true, Ttrue, Ptrue);

        THEN("The solver converges") { REQUIRE(success); }
      }
    }

    He_eos_h.Finalize();
    He_eos.Finalize();
    foam_eos_h.Finalize();
    foam_eos.Finalize();
  }
}
#endif // SINGULARITY_USE_SPINER_WITH_HDF5
#endif // SINGULARITY_TEST_SESAME
#endif // SINGULARITY_BUILD_CLOSURE
