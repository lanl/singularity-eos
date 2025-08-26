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

// This example demonstrates how to call the PTE solver for two
// materials. The core idea is that you have two materials who's
// volume averaged densities are
//
// rhobar1 = alpha1 rho1
// rhobar2 = alpha2 rho2
//
// where alpha1 and alpha2 are volume fractions such that
//
// alpha1 + alpha2 = 1
//
// and rho1 and rho2 are the microphysical densities.
//
// We seek to find the alpha1, alpha2 and temperature T such that the
// pressures P1 and P2 (of each materal 1 and 2) are equal:
//
// P1(rho1, T) = P2(rho2, T)
//
// and that the total volume-averaged energy density equals some
// predefined total:
//
// utot = alpha1 rho1 e(rho1, T) + alpha2 rho2 e(rho2, T)
//
// where here e1 and e2 are the specific internal energies opf each
// material. We also define
//
// rhobartot = rhobar1 + rhobar2
//
// and
//
// etot = utot / rhotot
//
// This problem is claled pressure-temperature equilibrium. And it
// generalizes to an arbitrary number of materials. Here we show how
// to do it with just two.

// C headers
#include <cmath>
#include <cstdio>

// C++ headers
#include <chrono>
#include <iostream>
#include <memory>
#include <string>

// This library contains portable utilities
#include <ports-of-call/portability.hpp>

// This contains logic for indexers
#include <singularity-eos/base/indexable_types.hpp>
// This contains useful tools for preventing things like divide by zero
#include <singularity-eos/base/robust_utils.hpp>
// The PTE closures
#include <singularity-eos/closure/mixed_cell_models.hpp>
// Needed to import the eos models
#include <singularity-eos/eos/eos.hpp>

using duration = std::chrono::microseconds;
inline auto now() { return std::chrono::high_resolution_clock::now(); }

// Number of equations of state to use. For this example, it's 2.
constexpr std::size_t NEOS = 2;
// Number of PTE solves to do
constexpr std::size_t NTRIAL = 1;

// Set the EOSs you want to use here.
using EOS =
    singularity::Variant<singularity::SpinerEOSDependsRhoT, singularity::IdealGas>;

// A LAMBDA contains additional arguments to an EOS, such as cached
// variables.  It accepts TYPES rather than integers, kind of like a
// map, as described in `IndexableTypes`. Here we use it and create a
// pair of these.
template <typename... Ts>
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
using MyLambdaIndexer = LambdaIndexer<
    singularity::IndexableTypes::LogDensity, singularity::IndexableTypes::LogTemperature,
    singularity::IndexableTypes::RootStatus, singularity::IndexableTypes::TableStatus>;

int main(int argc, char *argv[]) {
  if (argc < 8) {
    std::cerr << "Usage: " << argv[0]
              << " rhobar1 rhobar2 vfrac_guess 1 vfrac_guess_2 sietot Tguess eosargs..."
              << std::endl;
    std::exit(1);
  }

  // This is needed for Kokkos
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
#endif
  {
    const Real rhobar1 = atof(argv[1]);
    const Real rhobar2 = atof(argv[2]);
    const Real alpha_guess1 = atof(argv[3]);
    const Real alpha_guess2 = atof(argv[4]);
    const Real sietot = atof(argv[5]);
    const Real Tguess = atof(argv[6]);

    const Real rhotot = rhobar1 + rhobar2;
    const Real utot = sietot * rhotot;

    // modify for your needs.
    const std::string materialsfile = argv[7];
    const Real gm1 = atof(argv[8]);
    const Real Cv = atof(argv[9]);
    const int matid = atoi(argv[10]);

    // Now let's load up the EOS's.
    EOS eos1 = singularity::IdealGas(gm1, Cv);
    EOS eos2 = singularity::SpinerEOSDependsRhoT(materialsfile, matid);

    // and move them to device
    eos1 = eos1.GetOnDevice();
    eos2 = eos2.GetOnDevice();

    // We'll also need a device pointer for our EOS objects as the PTE
    // solver expects an array-like data structure. We're only going
    // to modify these arrays on device, however.
    EOS *eos = (EOS *)PORTABLE_MALLOC(NEOS * sizeof(EOS));

    // The volume fractions and microphysical material densities will
    // be modified "in place" on device. So lets create some data
    // structures for that.
    // PORTABLE_MALLOC is just like Malloc but for GPUs
    Real *rho = (Real *)PORTABLE_MALLOC(NEOS * sizeof(Real *));
    Real *alpha = (Real *)PORTABLE_MALLOC(NEOS * sizeof(Real *));
    // we need sie, temperature, pressure arrays
    Real *sie = (Real *)PORTABLE_MALLOC(NEOS * sizeof(Real *));
    Real *temp = (Real *)PORTABLE_MALLOC(NEOS * sizeof(Real *));
    Real *press = (Real *)PORTABLE_MALLOC(NEOS * sizeof(Real *));

    // The memory for the lambda array
    Real *plambda = (Real *)PORTABLE_MALLOC(MyLambdaIndexer::size() * sizeof(Real));

    // PTE solvers require internal scratch space. However, the solver
    // doesn't manage memory. We must provide it ourselves.
    const std::size_t pte_scratch_size = singularity::PTESolverRhoTRequiredScratch(NEOS);
    Real *pscratch = (Real *)PORTABLE_MALLOC(pte_scratch_size * sizeof(Real));

    // The pte_params object contains a number of settings you can
    // modify for the PTE solver, such as tolerances.
    singularity::MixParams pte_params;
    pte_params.pte_rel_tolerance_p = 1e-10;
    pte_params.pte_rel_tolerance_e = 1e-10;
    pte_params.pte_rel_tolerance_t = 1e-10;
    pte_params.pte_rel_tolerance_p_sufficient = 1e-10;
    pte_params.pte_rel_tolerance_e_sufficient = 1e-10;
    pte_params.pte_rel_tolerance_t_sufficient = 1e-10;
    pte_params.pte_abs_tolerance_p = 0;
    pte_params.pte_abs_tolerance_t = 0;
    pte_params.pte_abs_tolerance_p_sufficient = 0;
    pte_params.pte_abs_tolerance_t_sufficient = 0;
    pte_params.pte_reduced_system_exclude_idx = -1;

    // Now we call the device kernel
    auto start = now();
    portableFor(
        "Run a PTE solve", 0, NTRIAL, PORTABLE_LAMBDA(const int i) {
          // Prepare microphysical densities and initial guesses.  We
          // want to normalize initial guess volume fractions so they
          // sum to 1.
          Real vfrac_sum = alpha_guess1 + alpha_guess2;
          alpha[0] = alpha_guess1 / vfrac_sum;
          alpha[1] = alpha_guess2 / vfrac_sum;
          rho[0] = vfrac_sum * rhobar1 / alpha_guess1;
          rho[1] = vfrac_sum * rhobar2 / alpha_guess2;

          // Asign the EOS objects to the array
          eos[0] = eos1;
          eos[1] = eos2;

          // Repeare lambda indexer
          MyLambdaIndexer lambda(plambda);

          // Create PTE solver object
          // There are many different solvers available. Here we use the RhoT solver.
          singularity::PTESolverRhoT<EOS *, Real *, MyLambdaIndexer> method(
              NEOS, eos,
              1.0, // the 1 here is because we normalized volume fractions to sum to 1
              sietot, rho,
              alpha, // alpha is just an initial guess, but rho * alpha MUST equal rhobar
              sie, temp, press, // The PTE solver sets these. They're outputs.
              lambda, pscratch, // scratch and lambda
              Tguess,           // Initial guess for temperature. Also used
                                // as a normalization factor internal to
                                // the solver.
              pte_params);
          // Run the solver
          auto status = singularity::PTESolver(method);
          printf("Solver succeeded:      %d\n"
                 "Slow descent detected: %d\n"
                 "Num iterations:        %ld\n"
                 "Volume fractions:      %.14e %.14e\n"
                 "Temperature:           %.14e\n"
                 "Pressure:              %.14e\n",
                 status.converged, status.slow_convergence_detected, status.max_niter,
                 alpha[0], alpha[1], temp[0], press[0]);
        });
    auto stop = now();
    printf("Time in us = %ld\n",
           std::chrono::duration_cast<duration>(stop - start).count());

    eos1.Finalize();
    eos2.Finalize();
    PORTABLE_FREE(eos);
    PORTABLE_FREE(rho);
    PORTABLE_FREE(alpha);
    PORTABLE_FREE(sie);
    PORTABLE_FREE(temp);
    PORTABLE_FREE(press);
    PORTABLE_FREE(plambda);
    PORTABLE_FREE(pscratch);
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif
  return 0;
}
