//------------------------------------------------------------------------------
// © 2021-2022. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_CLOSURE_MIXED_CELL_MODELS_
#define _SINGULARITY_EOS_CLOSURE_MIXED_CELL_MODELS_

#include <ports-of-call/portability.hpp>
#include <singularity-eos/eos/eos.hpp>

#include <cmath>

#ifdef SINGULARITY_USE_KOKKOSKERNELS
#include <KokkosBatched_ApplyQ_Decl.hpp>
#include <KokkosBatched_QR_Decl.hpp>
#include <KokkosBatched_Trsv_Decl.hpp>
#else
#include <Eigen/Dense>
#endif // SINGULARITY_USE_KOKKOSKERNELS

namespace singularity {

/*
  This is the PTE code based on FLAG's matrix form in their mixgas treatment
  I have modified it to include a line search which improves robustness and
  quality of non-converged solutions, at the likely cost of speed (due to
  additional EOS calls).
  It is templated on nmat and takes an array of Indexers to the EOSs,
  Volume (total volume of materials to be equilibrated), total SIE,
  and an array of masses of each component as inputs, and returns component
  volumes and energies (SIEs) as output.

  EOSIndexer must have an operator[](int) that returns an EOS. e.g., EOS*
  RealIndexer must have an operator[](int) that returns a Real. e.g., Real*
  ConstRealIndexer is as RealIndexer, but assumed const type.
  LambdaIndexer must have an operator[](int) that returns a Real*. e.g., Real**
*/
// Version templated on nmat
// niter version produces histogram
template <int nmat, typename EOSIndexer, typename RealIndexer, typename ConstRealIndexer,
          typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool pte_closure_flag_with_line(
    EOSIndexer &&eoss, const Real Volume, const Real TotalSIE,
    ConstRealIndexer &&ComponentMasses, RealIndexer &&ComponentVolumes,
    RealIndexer &&ComponentEnergies, LambdaIndexer &&lambdas, int &niter);
template <int nmat, typename EOSIndexer, typename RealIndexer, typename ConstRealIndexer,
          typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool
pte_closure_flag_with_line(EOSIndexer &&eoss, const Real Volume, const Real TotalSIE,
                           ConstRealIndexer &&ComponentMasses,
                           RealIndexer &&ComponentVolumes,
                           RealIndexer &&ComponentEnergies, LambdaIndexer &&lambdas) {
  int niter;
  return pte_closure_flag_with_line(eoss, Volume, TotalSIE, ComponentMasses,
                                    ComponentVolumes, ComponentEnergies, lambdas, niter);
}
// Version with nmat available at runtime
template <typename EOSIndexer, typename RealIndexer, typename ConstRealIndexer,
          typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool
pte_closure_flag(int nmat, EOSIndexer &&eoss, const Real Volume, const Real TotalSIE,
                 ConstRealIndexer &&ComponentMasses, RealIndexer &&ComponentVolumes,
                 RealIndexer &&ComponentEnergies, LambdaIndexer &&lambdas, int &niter);
template <typename EOSIndexer, typename RealIndexer, typename ConstRealIndexer,
          typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool
pte_closure_flag(int nmat, EOSIndexer &&eoss, const Real Volume, const Real TotalSIE,
                 ConstRealIndexer &&ComponentMasses, RealIndexer &&ComponentVolumes,
                 RealIndexer &&ComponentEnergies, LambdaIndexer &&lambdas) {
  int niter;
  return pte_closure_flag(nmat, eoss, Volume, TotalSIE, ComponentMasses, ComponentVolumes,
                          ComponentEnergies, lambdas, niter);
}
// Pointer-only version with offset for EOS indexing
PORTABLE_INLINE_FUNCTION bool
pte_closure_flag_offset(int nmat, EOS *eoss, const Real Volume, const Real TotalSIE,
                        int const *const ComponentMats, const Real *ComponentMasses,
                        Real *ComponentVolumes, Real *ComponentEnergies, Real **lambdas,
                        int &niter);
PORTABLE_INLINE_FUNCTION bool
pte_closure_flag_offset(int nmat, EOS *eoss, const Real Volume, const Real TotalSIE,
                        int const *const ComponentMats, const Real *ComponentMasses,
                        Real *ComponentVolumes, Real *ComponentEnergies, Real **lambdas) {
  int niter;
  return pte_closure_flag_offset(nmat, eoss, Volume, TotalSIE, ComponentMats,
                                 ComponentMasses, ComponentVolumes, ComponentEnergies,
                                 lambdas, niter);
}
/*
  This solver was developed by Josh Dolence, and is guaranteed to
  return a state where energies and volume fractions add up to the
  total energy and 1 respectively, even if the solve fails.

  EOSIndexer must have an operator[](int) that returns an EOS. e.g., EOS*
  RealIndexer must have an operator[](int) that returns a Real. e.g., Real*
  LambdaIndexer must have an operator[](int) that returns a Real*. e.g., Real**
*/
// Version templated on nmat.
template <int nmat, typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool pte_closure_josh(EOSIndexer &&eos, const Real vfrac_tot,
                                               const Real sie_tot, RealIndexer &&rho,
                                               RealIndexer &&vfrac, RealIndexer &&sie,
                                               RealIndexer &&temp, RealIndexer &&press,
                                               LambdaIndexer &&lambda, int &niter);
template <int nmat, typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool
pte_closure_josh(EOSIndexer &&eos, const Real vfrac_tot, const Real sie_tot,
                 RealIndexer &&rho, RealIndexer &&vfrac, RealIndexer &&sie,
                 RealIndexer &&temp, RealIndexer &&press, LambdaIndexer &&lambda) {
  int niter;
  return pte_closure_josh(eos, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press, lambda,
                          niter);
}
// Version with nmat available at runtime.
template <typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool
pte_closure_josh(int nmat, EOSIndexer &&eoss, const Real vfrac_tot, const Real sie_tot,
                 RealIndexer &&rho, RealIndexer &&vfrac, RealIndexer &&sie,
                 RealIndexer &&temp, RealIndexer &&press, LambdaIndexer &&lambdas,
                 int &niter);
template <typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool
pte_closure_josh(int nmat, EOSIndexer &&eoss, const Real vfrac_tot, const Real sie_tot,
                 RealIndexer &&rho, RealIndexer &&vfrac, RealIndexer &&sie,
                 RealIndexer &&temp, RealIndexer &&press, LambdaIndexer &&lambdas) {
  int niter;
  return pte_closure_josh(nmat, eoss, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press,
                          lambdas, niter);
}
// Pointer-only version with offset for EOS indexing
PORTABLE_INLINE_FUNCTION bool
pte_closure_josh_offset(int nmat, EOS *eoss, const Real vfrac_tot, const Real sie_tot,
                        int const *const Mats, Real *rho, Real *vfrac, Real *sie,
                        Real *temp, Real *press, Real **lambdas, int &niter);
PORTABLE_INLINE_FUNCTION bool
pte_closure_josh_offset(int nmat, EOS *eoss, const Real vfrac_tot, const Real sie_tot,
                        int const *const Mats, Real *rho, Real *vfrac, Real *sie,
                        Real *temp, Real *press, Real **lambdas) {
  int niter;
  return pte_closure_josh_offset(nmat, eoss, vfrac_tot, sie_tot, Mats, rho, vfrac, sie,
                                 temp, press, lambdas, niter);
}

/*
    This is a formulation of the PTE closure equations which tries to minimize
    differences between T and P/T.  Since T (should be) positive definite, this
    should be well-founded.  In mixtures of ideal gasses, this should be linear
    (so for hot plasmas, this should be linear-ish).
    Equations 1-nmat will be:
    T1-T0 = 0, T2-T1=0, [...], Tnmat-Tnmatm1=0, sum(Ei) = Etot
    and Equations nmat+1-2*nmat will be:
    P1/T1-P0/T0 = 0, P2/T2-P1/T1 = 0 [...], Sum(volumes)=Volume

    EOSIndexer must have an operator[](int) that returns an EOS. e.g., EOS*
    RealIndexer must have an operator[](int) that returns a Real. e.g., Real*
    ConstRealIndexer is as RealIndexer, but assumed const type.
    LambdaIndexer must have an operator[](int) that returns a Real*. e.g., Real**
*/
template <typename EOSIndexer, typename RealIndexer, typename ConstRealIndexer,
          typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool
pte_closure_ideal(int nmat, EOSIndexer &&eoss, const Real Volume, const Real TotalSIE,
                  ConstRealIndexer &&ComponentMasses, RealIndexer &&ComponentVolumes,
                  RealIndexer &&ComponentEnergies, LambdaIndexer &&lambdas, int &niter);
template <typename EOSIndexer, typename RealIndexer, typename ConstRealIndexer,
          typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool
pte_closure_ideal(int nmat, EOSIndexer &&eoss, const Real Volume, const Real TotalSIE,
                  ConstRealIndexer &&ComponentMasses, RealIndexer &&ComponentVolumes,
                  RealIndexer &&ComponentEnergies, LambdaIndexer &&lambdas) {
  int niter;
  return pte_closure_ideal(nmat, eoss, Volume, TotalSIE, ComponentMasses,
                           ComponentVolumes, ComponentEnergies, lambdas, niter);
}

// ======================================================================
// Implementation details below
// ======================================================================

namespace mix_params {
constexpr Real derivative_eps = 3.0e-6;
constexpr Real pte_rel_tolerance_p = 1.e-4;
constexpr Real pte_rel_tolerance_t = 1.e-4;
constexpr Real pte_abs_tolerance_p = 0.0;
constexpr Real pte_abs_tolerance_t = 0.0;
constexpr int pte_max_iter_per_mat = 16;
constexpr Real line_search_alpha = 1.e-4;
constexpr int line_search_max_iter = 3;
constexpr Real line_search_fac = 0.05;
} // namespace mix_params

namespace mix_impl {
template <typename T,
          typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
constexpr bool isfinite(const T &a) {
  return (a == a) && ((a == 0) || (a != 2 * a));
}

constexpr Real square(const Real x) { return x * x; }

PORTABLE_INLINE_FUNCTION
bool check_nans(Real const *const a, const int n, const bool verbose = false) {
  bool retval = true;
  for (int i = 0; i < n; ++i)
    if (!isfinite(a[i])) {
      retval = false;
#ifndef KOKKOS_ENABLE_CUDA
      if (verbose) {
        printf("bad val in element %i/%i\n", i, n);
      }
#endif // KOKKOS_ENABLE_CUDA
    }
  return retval;
}

template <int n>
PORTABLE_INLINE_FUNCTION bool solve_Ax_b(Real *a, Real *b) {
#ifdef SINGULARITY_USE_KOKKOSKERNELS
#ifndef PORTABILITY_STRATEGY_KOKKOS
#error "Kokkos Kernels requires Kokkos."
#endif
  Real t_[n], w_[n];
  // aliases for kokkos views
  using Unmgd = Kokkos::MemoryTraits<Kokkos::Unmanaged>;
  using Lrgt = Kokkos::LayoutRight;
  using vec_t = Kokkos::View<Real *, Unmgd>;
  // aliases for QR solve template params
  using QR_alg = KokkosBatched::Algo::QR::Unblocked;
  using Lft = KokkosBatched::Side::Left;
  using Trs = KokkosBatched::Trans::Transpose;
  using nTrs = KokkosBatched::Trans::NoTranspose;
  using ApQ_alg = KokkosBatched::Algo::ApplyQ::Unblocked;
  using UP = KokkosBatched::Uplo::Upper;
  using NonU = KokkosBatched::Diag::NonUnit;
  using Tr_alg = KokkosBatched::Algo::Trsv::Unblocked;
  // aliases for solver structs ('invoke' member of the struct is the
  // actual function call)
  using QR_factor = KokkosBatched::SerialQR<QR_alg>;
  using ApplyQ_transpose = KokkosBatched::SerialApplyQ<Lft, Trs, ApQ_alg>;
  using InvertR = KokkosBatched::SerialTrsv<UP, nTrs, NonU, Tr_alg>;
  // view of matrix
  Kokkos::View<Real **, Lrgt, Unmgd> A(a, n, n);
  // view of RHS
  vec_t B(b, n);
  // view of reflectors
  vec_t t(t_, n);
  // view of workspace
  vec_t w(w_, n);
  // QR factor A, A x = B -> Q R x = B
  // store result in A and t
  QR_factor::invoke(A, t, w);
  // Apply Q^T from the left to both sides
  // Q^T Q R x = Q^T B -> R x = Q^T B
  // store result of Q^T B in B
  ApplyQ_transpose::invoke(A, t, B, w);
  // Apply R^-1 from the left to both sides
  // R^-1 R x = R^-1 Q^T B -> x = R^-1 Q^T B
  // store solution vector x in B
  InvertR::invoke(1.0, A, B);
#else
#ifdef PORTABILITY_STRATEGY_KOKKOS
#error "Eigen should not be used with Kokkos."
#endif
  // Eigen VERSION
  Eigen::Map<Eigen::Matrix<Real, n, n, Eigen::RowMajor>> A(a);
  Eigen::Map<Eigen::Matrix<Real, n, 1>> B(b);
  Eigen::Matrix<Real, n, 1> X;
  X = A.lu().solve(B);
  B = X;
#endif // SINGULARITY_USE_KOKKOSKERNELS
  bool retval = check_nans(b, n);
  return retval;
}

template <int nmat, typename EOSIndexer, typename RealIndexer, typename ConstRealIndexer,
          typename LambdaIndexer, typename OFFSETTER>
PORTABLE_INLINE_FUNCTION bool
pte_closure_flag_with_line_impl(EOSIndexer &&eoss, const Real Volume, const Real TotalSIE,
                                ConstRealIndexer &&ComponentMasses,
                                RealIndexer &&ComponentVolumes,
                                RealIndexer &&ComponentEnergies, LambdaIndexer &&lambdas,
                                const OFFSETTER ofst, int &iter) {
  using namespace mix_params;

  Real TotalMass = 0.0, Vsum = 0.0, Esum = 0.0;
  Real Densities[nmat], Pressures[nmat], Temperatures[nmat], dpdr[nmat];
  Real dpde[nmat], dtdr[nmat], dtde[nmat], b[2 * nmat], A[4 * nmat * nmat];
  Real rtemp[nmat], sietemp[nmat];
  // First, normalize the energies and volumes.
  // This may or may not be close to reality since SIE isn't strictly positive.
  for (int i = 0; i < nmat; ++i)
    TotalMass += ComponentMasses[i];
  bool iszero = true;
  for (int i = 0; i < nmat; ++i)
    iszero &= ComponentVolumes[i] <= 1.0e-16;
  if (iszero)
    for (int i = 0; i < nmat; ++i)
      ComponentVolumes[i] = Volume * ComponentMasses[i] / TotalMass;
  for (int i = 0; i < nmat; ++i)
    Vsum += ComponentVolumes[i];
  for (int i = 0; i < nmat; ++i)
    ComponentVolumes[i] *=
        Volume / Vsum; // Are compilers smart enough to calculate this RHS once?
  for (int i = 0; i < nmat; ++i)
    Esum += ComponentEnergies[i] * ComponentMasses[i];
  if (std::abs(Esum) <
      0.0001) { // Then we punt and just distribute based on mass fractions
    for (int i = 0; i < nmat; ++i)
      ComponentEnergies[i] = TotalSIE * TotalMass / ComponentMasses[i];
  } else { // Otherwise, we scale
    for (int i = 0; i < nmat; ++i)
      ComponentEnergies[i] *= (TotalSIE * TotalMass / Esum);
  }

  // This is the main loop
  constexpr const int pte_max_iter = nmat * pte_max_iter_per_mat;
  for (iter = 0; iter < pte_max_iter; ++iter) {
    // Call the EOSs
    for (int mat = 0; mat < nmat; ++mat) {
      Densities[mat] = ComponentMasses[mat] / ComponentVolumes[mat];
      eoss[ofst(mat)].PTofRE(Densities[mat], ComponentEnergies[mat], lambdas[mat],
                             Pressures[mat], Temperatures[mat], dpdr[mat], dpde[mat],
                             dtdr[mat], dtde[mat]);
    }
    // Now, zero the matrix and vector
    for (int i = 0; i < 2 * nmat * 2 * nmat; ++i)
      A[i] = 0.0;
    for (int i = 0; i < 2 * nmat; ++i)
      b[i] = 0.0;
    // print_matrix(2*nmat,2*nmat,A);
    // print_vector(2*nmat,b);
    Real vsum = 0.0;
    Real esum = 0.0;
    Real err = 0.0;
    // build matrix and vector
    for (int n = 0; n < nmat - 1; ++n) {
      const int i = n + nmat;
      // First, do the energy side
      A[n * 2 * nmat + n] = dtde[n];
      A[n * 2 * nmat + n + 1] = -dtde[n + 1];
      A[n * 2 * nmat + i] = dtdr[n];
      A[n * 2 * nmat + i + 1] = -dtdr[n + 1];
      A[(nmat - 1) * 2 * nmat + n] = ComponentMasses[n];
      b[n] = Temperatures[n + 1] - Temperatures[n];
      esum += ComponentMasses[n] * ComponentEnergies[n];
      // Now the density half
      A[i * 2 * nmat + n] = dpde[n];
      A[i * 2 * nmat + n + 1] = -dpde[n + 1];
      A[i * 2 * nmat + i] = dpdr[n];
      A[i * 2 * nmat + i + 1] = -dpdr[n + 1];
      A[(2 * nmat - 1) * 2 * nmat + i] = ComponentMasses[n] / square(Densities[n]);
      b[i] = Pressures[n + 1] - Pressures[n];
      vsum += ComponentMasses[n] / Densities[n];
    }
    A[(nmat - 1) * 2 * nmat + nmat - 1] = ComponentMasses[nmat - 1];
    esum += ComponentMasses[nmat - 1] * ComponentEnergies[nmat - 1];
    b[nmat - 1] = TotalSIE * TotalMass - esum; // esum - TotalSIE*TotalMass;
    A[4 * nmat * nmat - 1] = ComponentMasses[nmat - 1] / square(Densities[nmat - 1]);
    vsum += ComponentMasses[nmat - 1] / Densities[nmat - 1];
    b[2 * nmat - 1] = vsum - Volume;
    for (int i = 0; i < 2 * nmat; ++i)
      err += square(b[i]);
    err *= 0.5;
    if (!solve_Ax_b<2 * nmat>(A, b)) {
      // do something to crash out?  Tell folks what happened?
#ifndef KOKKOS_ENABLE_CUDA
      printf("Crashing out on iteration %i\n", iter);
#endif // KOKKOS_ENABLE_CUDA
      break;
    }
    // LINE SEARCH

    // Now, get overall scaling limit
    Real scale = 1.0;
    for (int m = 0; m < nmat; ++m) {
      Real rt = Densities[m] + scale * b[m + nmat];
      if (rt < 0.0) {
        scale = -0.9 * Densities[m] / b[m + nmat];
      }
      const Real dt = (dtdr[m] * b[m + nmat] + dtde[m] * b[m]);
      const Real tt = Temperatures[m] + scale * dt;
      if (tt < 0.0) scale = -0.9 * Temperatures[m] / dt;
    }
    // Now apply the overall pre-scaling
    for (int i = 0; i < 2 * nmat; ++i)
      b[i] *= scale;
    // Line search
    Real gradfdx = -2.0 * scale * err;
    scale = 1.0;
    Real err_old = err;
    int line_iter = 0;
    Real err_p, err_t;
    do {
      for (int m = 0; m < nmat; ++m) {
        rtemp[m] = Densities[m] + scale * b[m + nmat];
        sietemp[m] = ComponentEnergies[m] + scale * b[m];
        Temperatures[m] = eoss[ofst(m)].TemperatureFromDensityInternalEnergy(
            rtemp[m], sietemp[m], lambdas[m]);
        if (eoss[ofst(m)].PreferredInput() ==
            (thermalqs::density | thermalqs::specific_internal_energy)) {
          Pressures[m] = eoss[ofst(m)].PressureFromDensityInternalEnergy(
              rtemp[m], sietemp[m], lambdas[m]);
        } else if (eoss[ofst(m)].PreferredInput() ==
                   (thermalqs::density | thermalqs::temperature)) {
          Pressures[m] = eoss[ofst(m)].PressureFromDensityTemperature(
              rtemp[m], Temperatures[m], lambdas[m]);
        }
      }
      err_p = 0.0;
      err_t = 0.0;
      for (int n = 0; n < nmat - 1; ++n) {
        err_p += square(Pressures[n + 1] - Pressures[n]);
        err_t += square(Temperatures[n + 1] - Temperatures[n]);
      }
      err = 0.5 * (err_p + err_t);
      line_iter++;
      if (line_iter > line_search_max_iter ||
          err < err_old + line_search_alpha * scale * gradfdx)
        break;
      scale *= line_search_fac;
    } while (true);

    // END LINE SEARCH

    Real mean_p = 0, mean_t = 0;

    for (int i = 0; i < 2 * nmat; ++i)
      b[i] *= scale;
    for (int n = 0; n < nmat; ++n) {
      const int i = nmat + n;
      ComponentEnergies[n] += b[n];
      Densities[n] += b[i];
      ComponentVolumes[n] = ComponentMasses[n] / Densities[n];
      mean_p += ComponentVolumes[n] * Pressures[n];
      mean_t += ComponentMasses[n] * Temperatures[n];
    }
    mean_p /= Volume;
    mean_t /= TotalMass;
    err_p = std::sqrt(err_p);
    err_t = std::sqrt(err_t);
    bool converged_p =
        (err_p < pte_rel_tolerance_p * std::abs(mean_p) || err_p < pte_abs_tolerance_p);
    bool converged_t =
        (err_t < pte_rel_tolerance_t * std::abs(mean_t) || err_t < pte_abs_tolerance_t);
    if (converged_p && converged_t) {
      return true; // FIXME
    }
  }
  return false; // FIXME
}

// RealIndexer types may be different because some might be arrays and
// some might be pointers.
template <int nmat, int res_size, typename T1, typename T2, typename T3>
PORTABLE_INLINE_FUNCTION static void pte_residual(const Real vfrac_tot, const Real utot,
                                                  T1 &&vfrac, Real *u, T2 &&temp,
                                                  T3 &&press, Real *residual) {
  Real vsum = 0.0;
  Real esum = 0.0;
  for (int m = 0; m < nmat; ++m) {
    vsum += vfrac[m];
    esum += u[m];
  }
  residual[0] = vfrac_tot - vsum;
  residual[1] = utot - esum;
  for (int m = 0; m < nmat - 1; ++m) {
    residual[2 + m] = press[m + 1] - press[m];
  }
  for (int m = nmat+1; m < res_size; m++) {
    residual[m] = temp[m - nmat] - temp[m - nmat - 1];
  }
}

template <int nmat, typename RealIndexer>
PORTABLE_INLINE_FUNCTION bool CheckPTE2(const Real rho_total, RealIndexer &&vfrac,
                                        const Real rhobar[nmat], RealIndexer &&press,
                                        Real residual[nmat+1]) {
  using namespace mix_params;
  Real mean_p = vfrac[0] * press[0];
  Real error_p = 0.0;
  for (int m = 1; m < nmat; ++m) {
    mean_p += vfrac[m] * press[m];
    error_p += residual[m + 1] * residual[m + 1];
  }
  error_p = std::sqrt(error_p);
  // Check for convergence
  bool converged_p =
      (error_p < pte_rel_tolerance_p * std::abs(mean_p) || error_p < pte_abs_tolerance_p);
  return converged_p;
}

template <int nmat, typename RealIndexer>
PORTABLE_INLINE_FUNCTION bool CheckPTE(const Real rho_total, RealIndexer &&vfrac,
                                       const Real rhobar[nmat], RealIndexer &&press,
                                       RealIndexer &&temp, Real residual[2 * nmat]) {
  using namespace mix_params;
  Real mean_p = vfrac[0] * press[0];
  Real mean_t = rhobar[0] * temp[0];
  Real error_p = 0.0;
  Real error_t = 0.0;
  for (int m = 1; m < nmat; ++m) {
    mean_p += vfrac[m] * press[m];
    mean_t += rhobar[m] * temp[m];
    error_p += residual[m + 1] * residual[m + 1];
    error_t += residual[m + nmat] * residual[m + nmat];
  }
  mean_t /= rho_total;
  error_p = std::sqrt(error_p);
  error_t = std::sqrt(error_t);
  // Check for convergence
  bool converged_p =
      (error_p < pte_rel_tolerance_p * std::abs(mean_p) || error_p < pte_abs_tolerance_p);
  bool converged_t =
      (error_t < pte_rel_tolerance_t * mean_t || error_t < pte_abs_tolerance_t);
  return (converged_p && converged_t);
}

template <int nmat, typename RealIndexer>
PORTABLE_INLINE_FUNCTION void
get_ideal_pte(const Real vfrac_tot, const Real utot, const Real rho[nmat],
              RealIndexer &&vfrac, RealIndexer &&sie, RealIndexer &&temp,
              RealIndexer &&press, Real &Pequil, Real &Tequil) {
  Real rhoBsum = 0.0;
  Real Asum = 0.0;
  for (int m = 0; m < nmat; m++) {
    // A[m] = vfrac[m] * press[m]/temp[m];
    Asum += vfrac[m] * press[m] / temp[m];
    // B[m] = sie[m]/temp[m];
    rhoBsum += rho[m] * sie[m] / temp[m];
  }

  Tequil = utot / rhoBsum;
  Pequil = Tequil * Asum / vfrac_tot;
}

template <int nmat, int res_size, typename EOSIndexer, typename RealIndexer, typename LambdaIndexer,
          typename OFFSETTER>
PORTABLE_INLINE_FUNCTION void
try_ideal_pte(EOSIndexer &&eos, const Real vfrac_tot, const Real utot,
              const Real rhobar[nmat], RealIndexer &&vfrac, RealIndexer &&sie,
              RealIndexer &&temp, RealIndexer &&press, LambdaIndexer &&lambda,
              const OFFSETTER ofst, Real Cache[nmat][MAX_NUM_LAMBDAS],
              Real residual[res_size]) {
  Real Pequil, Tequil;
  get_ideal_pte<nmat, RealIndexer>(vfrac_tot, utot, rhobar, vfrac, sie, temp, press,
                                   Pequil, Tequil);
  const Real alpha = Pequil / Tequil;

  Real etrial[nmat], vtrial[nmat], ttrial[nmat], ptrial[nmat], utrial[nmat];
  for (int m = 0; m < nmat; m++) {
    etrial[m] = sie[m] / temp[m] * Tequil;
    utrial[m] = rhobar[m] * etrial[m];
    vtrial[m] = vfrac[m] * press[m] / (temp[m] * alpha);
  }
  for (int m = 0; m < nmat; m++) {
    ttrial[m] = eos[ofst(m)].TemperatureFromDensityInternalEnergy(rhobar[m] / vtrial[m],
                                                                  etrial[m], Cache[m]);
    if (eos[ofst(m)].PreferredInput() ==
        (thermalqs::density | thermalqs::specific_internal_energy)) {
      ptrial[m] = eos[ofst(m)].PressureFromDensityInternalEnergy(rhobar[m] / vtrial[m],
                                                                 etrial[m], Cache[m]);
    } else if (eos[ofst(m)].PreferredInput() ==
               (thermalqs::density | thermalqs::temperature)) {
      ptrial[m] = eos[ofst(m)].PressureFromDensityTemperature(rhobar[m] / vtrial[m],
                                                              ttrial[m], Cache[m]);
    }
  }

  Real res_new[res_size];
  pte_residual<nmat, res_size>(vfrac_tot, utot, vtrial, utrial, ttrial, ptrial, res_new);

  Real res0 = 0.0;
  Real res1 = 0.0;
  for (int m = 0; m < res_size; m++) {
    res0 += residual[m] * residual[m];
    res1 += res_new[m] * res_new[m];
  }

  if (res1 < res0) {
    for (int m = 0; m < nmat; m++) {
      press[m] = ptrial[m];
      temp[m] = ttrial[m];
      sie[m] = etrial[m];
      vfrac[m] = vtrial[m];
    }
    for (int m = 0; m < res_size; m++) {
      residual[m] = res_new[m];
    }
  }

  return;
}

template <int size>
PORTABLE_INLINE_FUNCTION int
MatIndex(const int &i, const int &j) {
  return i*size + j;
}

template <int nmat, typename EOSIndexer, typename RealIndexer, typename LambdaIndexer,
          typename OFFSETTER>
PORTABLE_INLINE_FUNCTION bool
pte_closure_josh2_impl(EOSIndexer &&eos, const Real vfrac_tot, const Real sie_tot,
                      RealIndexer &&rho, RealIndexer &&vfrac, RealIndexer &&sie,
                      RealIndexer &&temp, RealIndexer &&press, LambdaIndexer &&lambda,
                      const OFFSETTER ofst, int &niter) {
  using namespace mix_params;
  const Real ilog2 = 1.0 / std::log(2.0);
  Real vsum = 0.0;
  // Normalize vfrac
  for (int m = 0; m < nmat; ++m)
    vsum += vfrac[m];
  for (int m = 0; m < nmat; ++m)
    vfrac[m] /= vsum;

  Real rhobar[nmat]; // This is a fixed quantity: the average density of
                     // material m averaged over the full PTE volume
  for (int m = 0; m < nmat; ++m)
    rhobar[m] = rho[m] * vfrac[m];
  Real rho_total = 0.0;
  for (int m = 0; m < nmat; ++m)
    rho_total += rhobar[m];
  // Renormalize energies as well
  Real utot = rho_total * sie_tot;

  Real Cache[nmat][MAX_NUM_LAMBDAS];
  // set some options and make the initial EOS calls
  enum class EosPreference { RhoT, Rhoe };
  EosPreference eos_choice[nmat];
  for (int m = 0; m < nmat; m++) {
    temp[m] = eos[ofst(m)].TemperatureFromDensityInternalEnergy(rho[m], sie[m], Cache[m]);
    if (eos[ofst(m)].PreferredInput() ==
        (thermalqs::density | thermalqs::specific_internal_energy)) {
      eos_choice[m] = EosPreference::Rhoe;
      press[m] = eos[ofst(m)].PressureFromDensityInternalEnergy(rho[m], sie[m], Cache[m]);
    } else if (eos[ofst(m)].PreferredInput() ==
               (thermalqs::density | thermalqs::temperature)) {
      eos_choice[m] = EosPreference::RhoT;
      press[m] = eos[ofst(m)].PressureFromDensityTemperature(rho[m], temp[m], Cache[m]);
    }
  }

  Real u[nmat];
  Real esum = 0.0;
  for (int m = 0; m < nmat; ++m)
    u[m] = sie[m] * rhobar[m];
  for (int m = 0; m < nmat; ++m)
    esum += u[m];
  for (int m = 0; m < nmat; ++m)
    u[m] *= utot / esum;
  Real jacobian[(nmat+1) * (nmat+1)];
  Real dx[nmat+1];
  Real residual[nmat+1];
  // Real* Cache[nmat];
  // for (int m {0}; m < nmat; ++m) Cache[m] = nullptr;

  pte_residual<nmat,nmat+1>(vfrac_tot, utot, vfrac, u, temp, press, residual);

  // at this point we have initial guesses for rho, vfrac, sie, pressure, temperature
  // try to solve for PTE assuming an ideal gas to reset initial guess
  try_ideal_pte<nmat, nmat+1, EOSIndexer, RealIndexer, LambdaIndexer, OFFSETTER>(
      eos, vfrac_tot, utot, rhobar, vfrac, sie, temp, press, lambda, ofst, Cache,
      residual);
  
  // compute a mass-weighted avg temperature for an initial guess
  Real Tequil = 0.0;
  for (int m = 0; m < nmat; m++) {
    Tequil += rhobar[m] * temp[m];
  }
  Tequil /= rho_total;

  Real vtemp[nmat], rtemp[nmat], utemp[nmat];
  Real dpdT[nmat], dpdv[nmat], dedv[nmat];

  bool converged_p = false;
  bool converged = true;
  niter = 0;
  const int neq = nmat + 1;
  constexpr const int pte_max_iter = nmat * pte_max_iter_per_mat;
  // get the initial residual
  for (niter = 0; niter < pte_max_iter; ++niter) {
    // Calculate errors and break of converged
    converged =
        CheckPTE2<nmat, RealIndexer>(rho_total, vfrac, rhobar, press, residual);
    if (converged) break;

    Real dedT_sum = 0.0;
    for (int m = 0; m < nmat; m++) {
      //////////////////////////////
      // perturb volume fractions
      //////////////////////////////
      const Real deriv_mult = vfrac[m] > 0.01 ? vfrac[m] : 0.01;
      const Real ldv = std::log(derivative_eps * deriv_mult) * ilog2;
      Real dv = std::pow(2.0, std::round(ldv));
      dv *= (vfrac[m] < 0.5 ? 1.0 : -1.0);
      const Real vf_pert = vfrac[m] + dv;
      const Real rho_pert = rhobar[m] / vf_pert;

      Real p_pert = eos[ofst(m)].PressureFromDensityTemperature(rho_pert, Tequil, Cache[m]);
      Real e_pert = eos[ofst(m)].InternalEnergyFromDensityTemperature(rho_pert, Tequil, Cache[m]);
      dpdv[m] = (p_pert - press[m]) / dv;
      dedv[m] = (rhobar[m] * e_pert - u[m]) / dv;
      //////////////////////////////
      // perturb temperature
      //////////////////////////////
      Real dT = Tequil*derivative_eps;

      p_pert = eos[ofst(m)].PressureFromDensityTemperature(rho[m], Tequil + dT, Cache[m]);
      e_pert = eos[ofst(m)].InternalEnergyFromDensityTemperature(rho[m], Tequil + dT, Cache[m]);
      dpdT[m] = (p_pert - press[m]) / dT;
      dedT_sum += (rhobar[m] * e_pert - u[m]) / dT;
    }
    // Fill in the residual
    Real err = 0;
    for (int i = 0; i < neq; ++i)
      err += residual[i] * residual[i];
    err *= 0.5;
    // Fill in the Jacobian
    for (int i = 0; i < neq*neq; ++i)
      jacobian[i] = 0.0;
    for (int m = 0; m < nmat; ++m) {
      jacobian[m] = 1.0;
      jacobian[neq + m] = dedv[m];
    }
    jacobian[neq + nmat] = dedT_sum;
    for (int m = 0; m < nmat - 1; m++) {
      const int ind = MatIndex<nmat+1>(2+m,m);
      jacobian[ind] = dpdv[m];
      jacobian[ind+1] = -dpdv[m+1];
      jacobian[MatIndex<nmat+1>(2+m,nmat)] = dpdT[m] - dpdT[m+1];
    }

    for (int i = 0; i < nmat+1; ++i)
      dx[i] = residual[i];
    if (!solve_Ax_b<nmat+1>(jacobian, dx)) {
      // do something to crash out?  Tell folks what happened?
#ifndef KOKKOS_ENABLE_CUDA
      printf("crashing out at iteration: %i\n", niter);
#endif // KOKKOS_ENABLE_CUDA
      converged = false;
      // std::cout << "Crashing out on iteration "<< niter << std::endl;
      break;
    }
    // Now, get overall scaling limit
    Real scale = 1.0;
    for (int m = 0; m < nmat; ++m) {
      Real vt = vfrac[m] + scale * dx[m];
      if (vt < 0.0) {
        scale = -0.1 * vfrac[m] / dx[m];
      } else if (vt > 1.0) {
        scale = 0.1 * (1.0 - vfrac[m]) / dx[m];
      }
    }
    if (Tequil + scale * dx[nmat] < 0.0) {
      scale = -0.2 * Tequil / dx[nmat];
    }
    // Now apply the overall scaling
    for (int i = 0; i < neq; ++i)
      dx[i] *= scale;
    // Line search
    Real gradfdx = -2.0 * scale * err;
    scale = 1.0;
    Real err_old = err;
    int line_iter = 0;
    Real Ttemp;
    do {
      Ttemp = Tequil + scale * dx[nmat];
      for (int m = 0; m < nmat; ++m) {
        vtemp[m] = vfrac[m] + scale * dx[m];
        rtemp[m] = rhobar[m] / vtemp[m];
        utemp[m] = rhobar[m] * eos[ofst(m)].InternalEnergyFromDensityTemperature(rtemp[m], Ttemp, Cache[m]);
        sie[m] = utemp[m] / rhobar[m];
        temp[m] = Ttemp;
        switch (eos_choice[m]) {
        case EosPreference::Rhoe:
          press[m] =
              eos[ofst(m)].PressureFromDensityInternalEnergy(rtemp[m], sie[m], Cache[m]);
          break;
        case EosPreference::RhoT:
          press[m] =
              eos[ofst(m)].PressureFromDensityTemperature(rtemp[m], temp[m], Cache[m]);
          break;
        }
      }
      pte_residual<nmat,nmat+1>(vfrac_tot, utot, vtemp, utemp, temp, press, residual);
      Real err = 0;
      for (int i = 0; i < neq; ++i)
        err += residual[i] * residual[i];
      err *= 0.5;
      line_iter++;
      if (line_iter > line_search_max_iter ||
          err < err_old + line_search_alpha * scale * gradfdx)
        break;
      scale *= 0.5;
    } while (true);

    // Update values
    Tequil = Ttemp;
    for (int m = 0; m < nmat; ++m) {
      vfrac[m] = vtemp[m];
      rho[m] = rhobar[m] / vfrac[m];
      u[m] = utemp[m];
      sie[m] = u[m] / rhobar[m];
    }
    // niter++;
  } // while (niter < pte_max_iter && !converged);
  for (int m = 0; m < nmat; ++m)
    vfrac[m] *= vfrac_tot;
  return converged;
}


template <int nmat, typename EOSIndexer, typename RealIndexer, typename LambdaIndexer,
          typename OFFSETTER>
PORTABLE_INLINE_FUNCTION bool
pte_closure_josh_impl(EOSIndexer &&eos, const Real vfrac_tot, const Real sie_tot,
                      RealIndexer &&rho, RealIndexer &&vfrac, RealIndexer &&sie,
                      RealIndexer &&temp, RealIndexer &&press, LambdaIndexer &&lambda,
                      const OFFSETTER ofst, int &niter) {
  using namespace mix_params;
  const Real ilog2 = 1.0 / std::log(2.0);
  Real vsum = 0.0;
  // Normalize vfrac
  for (int m = 0; m < nmat; ++m)
    vsum += vfrac[m];
  for (int m = 0; m < nmat; ++m)
    vfrac[m] /= vsum;

  Real rhobar[nmat]; // This is a fixed quantity: the average density of
                     // material m averaged over the full PTE volume
  for (int m = 0; m < nmat; ++m)
    rhobar[m] = rho[m] * vfrac[m];
  Real rho_total = 0.0;
  for (int m = 0; m < nmat; ++m)
    rho_total += rhobar[m];
  // Renormalize energies as well
  Real utot = rho_total * sie_tot;

  Real Cache[nmat][MAX_NUM_LAMBDAS];
  // set some options and make the initial EOS calls
  enum class EosPreference { RhoT, Rhoe };
  EosPreference eos_choice[nmat];
  for (int m = 0; m < nmat; m++) {
    temp[m] = eos[ofst(m)].TemperatureFromDensityInternalEnergy(rho[m], sie[m], Cache[m]);
    if (eos[ofst(m)].PreferredInput() ==
        (thermalqs::density | thermalqs::specific_internal_energy)) {
      eos_choice[m] = EosPreference::Rhoe;
      press[m] = eos[ofst(m)].PressureFromDensityInternalEnergy(rho[m], sie[m], Cache[m]);
    } else if (eos[ofst(m)].PreferredInput() ==
               (thermalqs::density | thermalqs::temperature)) {
      eos_choice[m] = EosPreference::RhoT;
      press[m] = eos[ofst(m)].PressureFromDensityTemperature(rho[m], temp[m], Cache[m]);
    }
  }

  Real u[nmat];
  Real esum = 0.0;
  for (int m = 0; m < nmat; ++m)
    u[m] = sie[m] * rhobar[m];
  for (int m = 0; m < nmat; ++m)
    esum += u[m];
  for (int m = 0; m < nmat; ++m)
    u[m] *= utot / esum;
  Real jacobian[4 * nmat * nmat];
  Real dx[2 * nmat];
  Real residual[2 * nmat];
  // Real* Cache[nmat];
  // for (int m {0}; m < nmat; ++m) Cache[m] = nullptr;

  pte_residual<nmat,2*nmat>(vfrac_tot, utot, vfrac, u, temp, press, residual);

  // at this point we have initial guesses for rho, vfrac, sie, pressure, temperature
  // try to solve for PTE assuming an ideal gas to reset initial guess
  try_ideal_pte<nmat, 2*nmat, EOSIndexer, RealIndexer, LambdaIndexer, OFFSETTER>(
      eos, vfrac_tot, utot, rhobar, vfrac, sie, temp, press, lambda, ofst, Cache,
      residual);

  Real vtemp[nmat], rtemp[nmat], utemp[nmat];
  Real dpde[nmat], dtde[nmat], dpdv[nmat], dtdv[nmat];

  bool converged_p = false;
  bool converged_t = false;
  bool converged = true;
  niter = 0;
  constexpr const int pte_max_iter = nmat * pte_max_iter_per_mat;
  // get the initial residual
  for (niter = 0; niter < pte_max_iter; ++niter) {
    // Calculate errors and break of converged
    converged =
        CheckPTE<nmat, RealIndexer>(rho_total, vfrac, rhobar, press, temp, residual);
    if (converged) break;

    for (int m = 0; m < nmat; m++) {
      //////////////////////////////
      // perturb volume fractions
      //////////////////////////////
      const Real deriv_mult = vfrac[m] > 0.01 ? vfrac[m] : 0.01;
      const Real ldv = std::log(derivative_eps * deriv_mult) * ilog2;
      Real dv = std::pow(2.0, std::round(ldv));
      dv *= (vfrac[m] < 0.5 ? 1.0 : -1.0);
      const Real vf_pert = vfrac[m] + dv;
      const Real rho_pert = rhobar[m] / vf_pert;

      Real p_pert;
      Real t_pert =
          eos[ofst(m)].TemperatureFromDensityInternalEnergy(rho_pert, sie[m], Cache[m]);
      switch (eos_choice[m]) {
      case EosPreference::Rhoe:
        p_pert =
            eos[ofst(m)].PressureFromDensityInternalEnergy(rho_pert, sie[m], Cache[m]);
        break;
      case EosPreference::RhoT:
        p_pert = eos[ofst(m)].PressureFromDensityTemperature(rho_pert, t_pert, Cache[m]);
        break;
      }
      dpdv[m] = (p_pert - press[m]) / dv;
      dtdv[m] = (t_pert - temp[m]) / dv;
      //////////////////////////////
      // perturb energies
      //////////////////////////////
      Real lde = std::log(derivative_eps * std::abs(u[m])) * ilog2;
      const Real de = std::pow(2.0, std::round(lde));
      Real e_pert = (u[m] + de) / rhobar[m];

      t_pert =
          eos[ofst(m)].TemperatureFromDensityInternalEnergy(rho[m], e_pert, Cache[m]);
      switch (eos_choice[m]) {
      case EosPreference::Rhoe:
        p_pert = eos[ofst(m)].PressureFromDensityInternalEnergy(rho[m], e_pert, Cache[m]);
        break;
      case EosPreference::RhoT:
        p_pert = eos[ofst(m)].PressureFromDensityTemperature(rho[m], t_pert, Cache[m]);
        break;
      }
      dpde[m] = (p_pert - press[m]) / de;
      dtde[m] = (t_pert - temp[m]) / de;
      if (std::abs(dtde[m]) < 1.e-16) { // must be on the cold curve
        dtde[m] = derivative_eps;
      }
    }
    // Fill in the residual
    Real err = 0;
    for (int i = 0; i < 2 * nmat; ++i)
      err += residual[i] * residual[i];
    err *= 0.5;
    // Fill in the Jacobian
    for (int i = 0; i < 4 * nmat * nmat; ++i)
      jacobian[i] = 0.0;
    for (int m = 0; m < nmat; ++m) {
      jacobian[m] = 1.0;
      jacobian[2 * nmat + nmat + m] = 1.0;
    }
    jacobian[2 * 2 * nmat] = dpdv[0];
    jacobian[nmat * 2 * nmat + nmat - 1] = -dpdv[nmat - 1];
    jacobian[(nmat + 1) * 2 * nmat] = dtdv[0];
    jacobian[(2 * nmat - 1) * 2 * nmat + nmat - 1] = -dtdv[nmat - 1];
    jacobian[2 * 2 * nmat + nmat] = dpde[0];
    jacobian[nmat * 2 * nmat + 2 * nmat - 1] = -dpde[nmat - 1];
    jacobian[(nmat + 1) * 2 * nmat + nmat] = dtde[0];
    jacobian[(2 * nmat - 1) * 2 * nmat + 2 * nmat - 1] = -dtde[nmat - 1];
    for (int m = 1; m < nmat - 1; ++m) {
      jacobian[(1 + m) * 2 * nmat + m] = -dpdv[m];
      jacobian[(2 + m) * 2 * nmat + m] = dpdv[m];
      jacobian[(nmat + m) * 2 * nmat + m] = -dtdv[m];
      jacobian[(nmat + m + 1) * 2 * nmat + m] = dtdv[m];
      jacobian[(1 + m) * 2 * nmat + nmat + m] = -dpde[m];
      jacobian[(2 + m) * 2 * nmat + nmat + m] = dpde[m];
      jacobian[(nmat + m) * 2 * nmat + nmat + m] = -dtde[m];
      jacobian[(nmat + m + 1) * 2 * nmat + nmat + m] = dtde[m];
    }
    for (int i = 0; i < 2 * nmat; ++i)
      dx[i] = residual[i];
    if (!solve_Ax_b<2 * nmat>(jacobian, dx)) {
      // do something to crash out?  Tell folks what happened?
#ifndef KOKKOS_ENABLE_CUDA
      printf("crashing out at iteration: %i\n", niter);
#endif // KOKKOS_ENABLE_CUDA
      converged = false;
      // std::cout << "Crashing out on iteration "<< niter << std::endl;
      break;
    }
    // Now, get overall scaling limit
    Real scale = 1.0;
    for (int m = 0; m < nmat; ++m) {
      Real vt = vfrac[m] + scale * dx[m];
      if (vt < 0.0) {
        scale = -0.1 * vfrac[m] / dx[m];
      } else if (vt > 1.0) {
        scale = 0.1 * (1.0 - vfrac[m]) / dx[m];
      }
      // maybe the below is dangerous??
      // const Real dt = (dtdv[m] * dx[m] + dtde[m] * dx[m + nmat]);
      // const Real tt = temp[m] + scale * dt;
      // if (tt < 0.0)
      //  scale = -0.1 * temp[m] / dt;
    }
    // Now apply the overall scaling
    for (int i = 0; i < 2 * nmat; ++i)
      dx[i] *= scale;
    // Line search
    Real gradfdx = -2.0 * scale * err;
    scale = 1.0;
    Real err_old = err;
    int line_iter = 0;
    do {
      for (int m = 0; m < nmat; ++m) {
        vtemp[m] = vfrac[m] + scale * dx[m];
        rtemp[m] = rhobar[m] / vtemp[m];
        utemp[m] = u[m] + scale * dx[nmat + m];
        sie[m] = utemp[m] / rhobar[m];
        temp[m] =
            eos[ofst(m)].TemperatureFromDensityInternalEnergy(rtemp[m], sie[m], Cache[m]);
        switch (eos_choice[m]) {
        case EosPreference::Rhoe:
          press[m] =
              eos[ofst(m)].PressureFromDensityInternalEnergy(rtemp[m], sie[m], Cache[m]);
          break;
        case EosPreference::RhoT:
          press[m] =
              eos[ofst(m)].PressureFromDensityTemperature(rtemp[m], temp[m], Cache[m]);
          break;
        }
      }
      pte_residual<nmat,2*nmat>(vfrac_tot, utot, vtemp, utemp, temp, press, residual);
      Real err = 0;
      for (int i = 0; i < 2 * nmat; ++i)
        err += residual[i] * residual[i];
      err *= 0.5;
      line_iter++;
      if (line_iter > line_search_max_iter ||
          err < err_old + line_search_alpha * scale * gradfdx)
        break;
      scale *= 0.5;
    } while (true);

    // Update values
    for (int m = 0; m < nmat; ++m) {
      vfrac[m] = vtemp[m];
      rho[m] = rhobar[m] / vfrac[m];
      u[m] = utemp[m];
      sie[m] = u[m] / rhobar[m];
    }
    // niter++;
  } // while (niter < pte_max_iter && !converged);
  for (int m = 0; m < nmat; ++m)
    vfrac[m] *= vfrac_tot;
  return converged;
}

struct NullPtrIndexer {
  PORTABLE_INLINE_FUNCTION Real *operator[](const int i) { return nullptr; }
};

} // namespace mix_impl

template <int nmat, typename EOSIndexer, typename RealIndexer, typename ConstRealIndexer,
          typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool pte_closure_flag_with_line(
    EOSIndexer &&eoss, const Real Volume, const Real TotalSIE,
    ConstRealIndexer &&ComponentMasses, RealIndexer &&ComponentVolumes,
    RealIndexer &&ComponentEnergies, LambdaIndexer &&lambdas, int &niter) {
  using namespace mix_impl;
  return pte_closure_flag_with_line_impl<nmat>(
      eoss, Volume, TotalSIE, ComponentMasses, ComponentVolumes, ComponentEnergies,
      lambdas, [](const int m) { return m; }, niter);
}

template <int nmat>
PORTABLE_INLINE_FUNCTION bool
pte_closure_flag_with_line_offset(EOS *eoss, const Real Volume, const Real TotalSIE,
                                  int const *const ComponentMats,
                                  const Real *ComponentMasses, Real *ComponentVolumes,
                                  Real *ComponentEnergies, Real **lambdas, int &niter) {
  using namespace mix_impl;
  if (lambdas == nullptr) {
    NullPtrIndexer lambda_indexer;
    return pte_closure_flag_with_line_impl<nmat>(
        eoss, Volume, TotalSIE, ComponentMasses, ComponentVolumes, ComponentEnergies,
        lambda_indexer, [&ComponentMats](const int m) { return ComponentMats[m]; },
        niter);
  }
  return pte_closure_flag_with_line_impl<nmat>(
      eoss, Volume, TotalSIE, ComponentMasses, ComponentVolumes, ComponentEnergies,
      lambdas, [&ComponentMats](const int m) { return ComponentMats[m]; }, niter);
}

template <int nmat, typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool pte_closure_josh2(EOSIndexer &&eos, const Real vfrac_tot,
                                               const Real sie_tot, RealIndexer &&rho,
                                               RealIndexer &&vfrac, RealIndexer &&sie,
                                               RealIndexer &&temp, RealIndexer &&press,
                                               LambdaIndexer &&lambda, int &niter) {
  using namespace mix_impl;
  return pte_closure_josh2_impl<nmat>(
      eos, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press, lambda,
      [&](const int m) { return m; }, niter);
}

template <int nmat>
PORTABLE_INLINE_FUNCTION bool
pte_closure_josh2_offset(EOS *eos, const Real vfrac_tot, const Real sie_tot,
                        int const *const Mats, Real *rho, Real *vfrac, Real *sie,
                        Real *temp, Real *press, Real **lambda, int &niter) {
  using namespace mix_impl;
  if (lambda == nullptr) {
    NullPtrIndexer lambda_indexer;
    return pte_closure_josh2_impl<nmat>(
        eos, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press, lambda_indexer,
        [&](const int m) { return m; }, niter);
  }
  return pte_closure_josh2_impl<nmat>(
      eos, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press, lambda,
      [&Mats](const int m) { return Mats[m]; }, niter);
}
template <int nmat, typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool pte_closure_josh(EOSIndexer &&eos, const Real vfrac_tot,
                                               const Real sie_tot, RealIndexer &&rho,
                                               RealIndexer &&vfrac, RealIndexer &&sie,
                                               RealIndexer &&temp, RealIndexer &&press,
                                               LambdaIndexer &&lambda, int &niter) {
  using namespace mix_impl;
  return pte_closure_josh_impl<nmat>(
      eos, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press, lambda,
      [&](const int m) { return m; }, niter);
}

template <int nmat>
PORTABLE_INLINE_FUNCTION bool
pte_closure_josh_offset(EOS *eos, const Real vfrac_tot, const Real sie_tot,
                        int const *const Mats, Real *rho, Real *vfrac, Real *sie,
                        Real *temp, Real *press, Real **lambda, int &niter) {
  using namespace mix_impl;
  if (lambda == nullptr) {
    NullPtrIndexer lambda_indexer;
    return pte_closure_josh_impl<nmat>(
        eos, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press, lambda_indexer,
        [&](const int m) { return m; }, niter);
  }
  return pte_closure_josh_impl<nmat>(
      eos, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press, lambda,
      [&Mats](const int m) { return Mats[m]; }, niter);
}

template <int nmat, typename EOSIndexer, typename RealIndexer, typename ConstRealIndexer,
          typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool
pte_closure_ideal(EOSIndexer &&eoss, const Real Volume, const Real TotalSIE,
                  ConstRealIndexer &&ComponentMasses, RealIndexer &&ComponentVolumes,
                  RealIndexer &&ComponentEnergies, LambdaIndexer &&lambdas, int &iter) {
  using namespace mix_params;
  using namespace mix_impl;
  Real TotalMass = 0.0, Vsum = 0.0, Esum = 0.0;
  Real Densities[nmat], Pressures[nmat], Temperatures[nmat], dpdr[nmat], dpde[nmat],
      dtdr[nmat], dtde[nmat], dtpdv[nmat], dtpde[nmat], dtdv[nmat], b[2 * nmat],
      A[4 * nmat * nmat], rtemp[nmat], sietemp[nmat];
  // First, normalize the energies and volumes.
  // This may or may not be close to reality since SIE isn't strictly positive.
  for (int i = 0; i < nmat; ++i)
    TotalMass += ComponentMasses[i];
  bool iszero = false;
  for (int i = 0; i < nmat; ++i)
    iszero &= ComponentVolumes[i] <= 1.0e-16;
  if (iszero)
    for (int i = 0; i < nmat; ++i)
      ComponentVolumes[i] = Volume * ComponentMasses[i] / TotalMass;
  for (int i = 0; i < nmat; ++i)
    Vsum += ComponentVolumes[i];
  for (int i = 0; i < nmat; ++i)
    ComponentVolumes[i] *=
        Volume / Vsum; // Are compilers smart enough to calculate this RHS once?
  for (int i = 0; i < nmat; ++i)
    Esum += ComponentEnergies[i] * ComponentMasses[i];
  if (std::abs(Esum) <
      0.0001) { // Then we punt and just distribute based on mass fractions
    for (int i = 0; i < nmat; ++i)
      ComponentEnergies[i] = TotalSIE * TotalMass / ComponentMasses[i];
  } else { // Otherwise, we scale
    for (int i = 0; i < nmat; ++i)
      ComponentEnergies[i] *= (TotalSIE * TotalMass / Esum);
  }

  // This is the main loop
  constexpr const int pte_max_iter = nmat * pte_max_iter_per_mat;
  for (iter = 0; iter < pte_max_iter; ++iter) {
    // Call the EOSs
    for (int mat = 0; mat < nmat; ++mat) {
      Densities[mat] = ComponentMasses[mat] / ComponentVolumes[mat];
      eoss[mat].PTofRE(Densities[mat], ComponentEnergies[mat], lambdas[mat],
                       Pressures[mat], Temperatures[mat], dpdr[mat], dpde[mat], dtdr[mat],
                       dtde[mat]);
      // NOTE: here we are re-using this space for d(P/T)/dr and /de using the
      // chain rule If it is deemed useful, we can output this quantity from the
      // EOS call instead
      dtpdv[mat] = -(dtdr[mat] - Temperatures[mat] / Pressures[mat] * dpdr[mat]) /
                   Pressures[mat] * square(Densities[mat]);
      dtpde[mat] =
          (dtde[mat] - Temperatures[mat] / Pressures[mat] * dpde[mat]) / Pressures[mat];
      dtdv[mat] = -dtdr[mat] * square(Densities[mat]);
    }
    // Now, zero the matrix and vector
    for (int i = 0; i < 2 * nmat * 2 * nmat; ++i)
      A[i] = 0.0;
    for (int i = 0; i < 2 * nmat; ++i)
      b[i] = 0.0;
    Real vsum = 0.0;
    Real esum = 0.0;
    Real err = 0.0;
    // build matrix and vector
    for (int n = 0; n < nmat - 1; ++n) {
      const int i = n + nmat;
      // First, do the energy side
      A[n * 2 * nmat + n] = dtde[n];
      A[n * 2 * nmat + n + 1] = -dtde[n + 1];
      A[n * 2 * nmat + i] = dtdv[n];
      A[n * 2 * nmat + i + 1] = -dtdv[n + 1];
      A[(nmat - 1) * 2 * nmat + n] = ComponentMasses[n];
      b[n] = Temperatures[n + 1] - Temperatures[n];
      // Now the density half
      A[i * 2 * nmat + n] = dtpde[n];
      A[i * 2 * nmat + n + 1] = -dtpde[n + 1];
      A[i * 2 * nmat + i] = dtpdv[n];
      A[i * 2 * nmat + i + 1] = -dtpdv[n + 1];
      A[(2 * nmat - 1) * 2 * nmat + i] = 1;
      b[i] = Pressures[n + 1] / Temperatures[n + 1] - Pressures[n] / Temperatures[n];
      vsum += ComponentMasses[n] / Densities[n];
      esum += ComponentMasses[n] * ComponentEnergies[n];
    }
    A[(nmat - 1) * 2 * nmat + nmat - 1] = ComponentMasses[nmat - 1];
    b[nmat - 1] = esum - TotalMass * TotalSIE;
    A[4 * nmat * nmat - 1] = 1.0;
    vsum += ComponentMasses[nmat - 1] / Densities[nmat - 1];
    b[2 * nmat - 1] = vsum - Volume;
    for (int i = 0; i < 2 * nmat; ++i)
      err += square(b[i]);
    err *= 0.5;
    // Solve the thing
    if (!solve_Ax_b<2 * nmat>(A, b)) {
      // do something to crash out?  Tell folks what happened?
      break;
    }

    // LINE SEARCH

    // Now, get overall scaling limit
    Real scale = 1.0;
    for (int m = 0; m < nmat; ++m) {
      Real rt = ComponentMasses[m] / (ComponentVolumes[m] + scale * b[m + nmat]);
      if (rt < 0.0) {
        scale = -0.9 * ComponentVolumes[m] / b[m + nmat];
      }
      const Real dt = (dtdv[m] * b[m + nmat] + dtde[m] * b[m]);
      const Real tt = Temperatures[m] + scale * dt;
      if (tt < 0.0) scale = -0.9 * Temperatures[m] / dt;
    }
    // Now apply the overall pre-scaling
    for (int i = 0; i < 2 * nmat; ++i)
      b[i] *= scale;
    // Line search
    Real gradfdx = -2.0 * scale * err;
    scale = 1.0;
    Real err_old = err;
    int line_iter = 0;
    Real err_p, err_t;
    do {
      for (int m = 0; m < nmat; ++m) {
        rtemp[m] = ComponentMasses[m] / (ComponentVolumes[m] + scale * b[m + nmat]);
        sietemp[m] = ComponentEnergies[m] + scale * b[m];
        Temperatures[m] = eoss[m].TemperatureFromDensityInternalEnergy(
            rtemp[m], sietemp[m], lambdas[m]);
        if (eoss[m].PreferredInput() ==
            (thermalqs::density | thermalqs::specific_internal_energy)) {
          Pressures[m] =
              eoss[m].PressureFromDensityInternalEnergy(rtemp[m], sietemp[m], lambdas[m]);
        } else if (eoss[m].PreferredInput() ==
                   (thermalqs::density | thermalqs::temperature)) {
          Pressures[m] = eoss[m].PressureFromDensityTemperature(rtemp[m], Temperatures[m],
                                                                lambdas[m]);
        }
      }
      err_p = 0.0;
      err_t = 0.0;
      for (int n = 0; n < nmat - 1; ++n) {
        err_p += square(Pressures[n + 1] - Pressures[n]);
        err_t += square(Temperatures[n + 1] - Temperatures[n]);
      }
      err = 0.5 * (err_p + err_t);
      line_iter++;
      if (line_iter > line_search_max_iter ||
          err < err_old + line_search_alpha * scale * gradfdx)
        break;
      scale *= line_search_fac;
    } while (true);

    // END LINE SEARCH

    Real mean_p = 0, mean_t = 0;

    for (int i = 0; i < 2 * nmat; ++i)
      b[i] *= scale;
    for (int n = 0; n < nmat; ++n) {
      const int i = nmat + n;
      ComponentEnergies[n] += b[n];
      ComponentVolumes[n] += b[i];
      Densities[n] = ComponentMasses[n] / ComponentVolumes[n];
      mean_p += ComponentVolumes[n] * Pressures[n];
      mean_t += ComponentMasses[n] * Temperatures[n];
    }
    mean_p /= Volume;
    mean_t /= TotalMass;
    err_p = std::sqrt(err_p);
    err_t = std::sqrt(err_t);
    bool converged_p =
        (err_p < pte_rel_tolerance_p * std::abs(mean_p) || err_p < pte_abs_tolerance_p);
    bool converged_t =
        (err_t < pte_rel_tolerance_t * std::abs(mean_t) || err_t < pte_abs_tolerance_t);
    if (converged_p && converged_t) {
      return true; // FIXME
    }
  }
  return false; // FIXME
}

// TODO(JMM): The case statement below to switch between templated
// versions of these functions is still extremely gross and I hate it.
// Can we replace with a malloc+free? And if so, what's the
// consequences of that?
template <typename EOSIndexer, typename RealIndexer, typename ConstRealIndexer,
          typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool
pte_closure_flag(int nmat, EOSIndexer &&eoss, const Real Volume, const Real TotalSIE,
                 ConstRealIndexer &&ComponentMasses, RealIndexer &&ComponentVolumes,
                 RealIndexer &&ComponentEnergies, LambdaIndexer &&lambdas, int &niter) {
  switch (nmat) {
  case 1:
    return true; // (or should we call the EOS actually?)
  case 2:
    return pte_closure_flag_with_line<2>(eoss, Volume, TotalSIE, ComponentMasses,
                                         ComponentVolumes, ComponentEnergies, lambdas,
                                         niter);
    break;
  case 3:
    return pte_closure_flag_with_line<3>(eoss, Volume, TotalSIE, ComponentMasses,
                                         ComponentVolumes, ComponentEnergies, lambdas,
                                         niter);
    break;
  case 4:
    return pte_closure_flag_with_line<4>(eoss, Volume, TotalSIE, ComponentMasses,
                                         ComponentVolumes, ComponentEnergies, lambdas,
                                         niter);
    break;
  case 5:
    return pte_closure_flag_with_line<5>(eoss, Volume, TotalSIE, ComponentMasses,
                                         ComponentVolumes, ComponentEnergies, lambdas,
                                         niter);
    break;
  case 6:
    return pte_closure_flag_with_line<6>(eoss, Volume, TotalSIE, ComponentMasses,
                                         ComponentVolumes, ComponentEnergies, lambdas,
                                         niter);
    break;
  case 7:
    return pte_closure_flag_with_line<7>(eoss, Volume, TotalSIE, ComponentMasses,
                                         ComponentVolumes, ComponentEnergies, lambdas,
                                         niter);
    break;
  case 8:
    return pte_closure_flag_with_line<8>(eoss, Volume, TotalSIE, ComponentMasses,
                                         ComponentVolumes, ComponentEnergies, lambdas,
                                         niter);
    break;
  }
  return false;
}
PORTABLE_INLINE_FUNCTION
bool pte_closure_flag_offset(int nmat, EOS *eoss, const Real Volume, const Real TotalSIE,
                             int const *const ComponentMats, const Real *ComponentMasses,
                             Real *ComponentVolumes, Real *ComponentEnergies,
                             Real **lambdas, int &niter) {
  switch (nmat) {
  case 1:
    return true; // (or should we call the EOS actually?)
  case 2:
    return pte_closure_flag_with_line_offset<2>(eoss, Volume, TotalSIE, ComponentMats,
                                                ComponentMasses, ComponentVolumes,
                                                ComponentEnergies, lambdas, niter);
    break;
  case 3:
    return pte_closure_flag_with_line_offset<3>(eoss, Volume, TotalSIE, ComponentMats,
                                                ComponentMasses, ComponentVolumes,
                                                ComponentEnergies, lambdas, niter);
    break;
  case 4:
    return pte_closure_flag_with_line_offset<4>(eoss, Volume, TotalSIE, ComponentMats,
                                                ComponentMasses, ComponentVolumes,
                                                ComponentEnergies, lambdas, niter);
    break;
  case 5:
    return pte_closure_flag_with_line_offset<5>(eoss, Volume, TotalSIE, ComponentMats,
                                                ComponentMasses, ComponentVolumes,
                                                ComponentEnergies, lambdas, niter);
    break;
  case 6:
    return pte_closure_flag_with_line_offset<6>(eoss, Volume, TotalSIE, ComponentMats,
                                                ComponentMasses, ComponentVolumes,
                                                ComponentEnergies, lambdas, niter);
    break;
  case 7:
    return pte_closure_flag_with_line_offset<7>(eoss, Volume, TotalSIE, ComponentMats,
                                                ComponentMasses, ComponentVolumes,
                                                ComponentEnergies, lambdas, niter);
    break;
  case 8:
    return pte_closure_flag_with_line_offset<8>(eoss, Volume, TotalSIE, ComponentMats,
                                                ComponentMasses, ComponentVolumes,
                                                ComponentEnergies, lambdas, niter);
    break;
  }
  return false;
}
template <typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool
pte_closure_josh2(int nmat, EOSIndexer &&eoss, const Real vfrac_tot, const Real sie_tot,
                 RealIndexer &&rho, RealIndexer &&vfrac, RealIndexer &&sie,
                 RealIndexer &&temp, RealIndexer &&press, LambdaIndexer &&lambdas,
                 int &niter) {
  switch (nmat) {
  case 1:
    return true; // (or should we call the EOS actually?)
  case 2:
    return pte_closure_josh2<2>(eoss, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press,
                               lambdas, niter);
    break;
  case 3:
    return pte_closure_josh2<3>(eoss, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press,
                               lambdas, niter);
    break;
  case 4:
    return pte_closure_josh2<4>(eoss, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press,
                               lambdas, niter);
    break;
  case 5:
    return pte_closure_josh2<5>(eoss, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press,
                               lambdas, niter);
    break;
  case 6:
    return pte_closure_josh2<6>(eoss, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press,
                               lambdas, niter);
    break;
  case 7:
    return pte_closure_josh2<7>(eoss, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press,
                               lambdas, niter);
    break;
  case 8:
    return pte_closure_josh2<8>(eoss, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press,
                               lambdas, niter);
    break;
  }
  return false;
}
PORTABLE_INLINE_FUNCTION bool
pte_closure_josh2_offset(int nmat, EOS *eoss, const Real vfrac_tot, const Real sie_tot,
                        int const *const Mats, Real *rho, Real *vfrac, Real *sie,
                        Real *temp, Real *press, Real **lambdas, int &niter) {
  switch (nmat) {
  case 1:
    return true; // (or should we call the EOS actually?)
  case 2:
    return pte_closure_josh2_offset<2>(eoss, vfrac_tot, sie_tot, Mats, rho, vfrac, sie,
                                      temp, press, lambdas, niter);
    break;
  case 3:
    return pte_closure_josh2_offset<3>(eoss, vfrac_tot, sie_tot, Mats, rho, vfrac, sie,
                                      temp, press, lambdas, niter);
    break;
  case 4:
    return pte_closure_josh2_offset<4>(eoss, vfrac_tot, sie_tot, Mats, rho, vfrac, sie,
                                      temp, press, lambdas, niter);
    break;
  case 5:
    return pte_closure_josh2_offset<5>(eoss, vfrac_tot, sie_tot, Mats, rho, vfrac, sie,
                                      temp, press, lambdas, niter);
    break;
  case 6:
    return pte_closure_josh2_offset<6>(eoss, vfrac_tot, sie_tot, Mats, rho, vfrac, sie,
                                      temp, press, lambdas, niter);
    break;
  case 7:
    return pte_closure_josh2_offset<7>(eoss, vfrac_tot, sie_tot, Mats, rho, vfrac, sie,
                                      temp, press, lambdas, niter);
    break;
  case 8:
    return pte_closure_josh2_offset<8>(eoss, vfrac_tot, sie_tot, Mats, rho, vfrac, sie,
                                      temp, press, lambdas, niter);
    break;
  }
  return false;
}
template <typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool
pte_closure_josh(int nmat, EOSIndexer &&eoss, const Real vfrac_tot, const Real sie_tot,
                 RealIndexer &&rho, RealIndexer &&vfrac, RealIndexer &&sie,
                 RealIndexer &&temp, RealIndexer &&press, LambdaIndexer &&lambdas,
                 int &niter) {
  switch (nmat) {
  case 1:
    return true; // (or should we call the EOS actually?)
  case 2:
    return pte_closure_josh<2>(eoss, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press,
                               lambdas, niter);
    break;
  case 3:
    return pte_closure_josh<3>(eoss, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press,
                               lambdas, niter);
    break;
  case 4:
    return pte_closure_josh<4>(eoss, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press,
                               lambdas, niter);
    break;
  case 5:
    return pte_closure_josh<5>(eoss, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press,
                               lambdas, niter);
    break;
  case 6:
    return pte_closure_josh<6>(eoss, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press,
                               lambdas, niter);
    break;
  case 7:
    return pte_closure_josh<7>(eoss, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press,
                               lambdas, niter);
    break;
  case 8:
    return pte_closure_josh<8>(eoss, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press,
                               lambdas, niter);
    break;
  }
  return false;
}
PORTABLE_INLINE_FUNCTION bool
pte_closure_josh_offset(int nmat, EOS *eoss, const Real vfrac_tot, const Real sie_tot,
                        int const *const Mats, Real *rho, Real *vfrac, Real *sie,
                        Real *temp, Real *press, Real **lambdas, int &niter) {
  switch (nmat) {
  case 1:
    return true; // (or should we call the EOS actually?)
  case 2:
    return pte_closure_josh_offset<2>(eoss, vfrac_tot, sie_tot, Mats, rho, vfrac, sie,
                                      temp, press, lambdas, niter);
    break;
  case 3:
    return pte_closure_josh_offset<3>(eoss, vfrac_tot, sie_tot, Mats, rho, vfrac, sie,
                                      temp, press, lambdas, niter);
    break;
  case 4:
    return pte_closure_josh_offset<4>(eoss, vfrac_tot, sie_tot, Mats, rho, vfrac, sie,
                                      temp, press, lambdas, niter);
    break;
  case 5:
    return pte_closure_josh_offset<5>(eoss, vfrac_tot, sie_tot, Mats, rho, vfrac, sie,
                                      temp, press, lambdas, niter);
    break;
  case 6:
    return pte_closure_josh_offset<6>(eoss, vfrac_tot, sie_tot, Mats, rho, vfrac, sie,
                                      temp, press, lambdas, niter);
    break;
  case 7:
    return pte_closure_josh_offset<7>(eoss, vfrac_tot, sie_tot, Mats, rho, vfrac, sie,
                                      temp, press, lambdas, niter);
    break;
  case 8:
    return pte_closure_josh_offset<8>(eoss, vfrac_tot, sie_tot, Mats, rho, vfrac, sie,
                                      temp, press, lambdas, niter);
    break;
  }
  return false;
}
template <typename EOSIndexer, typename RealIndexer, typename ConstRealIndexer,
          typename LambdaIndexer>
PORTABLE_INLINE_FUNCTION bool
pte_closure_ideal(int nmat, EOSIndexer &&eoss, const Real Volume, const Real TotalSIE,
                  ConstRealIndexer &&ComponentMasses, RealIndexer &&ComponentVolumes,
                  RealIndexer &&ComponentEnergies, LambdaIndexer &&lambdas, int &niter) {
  switch (nmat) {
  case 1:
    return true; // (or should we call the EOS actually?)
  case 2:
    return pte_closure_ideal<2>(eoss, Volume, TotalSIE, ComponentMasses, ComponentVolumes,
                                ComponentEnergies, lambdas, niter);
    break;
  case 3:
    return pte_closure_ideal<3>(eoss, Volume, TotalSIE, ComponentMasses, ComponentVolumes,
                                ComponentEnergies, lambdas, niter);
    break;
  case 4:
    return pte_closure_ideal<4>(eoss, Volume, TotalSIE, ComponentMasses, ComponentVolumes,
                                ComponentEnergies, lambdas, niter);
    break;
  case 5:
    return pte_closure_ideal<5>(eoss, Volume, TotalSIE, ComponentMasses, ComponentVolumes,
                                ComponentEnergies, lambdas, niter);
    break;
  case 6:
    return pte_closure_ideal<6>(eoss, Volume, TotalSIE, ComponentMasses, ComponentVolumes,
                                ComponentEnergies, lambdas, niter);
    break;
  case 7:
    return pte_closure_ideal<7>(eoss, Volume, TotalSIE, ComponentMasses, ComponentVolumes,
                                ComponentEnergies, lambdas, niter);
    break;
  case 8:
    return pte_closure_ideal<8>(eoss, Volume, TotalSIE, ComponentMasses, ComponentVolumes,
                                ComponentEnergies, lambdas, niter);
    break;
  }
  return false;
}
} // namespace singularity

#endif // _SINGULARITY_EOS_CLOSURE_MIXED_CELL_MODELS_
