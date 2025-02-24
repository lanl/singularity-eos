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
//------------------------------------------------------------------------------

#ifndef _SINGULARITY_EOS_CLOSURE_MIXED_CELL_MODELS_
#define _SINGULARITY_EOS_CLOSURE_MIXED_CELL_MODELS_

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/fast-math/logs.hpp>
#include <singularity-eos/base/robust_utils.hpp>
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

// ======================================================================
// Implementation details below
// ======================================================================

struct MixParams {
  bool verbose = false;
  Real derivative_eps = 3.0e-6;
  Real pte_rel_tolerance_p = 1.e-6;
  Real pte_rel_tolerance_e = 1.e-6;
  Real pte_rel_tolerance_t = 1.e-4;
  Real pte_abs_tolerance_p = 0.0;
  Real pte_abs_tolerance_e = 1.e-4;
  Real pte_abs_tolerance_t = 0.0;
  Real pte_rel_tolerance_v = 1e-6;
  Real pte_abs_tolerance_v = 1e-6;
  Real pte_residual_tolerance = 1.e-8;
  std::size_t pte_max_iter_per_mat = 128;
  Real line_search_alpha = 1.e-2;
  std::size_t line_search_max_iter = 6;
  Real line_search_fac = 0.5;
  Real vfrac_safety_fac = 0.95;
  Real temperature_limit = 1.0e15;
  Real default_tguess = 300.;
  Real min_dtde = 1.0e-16;
};

struct SolverStatus {
  bool converged = false;
  std::size_t max_niter = 0;
  std::size_t max_line_niter = 0;
  Real residual;
};

namespace mix_impl {
template <typename T,
          typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
constexpr bool isfinite(const T &a) {
  return (a == a) && ((a == 0) || (a != 2 * a));
}

constexpr Real square(const Real x) { return x * x; }

PORTABLE_INLINE_FUNCTION
bool check_nans(Real const *const a, const std::size_t n, const bool verbose = false) {
  bool retval = true;
  for (std::size_t i = 0; i < n; ++i)
    if (!isfinite(a[i])) {
      retval = false;
#ifndef KOKKOS_ENABLE_CUDA
      if (verbose) {
        printf("bad val in element %ld/%ld\n", i, n);
      }
#endif // KOKKOS_ENABLE_CUDA
    }
  return retval;
}

PORTABLE_INLINE_FUNCTION
bool solve_Ax_b_wscr(const std::size_t n, Real *a, Real *b, Real *scr) {
#ifdef SINGULARITY_USE_KOKKOSKERNELS
#ifndef PORTABILITY_STRATEGY_KOKKOS
#error "Kokkos Kernels requires Kokkos."
#endif
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
  vec_t t(scr, n);
  // view of workspace
  vec_t w(scr + n, n);
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
#warning "Eigen should not be used with Kokkos."
#endif
  // Eigen VERSION
  Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> A(a, n,
                                                                                     n);
  Eigen::Map<Eigen::VectorXd> B(b, n);
  Eigen::Map<Eigen::VectorXd> X(scr, n);
  X = A.lu().solve(B);
  B = X;
#endif // SINGULARITY_USE_KOKKOSKERNELS
  bool retval = check_nans(b, n);
  return retval;
}

struct NullPtrIndexer {
  PORTABLE_INLINE_FUNCTION Real *operator[](const std::size_t i) { return nullptr; }
};

class CacheAccessor {
 public:
  CacheAccessor() = default;
  PORTABLE_INLINE_FUNCTION
  explicit CacheAccessor(Real *scr) : cache_(scr) {}
  PORTABLE_FORCEINLINE_FUNCTION
  Real *operator[](const std::size_t m) const { return cache_ + m * MAX_NUM_LAMBDAS; }

 private:
  Real *cache_;
};

// ======================================================================
// Base class
// ======================================================================
template <typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
class PTESolverBase {
 public:
  PTESolverBase() = delete;
  PORTABLE_INLINE_FUNCTION std::size_t Nmat() const { return nmat; }
  PORTABLE_INLINE_FUNCTION std::size_t &Niter() { return niter; }
  // Fixup is meant to be a hook for derived classes to provide arbitrary manipulations
  // after each iteration of the Newton solver.  This version just renormalizes the
  // volume fractions, which is useful to deal with roundoff error.
  PORTABLE_INLINE_FUNCTION
  virtual void Fixup() const {
    Real vsum = 0;
    for (std::size_t m = 0; m < nmat; ++m) {
      vsum += vfrac[m];
    }
    for (std::size_t m = 0; m < nmat; ++m)
      vfrac[m] *= robust::ratio(vfrac_total, vsum);
  }
  // Finalize restores the temperatures, energies, and pressures to unscaled values from
  // the internally scaled quantities used by the solvers
  PORTABLE_INLINE_FUNCTION
  void Finalize() {
    for (std::size_t m = 0; m < nmat; m++) {
      temp[m] *= Tnorm;
      u[m] *= uscale;
      press[m] *= uscale;
    }
  }
  // Solve the linear system for the update dx
  PORTABLE_INLINE_FUNCTION
  bool Solve() const {
    for (std::size_t m = 0; m < neq; m++)
      dx[m] = residual[m];
    return solve_Ax_b_wscr(neq, jacobian, dx, sol_scratch);
  }
  // Parameters for the solver
  PORTABLE_INLINE_FUNCTION
  const MixParams &GetParams() const { return params_; }

 protected:
  PORTABLE_INLINE_FUNCTION
  PTESolverBase(std::size_t nmats, std::size_t neqs, const EOSIndexer &eos_,
                const Real vfrac_tot, const Real sie_tot, const RealIndexer &rho_,
                const RealIndexer &vfrac_, const RealIndexer &sie_,
                const RealIndexer &temp_, const RealIndexer &press_,
                LambdaIndexer &lambda_, Real *&scratch, Real Tnorm,
                const MixParams &params = MixParams())
      : params_(params), nmat(nmats), neq(neqs), niter(0), vfrac_total(vfrac_tot),
        sie_total(sie_tot), eos(eos_), rho(rho_), vfrac(vfrac_), sie(sie_), temp(temp_),
        press(press_), lambda(lambda_), Tnorm(Tnorm) {
    jacobian = AssignIncrement(scratch, neq * neq);
    dx = AssignIncrement(scratch, neq);
    sol_scratch = AssignIncrement(scratch, 2 * neq);
    residual = AssignIncrement(scratch, neq);
    u = AssignIncrement(scratch, nmat);
    rhobar = AssignIncrement(scratch, nmat);
    if (needs_cache_) {
      Cache = CacheAccessor(AssignIncrement(scratch, nmat * MAX_NUM_LAMBDAS));
    }
  }

  PORTABLE_FORCEINLINE_FUNCTION
  void InitRhoBarandRho() {
    // rhobar is a fixed quantity: the average density of
    // material m averaged over the full PTE volume
    rho_total = 0.0;
    for (std::size_t m = 0; m < nmat; ++m) {
      PORTABLE_REQUIRE(vfrac[m] > 0.,
                       "Non-positive volume fraction provided to PTE solver");
      PORTABLE_REQUIRE(rho[m] > 0., "Non-positive density provided to PTE solver");
      rhobar[m] = rho[m] * vfrac[m];
      rho_total += rhobar[m];
    }
  }

  PORTABLE_FORCEINLINE_FUNCTION
  void SetVfracFromT(const Real T) {
    Real vsum = 0.0;
    // set volume fractions
    for (std::size_t m = 0; m < nmat; ++m) {
      const Real rho_min = eos[m].RhoPmin(T);
      const Real vmax = std::min(0.9 * robust::ratio(rhobar[m], rho_min), 1.0);
      vfrac[m] = (vfrac[m] > 0.0 ? std::min(vmax, vfrac[m]) : vmax);
      vsum += vfrac[m];
    }
    // Normalize vfrac
    for (std::size_t m = 0; m < nmat; ++m) {
      vfrac[m] *= robust::ratio(vfrac_total, vsum);
    }
  }

  PORTABLE_FORCEINLINE_FUNCTION
  Real GetTguess() {
    // guess some non-zero temperature to start
    // If a guess was passed in, it's saved in Tnorm
    Real Tguess{};
    if (Tnorm > 0.0) {
      Tguess = Tnorm;
    } else {
      Tguess = params_.default_tguess;
      for (std::size_t m = 0; m < nmat; ++m)
        Tguess = std::max(Tguess, temp[m]);
    }
    PORTABLE_REQUIRE(Tguess > 0., "Non-positive temperature guess for PTE");
    // check for sanity.  basically checks that the input temperatures weren't garbage
    PORTABLE_REQUIRE(Tguess < params_.temperature_limit,
                     "Very large input temperature or temperature guess");
    // iteratively increase temperature guess until all rho's are above rho_at_pmin
    const Real Tfactor = 10.0;
    bool rho_fail;
    for (std::size_t i = 0; i < 3; i++) {
      SetVfracFromT(Tguess);
      // check to make sure the normalization didn't put us below rho_at_pmin
      rho_fail = false;
      for (std::size_t m = 0; m < nmat; ++m) {
        const Real rho_min = eos[m].RhoPmin(Tguess);
        rho[m] = robust::ratio(rhobar[m], vfrac[m]);
        if (rho[m] < rho_min) {
          rho_fail = true;
          Tguess *= Tfactor;
          break;
        }
      }
      if (!rho_fail) break;
    }

    if (rho_fail && params_.verbose) {
      PORTABLE_ALWAYS_WARN(
          "rho < rho_min in PTE initialization!  Solver may not converge.\n");
    }
    return Tguess;
  }

  template <typename EOS_t, typename Indexer_t>
  PORTABLE_FORCEINLINE_FUNCTION static Real
  GetPressureFromPreferred(const EOS_t &eos, const Real rho, const Real T, Real sie,
                           Indexer_t lambda, const bool do_e_lookup) {
    Real P{};
    if (eos.PreferredInput() ==
        (thermalqs::density | thermalqs::specific_internal_energy)) {
      if (do_e_lookup) {
        sie = eos.InternalEnergyFromDensityTemperature(rho, T, lambda);
      }
      P = eos.PressureFromDensityInternalEnergy(rho, sie, lambda);
    } else if (eos.PreferredInput() == (thermalqs::density | thermalqs::temperature)) {
      P = eos.PressureFromDensityTemperature(rho, T, lambda);
    }
    return P;
  }

  // Initialize the volume fractions, avg densities, temperatures, energies, and
  // pressures of the materials.  Compute the total density and internal energy.
  // Scale the energy densities u = rho sie and pressures with the total internal energy
  // density.  This makes Sum(u(mat)) = 1.  Also scale the temperature with the initial
  // guess, so all the initial, scaled material temperatures are = 1.
  PORTABLE_INLINE_FUNCTION
  void InitBase() {

    // intialize rhobar array and final density
    InitRhoBarandRho();
    const Real utotal = rho_total * sie_total;
    uscale = std::abs(utotal);
    // TODO(): Consider edge case when utotal \simeq 0
    utotal_scale = robust::ratio(utotal, uscale);

    // guess some non-zero temperature to start
    const Real Tguess = GetTguess();
    // set the temperature normalization
    Tnorm = Tguess;

    for (std::size_t m = 0; m < nmat; m++) {
      // scaled initial guess for temperature is just 1
      temp[m] = 1.0;
      sie[m] = eos[m].InternalEnergyFromDensityTemperature(rho[m], Tguess, GetLambda(m));
      // note the scaling of pressure
      press[m] = robust::ratio(this->GetPressureFromPreferred(
                                   eos[m], rho[m], Tguess, sie[m], GetLambda(m), false),
                               uscale);
    }

    // note the scaling of the material internal energy densities
    for (std::size_t m = 0; m < nmat; ++m)
      u[m] = sie[m] * robust::ratio(rhobar[m], uscale);
  }

  PORTABLE_INLINE_FUNCTION
  Real ResidualNorm() const {
    Real norm = 0.0;
    for (std::size_t m = 0; m < neq; m++)
      norm += residual[m] * residual[m];
    return 0.5 * norm;
  }

  PORTABLE_FORCEINLINE_FUNCTION
  std::size_t MatIndex(const std::size_t &i, const std::size_t &j) const {
    return i * neq + j;
  }

  // Compute the equilibrium pressure and temperature assuming an ideal EOS for each
  // material.  Set Pideal and Tideal to the ***scaled*** solution.
  PORTABLE_INLINE_FUNCTION
  void GetIdealPTE(Real &Pideal, Real &Tideal) const {
    Real rhoBsum = 0.0;
    Real Asum = 0.0;
    for (std::size_t m = 0; m < nmat; m++) {
      Asum += vfrac[m] * robust::ratio(press[m], temp[m]);
      rhoBsum += rho[m] * vfrac[m] * robust::ratio(sie[m], temp[m]);
    }
    Asum *= uscale / Tnorm;
    rhoBsum /= Tnorm;
    Tideal = uscale / rhoBsum / Tnorm;
    Pideal = robust::ratio(Tnorm * Tideal * Asum, uscale * vfrac_total);
  }

  // Compute the ideal EOS PTE solution and replace the initial guess if it has a lower
  // residual.
  template <typename T>
  PORTABLE_INLINE_FUNCTION void TryIdealPTE(T *solver) {
    Real Pideal, Tideal;
    GetIdealPTE(Pideal, Tideal);

    // temporarily hijack some of the scratch space
    Real *etemp = jacobian;
    Real *ptemp = jacobian + nmat;
    Real *vtemp = jacobian + 2 * nmat;
    Real *rtemp = jacobian + 3 * nmat;
    Real *res = jacobian + 4 * nmat;
    // copy out the initial guess
    for (std::size_t m = 0; m < nmat; ++m) {
      etemp[m] = u[m];
      ptemp[m] = press[m];
      vtemp[m] = vfrac[m];
      rtemp[m] = rho[m];
    }
    Real res_norm_old = 0.0;
    for (std::size_t m = 0; m < neq; ++m) {
      res[m] = residual[m];
      res_norm_old += res[m] * res[m];
    }
    // check if the volume fractions are reasonable
    const Real alpha = robust::ratio(Pideal, Tideal);
    for (std::size_t m = 0; m < nmat; ++m) {
      vfrac[m] *= robust::ratio(press[m], (temp[m] * alpha));
      if (Tnorm * Tideal < 0 ||
          robust::ratio(rhobar[m], vfrac[m]) < eos[m].RhoPmin(Tnorm * Tideal)) {
        // abort because this is putting this material into a bad state
        for (std::size_t n = 0; n <= m; n++)
          vfrac[n] = vtemp[n];
        return;
      }
    }
    // fill in the rest of the state
    for (std::size_t m = 0; m < nmat; ++m) {
      rho[m] = robust::ratio(rhobar[m], vfrac[m]);

      const Real sie_m = eos[m].InternalEnergyFromDensityTemperature(
          rho[m], Tnorm * Tideal, GetLambda(m));
      u[m] = rhobar[m] * robust::ratio(sie_m, uscale);
      press[m] =
          robust::ratio(this->GetPressureFromPreferred(eos[m], rho[m], Tnorm * Tideal,
                                                       sie_m, GetLambda(m), false),
                        uscale);
    }
    // fill in the residual
    solver->Residual();
    Real res_norm_new = 0.0;
    for (std::size_t m = 0; m < neq; m++)
      res_norm_new += residual[m] * residual[m];

    if (res_norm_new > res_norm_old) {
      // this didn't work out so reset everything
      for (std::size_t m = 0; m < nmat; ++m) {
        vfrac[m] = vtemp[m];
        rho[m] = rtemp[m];
        u[m] = etemp[m];
        press[m] = ptemp[m];
      }
      for (std::size_t m = 0; m < neq; ++m) {
        residual[m] = res[m];
      }
    } else {
      // did work, fill in temp and energy density
      for (std::size_t m = 0; m < nmat; ++m) {
        temp[m] = Tideal;
        sie[m] = uscale * robust::ratio(u[m], rhobar[m]);
      }
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real *AssignIncrement(Real *&scratch, const std::size_t size) const {
    Real *p = scratch;
    scratch += size;
    return p;
  }

  PORTABLE_FORCEINLINE_FUNCTION
  auto GetLambda(std::size_t m) const {
    if constexpr (needs_cache_) {
      return Cache[m];
    } else {
      return lambda[m];
    }
  }

  static constexpr const bool needs_cache_ =
      std::is_same<LambdaIndexer, NullIndexer>::value;
  const MixParams params_;
  const std::size_t nmat, neq;
  std::size_t niter;
  const Real vfrac_total, sie_total;
  const EOSIndexer &eos;
  const RealIndexer &rho;
  const RealIndexer &vfrac;
  const RealIndexer &sie;
  const RealIndexer &temp;
  const RealIndexer &press;
  const LambdaIndexer &lambda;
  Real *jacobian, *dx, *sol_scratch, *residual, *u, *rhobar;
  CacheAccessor Cache;
  Real rho_total, uscale, utotal_scale, Tnorm;
};

} // namespace mix_impl

template <typename EOSIndexer, typename RealIndexer, typename LambdaIndexer = NullIndexer>
PORTABLE_INLINE_FUNCTION Real ApproxTemperatureFromRhoMatU(
    const std::size_t nmat, EOSIndexer &&eos, const Real u_tot, RealIndexer &&rho,
    RealIndexer &&vfrac, const Real Tguess = 0.0,
    LambdaIndexer &&lambda = NullIndexer()) {
  // should these be passed in?
  constexpr Real minimum_temperature = 1.e-9;
  constexpr Real maximum_temperature = 1.e9;

  // given material microphysical densities, volume fractions, and a total internal energy
  // density (rho e => erg/cm^3), solve for the temperature that gives the right sum
  // of material energies.  this should only be used for a rough guess since it has a
  // hard coded and fairly loose tolerance
  auto ufunc = [&](const Real T) {
    Real usum = 0.0;
    for (std::size_t m = 0; m < nmat; m++) {
      usum += rho[m] * vfrac[m] *
              eos[m].InternalEnergyFromDensityTemperature(rho[m], T, lambda[m]);
    }
    return usum;
  };

  Real ulo = ufunc(minimum_temperature);
  if (u_tot < ulo) return minimum_temperature;
  Real uhi = ufunc(maximum_temperature);
  if (u_tot > uhi) return maximum_temperature;
  Real lTlo = FastMath::lg(minimum_temperature);
  Real lThi = FastMath::lg(maximum_temperature);
  if (Tguess > minimum_temperature && Tguess < maximum_temperature) {
    const Real ug = ufunc(Tguess);
    if (ug < u_tot) {
      lTlo = FastMath::lg(Tguess);
      ulo = ug;
    } else {
      lThi = FastMath::lg(Tguess);
      uhi = ug;
    }
  }
  std::size_t iter = 0;
  constexpr std::size_t max_iter = 10;
  while (lThi - lTlo > 0.01 && iter < max_iter) {
    // apply bisection which is much better behaved
    // for materials that have a flat sie at low temperatures
    const Real lT = 0.5 * (lTlo + lThi);
    const Real uT = ufunc(FastMath::pow2(lT));
    if (uT < u_tot) {
      lTlo = lT;
      ulo = uT;
    } else {
      lThi = lT;
      uhi = uT;
    }
    iter++;
  }

  const Real alpha = robust::ratio((u_tot - ulo), (uhi - ulo));
  return FastMath::pow2((1.0 - alpha) * lTlo + alpha * lThi);
}

// ======================================================================
// PTE Solver RhoT
// ======================================================================
inline int PTESolverRhoTRequiredScratch(const std::size_t nmat, bool with_cache = true) {
  std::size_t neq = nmat + 1;
  return neq * neq                              // jacobian
         + 4 * neq                              // dx, residual, and sol_scratch
         + 6 * nmat                             // all the nmat sized arrays
         + with_cache * MAX_NUM_LAMBDAS * nmat; // the cache
}
inline size_t PTESolverRhoTRequiredScratchInBytes(const std::size_t nmat,
                                                  bool with_cache = true) {
  return PTESolverRhoTRequiredScratch(nmat, with_cache) * sizeof(Real);
}

template <typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
class PTESolverRhoT
    : public mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer> {
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::InitBase;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::AssignIncrement;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::nmat;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::neq;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::ResidualNorm;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::vfrac_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::sie_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::eos;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::rho;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::vfrac;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::sie;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::temp;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::press;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::rho_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::uscale;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::utotal_scale;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::MatIndex;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::TryIdealPTE;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::jacobian;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::residual;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::dx;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::u;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::rhobar;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::Cache;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::Tnorm;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::params_;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::lambda;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::GetLambda;

 public:
  // template the ctor to get type deduction/universal references prior to c++17
  template <typename EOS_t, typename Real_t, typename Lambda_t>
  PORTABLE_INLINE_FUNCTION
  PTESolverRhoT(const std::size_t nmat, EOS_t &&eos, const Real vfrac_tot,
                const Real sie_tot, Real_t &&rho, Real_t &&vfrac, Real_t &&sie,
                Real_t &&temp, Real_t &&press, Lambda_t &&lambda, Real *scratch,
                const Real Tnorm = 0.0, const MixParams &params = MixParams())
      : mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>(
            nmat, nmat + 1, eos, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press, lambda,
            scratch, Tnorm, params) {
    dpdv = AssignIncrement(scratch, nmat);
    dedv = AssignIncrement(scratch, nmat);
    dpdT = AssignIncrement(scratch, nmat);
    vtemp = AssignIncrement(scratch, nmat);
  }

  PORTABLE_INLINE_FUNCTION
  Real Init() {
    InitBase();
    Residual();
    // Leave this in for now, but comment out because I'm not sure it's a good idea
    // TryIdealPTE(this);
    // Set the current guess for the equilibrium temperature.  Note that this is already
    // scaled.
    Tequil = temp[0];
    return ResidualNorm();
  }

  PORTABLE_INLINE_FUNCTION
  void Residual() const {
    Real vsum = 0.0;
    Real esum = 0.0;
    for (std::size_t m = 0; m < nmat; ++m) {
      vsum += vfrac[m];
      esum += u[m];
    }
    residual[0] = vfrac_total - vsum;
    residual[1] = utotal_scale - esum;
    for (std::size_t m = 0; m < nmat - 1; ++m) {
      residual[2 + m] = press[m + 1] - press[m];
    }
  }

  PORTABLE_INLINE_FUNCTION
  bool CheckPTE() const {
    Real mean_p = vfrac[0] * press[0];
    Real error_p = 0.0;
    for (std::size_t m = 1; m < nmat; ++m) {
      mean_p += vfrac[m] * press[m];
      error_p += residual[m + 1] * residual[m + 1];
    }
    error_p = std::sqrt(error_p);
    Real error_u = std::abs(residual[1]);
    // Check for convergence
    bool converged_p = (error_p < params_.pte_rel_tolerance_p * std::abs(mean_p) ||
                        error_p < params_.pte_abs_tolerance_p);
    bool converged_u =
        (error_u < params_.pte_rel_tolerance_e || error_u < params_.pte_abs_tolerance_e);
    return converged_p && converged_u;
  }

  PORTABLE_INLINE_FUNCTION
  void Jacobian() const {
    Real dedT_sum = 0.0;
    for (std::size_t m = 0; m < nmat; m++) {
      //////////////////////////////
      // perturb volume fractions
      //////////////////////////////
      Real dv = (vfrac[m] < 0.5 ? 1.0 : -1.0) * vfrac[m] * params_.derivative_eps;
      const Real vf_pert = vfrac[m] + dv;
      const Real rho_pert = robust::ratio(rhobar[m], vf_pert);

      Real e_pert = eos[m].InternalEnergyFromDensityTemperature(rho_pert, Tnorm * Tequil,
                                                                GetLambda(m));
      Real p_pert =
          robust::ratio(this->GetPressureFromPreferred(eos[m], rho_pert, Tnorm * Tequil,
                                                       e_pert, GetLambda(m), false),
                        uscale);
      dpdv[m] = robust::ratio((p_pert - press[m]), dv);
      dedv[m] = robust::ratio(rhobar[m] * robust::ratio(e_pert, uscale) - u[m], dv);
      //////////////////////////////
      // perturb temperature
      //////////////////////////////
      Real dT = Tequil * params_.derivative_eps;
      e_pert = eos[m].InternalEnergyFromDensityTemperature(rho[m], Tnorm * (Tequil + dT),
                                                           GetLambda(m));
      p_pert = robust::ratio(this->GetPressureFromPreferred(eos[m], rho[m],
                                                            Tnorm * (Tequil + dT), e_pert,
                                                            GetLambda(m), false),
                             uscale);
      dpdT[m] = robust::ratio((p_pert - press[m]), dT);
      dedT_sum += robust::ratio(rhobar[m] * robust::ratio(e_pert, uscale) - u[m], dT);
    }

    // Fill in the Jacobian
    for (std::size_t i = 0; i < neq * neq; ++i)
      jacobian[i] = 0.0;
    for (std::size_t m = 0; m < nmat; ++m) {
      jacobian[m] = 1.0;
      jacobian[neq + m] = dedv[m];
    }
    jacobian[neq + nmat] = dedT_sum;
    for (std::size_t m = 0; m < nmat - 1; m++) {
      const std::size_t ind = MatIndex(2 + m, m);
      jacobian[ind] = dpdv[m];
      jacobian[ind + 1] = -dpdv[m + 1];
      jacobian[MatIndex(2 + m, nmat)] = dpdT[m] - dpdT[m + 1];
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real ScaleDx() const {
    // Each check reduces the scale further if necessary
    Real scale = 1.0;
    // control how big of a step toward vfrac = 0 is allowed
    for (std::size_t m = 0; m < nmat; ++m) {
      if (scale * dx[m] < -params_.vfrac_safety_fac * vfrac[m]) {
        scale = -params_.vfrac_safety_fac * robust::ratio(vfrac[m], dx[m]);
      }
    }
    const Real Tnew = Tequil + scale * dx[nmat];
    // control how big of a step toward rho = rho(Pmin) is allowed
    for (std::size_t m = 0; m < nmat; m++) {
      const Real rho_min =
          std::max(eos[m].RhoPmin(Tnorm * Tequil), eos[m].RhoPmin(Tnorm * Tnew));
      const Real alpha_max = robust::ratio(rhobar[m], rho_min);
      if (alpha_max < vfrac[m]) {
        // Despite our best efforts, we're already in the unstable regime (i.e.
        // dPdV_T > 0) so we would actually want to *increase* the step instead
        // of decreasing it. As a result, this code doesn't work as intended and
        // should be skipped.
        continue;
      }
      if (scale * dx[m] > 0.5 * (alpha_max - vfrac[m])) {
        scale = robust::ratio(0.5 * (alpha_max - vfrac[m]), dx[m]);
      }
    }
    // control how big of a step toward T = 0 is allowed
    if (scale * dx[nmat] < -0.95 * Tequil) {
      scale = robust::ratio(-0.95 * Tequil, dx[nmat]);
    }
    // Now apply the overall scaling
    for (std::size_t i = 0; i < neq; ++i)
      dx[i] *= scale;
    return scale;
  }

  // Update the solution and return new residual.  Possibly called repeatedly with
  // different scale factors as part of a line search
  PORTABLE_INLINE_FUNCTION
  Real TestUpdate(const Real scale, bool const cache_state = false) {
    if (cache_state) {
      // Store the current state in temp variables for first iteration of line
      // search
      Ttemp = Tequil;
      for (std::size_t m = 0; m < nmat; ++m)
        vtemp[m] = vfrac[m];
    }
    Tequil = Ttemp + scale * dx[nmat];
    for (std::size_t m = 0; m < nmat; ++m) {
      vfrac[m] = vtemp[m] + scale * dx[m];
      rho[m] = robust::ratio(rhobar[m], vfrac[m]);
      u[m] = rhobar[m] * eos[m].InternalEnergyFromDensityTemperature(
                             rho[m], Tnorm * Tequil, GetLambda(m));
      sie[m] = robust::ratio(u[m], rhobar[m]);
      u[m] = robust::ratio(u[m], uscale);
      temp[m] = Tequil;
      press[m] =
          robust::ratio(this->GetPressureFromPreferred(eos[m], rho[m], Tnorm * Tequil,
                                                       sie[m], GetLambda(m), false),
                        uscale);
    }
    Residual();
    return ResidualNorm();
  }

 private:
  Real *dpdv, *dedv, *dpdT, *vtemp;
  Real Tequil, Ttemp;
};

// ======================================================================
// PT space solver
// ======================================================================
inline int PTESolverPTRequiredScratch(const std::size_t nmat, bool with_cache = true) {
  constexpr int neq = 2;
  return neq * neq                              // jacobian
         + 4 * neq                              // dx, residual, and sol_scratch
         + 2 * nmat                             // all the nmat sized arrays
         + with_cache * MAX_NUM_LAMBDAS * nmat; // the cache
}
inline size_t PTESolverPTRequiredScratchInBytes(const std::size_t nmat,
                                                bool with_cache = true) {
  return PTESolverPTRequiredScratch(nmat, with_cache) * sizeof(Real);
}
template <typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
class PTESolverPT
    : public mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer> {
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::InitBase;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::AssignIncrement;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::nmat;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::neq;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::ResidualNorm;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::vfrac_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::sie_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::eos;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::rho;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::vfrac;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::sie;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::temp;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::press;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::rho_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::uscale;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::utotal_scale;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::MatIndex;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::TryIdealPTE;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::jacobian;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::residual;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::dx;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::u;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::rhobar;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::Cache;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::Tnorm;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::params_;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::lambda;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::GetLambda;

  enum RES { RV = 0, RSIE = 1 };

 public:
  // template the ctor to get type deduction/universal references prior to c++17
  template <typename EOS_t, typename Real_t, typename Lambda_t>
  PORTABLE_INLINE_FUNCTION
  PTESolverPT(const std::size_t nmat, EOS_t &&eos, const Real vfrac_tot,
              const Real sie_tot, Real_t &&rho, Real_t &&vfrac, Real_t &&sie,
              Real_t &&temp, Real_t &&press, Lambda_t &&lambda, Real *scratch,
              const Real Tnorm = 0.0, const MixParams &params = MixParams())
      : mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>(
            nmat, 2, eos, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press, lambda,
            scratch, Tnorm, params) {}

  PORTABLE_INLINE_FUNCTION
  Real Init() {
    InitBase();
    Residual();

    // TODO(JMM): Suggestion from Jeff:
    // I suspect that you could multiply the bulk modulus by the
    // volume fraction to produce a weighting that wasn't a bad
    // guess. Worth keeping in mind for the future maybe if the
    // iteration count becomes an issue.
    Pequil = 0;
    Real vsum = 0;
    for (std::size_t m = 0; m < nmat; ++m) {
      // always approach from >0 side
      Pequil += std::abs(press[m]) * vfrac[m];
      vsum += vfrac[m];
    }
    Pequil /= vsum;
    Tequil = 1; // Because it's = Tnorm = initial guess

    // Set the state based on the P/T chosen
    for (std::size_t m = 0; m < nmat; ++m) {
      eos[m].DensityEnergyFromPressureTemperature(Pequil * uscale, Tequil * Tnorm,
                                                  GetLambda(m), rho[m], sie[m]);
      vfrac[m] = robust::ratio(rhobar[m], rho[m]);
      u[m] = robust::ratio(sie[m] * rhobar[m], uscale);
      temp[m] = Tequil;
    }
    Residual();

    // Set the current guess for the equilibrium temperature.  Note
    // that this is already scaled.
    return ResidualNorm();
  }

  PORTABLE_INLINE_FUNCTION
  void Residual() const {
    Real vsum = 0.0;
    Real esum = 0.0;
    for (std::size_t m = 0; m < nmat; ++m) {
      vsum += vfrac[m];
      esum += u[m];
    }
    residual[RV] = vfrac_total - vsum;
    residual[RSIE] = utotal_scale - esum;
  }

  PORTABLE_INLINE_FUNCTION
  bool CheckPTE() const {
    Real error_v = std::abs(residual[RV]);
    Real error_u = std::abs(residual[RSIE]);
    // this may not be quite right as we may need an energy scaling factor
    // Check for convergence
    bool converged_u = (error_u < params_.pte_rel_tolerance_e ||
                        uscale * error_u < params_.pte_abs_tolerance_e);
    bool converged_v =
        (error_v < params_.pte_rel_tolerance_v || error_v < params_.pte_abs_tolerance_v);

    return converged_v && converged_u;
  }

  PORTABLE_INLINE_FUNCTION
  void Jacobian() const {
    // sum_m d u_m / dT )_P
    Real dudT_P_sum = 0.0;
    // sum_m d u_m / dP )_T
    Real dudP_T_sum = 0.0;
    // - sum_m rhobar / rho_m^2 * d rho_m / dT )_P
    Real rbor2_dr_dT_P_sum = 0.0;
    // - sum_m rhobar / rho_m^2 d rho_m / dP )_T
    Real rbor2_dr_dP_T_sum = 0.0;

    // JMM: Note rescaling for u, P, T

    // TODO(JMM): Should we use the thermodynamic derivatives rather
    // than finite differences?
    for (std::size_t m = 0; m < nmat; ++m) {
      Real r_pert = rho[m]; // provide initial guesses
      Real e_pert = sie[m];

      //////////////////////////////
      // perturb pressures
      //////////////////////////////
      Real dp = -Pequil * params_.derivative_eps; // always move towards phase transition
      eos[m].DensityEnergyFromPressureTemperature(uscale * (Pequil + dp), Tnorm * Tequil,
                                                  GetLambda(m), r_pert, e_pert);
      Real drdp = robust::ratio(r_pert - rho[m], dp);
      Real dudp = robust::ratio(robust::ratio(rhobar[m] * e_pert, uscale) - u[m], dp);

      rbor2_dr_dP_T_sum += robust::ratio(rhobar[m], rho[m] * rho[m]) * drdp;
      dudP_T_sum += dudp;

      //////////////////////////////
      // perturb temperatures
      //////////////////////////////
      Real dT = Tequil * params_.derivative_eps;
      eos[m].DensityEnergyFromPressureTemperature(uscale * Pequil, Tnorm * (Tequil + dT),
                                                  GetLambda(m), r_pert, e_pert);
      Real drdT = robust::ratio(r_pert - rho[m], dT);
      Real dudT = robust::ratio(robust::ratio(rhobar[m] * e_pert, uscale) - u[m], dT);

      rbor2_dr_dT_P_sum += robust::ratio(rhobar[m], rho[m] * rho[m]) * drdT;
      dudT_P_sum += dudT;
    }

    // Fill in the Jacobian
    jacobian[0] = -rbor2_dr_dT_P_sum; // TODO(JMM): Check positions
    jacobian[1] = -rbor2_dr_dP_T_sum;
    jacobian[2] = dudT_P_sum;
    jacobian[3] = dudP_T_sum;
  }

  PORTABLE_INLINE_FUNCTION
  Real ScaleDx() const {
    Real scale = 1.0;
    if (scale * dx[0] < -0.95 * Tequil) {
      scale = robust::ratio(-0.95 * Tequil, dx[0]);
    }
    auto bounded = [=](Real Pbound, Real delta) {
      return robust::ratio(robust::ratio(Pbound, uscale) - Pequil, delta);
    };

    for (std::size_t m = 0; m < nmat; ++m) {
      Real Pmin = eos[m].MinimumPressure();
      scale = std::min(std::abs(scale), std::abs(0.95 * bounded(Pmin, dx[1])));
    }
    for (std::size_t m = 0; m < nmat; ++m) {
      Real Ttest = (Tequil + scale * dx[0]) * Tnorm;
      Real Pmax = eos[m].MaximumPressureAtTemperature(Ttest);
      scale = std::min(std::abs(scale), std::abs(0.95 * bounded(Pmax, dx[0])));
    }

    for (std::size_t i = 0; i < neq; ++i) {
      dx[i] *= scale;
    }
    return scale;
  }

  // Update the solution and return new residual.  Possibly called repeatedly with
  // different scale factors as part of a line search
  PORTABLE_INLINE_FUNCTION
  Real TestUpdate(const Real scale, bool const cache_state = false) {
    if (cache_state) {
      // Store the current state in temp variables for first iteration of line
      // search
      Ttemp = Tequil;
      Ptemp = Pequil;
    }
    Tequil = Ttemp + scale * dx[0];
    Pequil = Ptemp + scale * dx[1];
    for (std::size_t m = 0; m < nmat; ++m) {
      eos[m].DensityEnergyFromPressureTemperature(Pequil * uscale, Tequil * Tnorm,
                                                  GetLambda(m), rho[m], sie[m]);
      vfrac[m] = robust::ratio(rhobar[m], rho[m]);
      u[m] = robust::ratio(sie[m] * rhobar[m], uscale);
      temp[m] = Tequil;
      press[m] = Pequil;
    }
    Residual();
    return ResidualNorm();
  }

  // Solve the linear system for the update dx
  // Customized because this system can invert its jacobian analytically
  PORTABLE_INLINE_FUNCTION
  bool Solve() const {
    dx[0] = robust::ratio(jacobian[3] * residual[0] - jacobian[1] * residual[1],
                          jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2]);
    dx[1] = robust::ratio(jacobian[2] * residual[0] - jacobian[0] * residual[1],
                          jacobian[1] * jacobian[2] - jacobian[0] * jacobian[3]);
    return mix_impl::check_nans(dx, 2);
  }

 private:
  // TODO(JMM): Should these have trailing underscores?
  // Current P, T state
  Real Pequil;
  Real Tequil;
  // Scratch states for test update
  Real Ptemp;
  Real Ttemp;
  // TODO(JMM): Should there be a P norm as well as a Tnorm?
};

// ======================================================================
// fixed temperature solver
// ======================================================================
inline std::size_t PTESolverFixedTRequiredScratch(const std::size_t nmat,
                                                  const bool with_cache = true) {
  std::size_t neq = nmat;
  return neq * neq                              // jacobian
         + 4 * neq                              // dx, residual, and sol_scratch
         + 2 * nmat                             // rhobar and u in base
         + 2 * nmat                             // nmat sized arrays in fixed T solver
         + with_cache * MAX_NUM_LAMBDAS * nmat; // the cache
}
inline size_t PTESolverFixedTRequiredScratchInBytes(const std::size_t nmat,
                                                    const bool with_cache = true) {
  return PTESolverFixedTRequiredScratch(nmat, with_cache) * sizeof(Real);
}

template <typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
class PTESolverFixedT
    : public mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer> {
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::AssignIncrement;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::nmat;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::neq;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::ResidualNorm;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::vfrac_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::eos;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::rho;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::vfrac;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::sie;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::temp;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::press;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::rho_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::uscale;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::MatIndex;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::jacobian;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::residual;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::dx;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::u;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::rhobar;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::Cache;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::Tnorm;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::params_;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::lambda;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::GetLambda;

 public:
  // template the ctor to get type deduction/universal references prior to c++17
  // allow the type of the temperature array to be different, potentially a const Real*
  template <typename EOS_t, typename Real_t, typename CReal_t, typename Lambda_t>
  PORTABLE_INLINE_FUNCTION
  PTESolverFixedT(const std::size_t nmat, EOS_t &&eos, const Real vfrac_tot,
                  const Real T_true, Real_t &&rho, Real_t &&vfrac, Real_t &&sie,
                  CReal_t &&temp, Real_t &&press, Lambda_t &lambda, Real *scratch,
                  const MixParams &params = MixParams())
      : mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>(
            nmat, nmat, eos, vfrac_tot, 1.0, rho, vfrac, sie, temp, press, lambda,
            scratch, T_true, params) {
    dpdv = AssignIncrement(scratch, nmat);
    vtemp = AssignIncrement(scratch, nmat);
    Tequil = T_true;
    Ttemp = T_true;
    Tnorm = 1.0;
  }

  PORTABLE_INLINE_FUNCTION
  Real Init() {
    // InitBase();
    // fake init base
    // rhobar is a fixed quantity: the average density of
    // material m averaged over the full PTE volume
    Tnorm = 1.0;
    this->InitRhoBarandRho();
    this->SetVfracFromT(Tequil);
    uscale = 0.0;
    for (std::size_t m = 0; m < nmat; m++) {
      // volume fractions have been potentially reset to ensure densitites are
      // larger than rho(Pmin(Tequil)); set the physical density to reflect
      // this change in volume fraction
      rho[m] = robust::ratio(rhobar[m], vfrac[m]);
      sie[m] = eos[m].InternalEnergyFromDensityTemperature(rho[m], Tequil, GetLambda(m));
      uscale += sie[m] * rho[m];
      // note the scaling of pressure
      press[m] = eos[m].PressureFromDensityTemperature(rho[m], Tequil, GetLambda(m));
    }
    for (std::size_t m = 0; m < nmat; ++m) {
      press[m] = robust::ratio(press[m], uscale);
      u[m] = sie[m] * robust::ratio(rhobar[m], uscale);
    }
    Residual();
    return ResidualNorm();
  }

  PORTABLE_INLINE_FUNCTION
  void Residual() const {
    Real vsum = 0.0;
    for (std::size_t m = 0; m < nmat; ++m) {
      vsum += vfrac[m];
    }
    residual[0] = vfrac_total - vsum;
    for (std::size_t m = 0; m < nmat - 1; ++m) {
      residual[1 + m] = press[m] - press[m + 1];
    }
  }

  PORTABLE_INLINE_FUNCTION
  bool CheckPTE() const {
    Real mean_p = vfrac[0] * press[0];
    Real error_p = 0;
    for (std::size_t m = 1; m < nmat; ++m) {
      mean_p += vfrac[m] * press[m];
      error_p += residual[m + 1] * residual[m + 1];
    }
    error_p = std::sqrt(error_p);
    Real error_v = std::abs(residual[0]);
    // Check for convergence
    bool converged_p = (error_p < params_.pte_rel_tolerance_p * std::abs(mean_p) ||
                        error_p < params_.pte_abs_tolerance_p);
    bool converged_v =
        (error_v < params_.pte_rel_tolerance_e || error_v < params_.pte_abs_tolerance_e);
    return converged_p && converged_v;
  }

  PORTABLE_INLINE_FUNCTION
  void Jacobian() const {
    for (std::size_t m = 0; m < nmat; m++) {
      //////////////////////////////
      // perturb volume fractions
      //////////////////////////////
      Real dv = (vfrac[m] < 0.5 ? 1.0 : -1.0) * vfrac[m] * params_.derivative_eps;
      const Real vf_pert = vfrac[m] + dv;
      const Real rho_pert = robust::ratio(rhobar[m], vf_pert);

      Real p_pert = robust::ratio(
          eos[m].PressureFromDensityTemperature(rho_pert, Tequil, GetLambda(m)), uscale);
      dpdv[m] = robust::ratio((p_pert - press[m]), dv);
    }

    // Fill in the Jacobian
    for (std::size_t i = 0; i < neq * neq; ++i)
      jacobian[i] = 0.0;
    for (std::size_t m = 0; m < nmat; ++m) {
      jacobian[m] = 1.0;
    }
    for (std::size_t m = 0; m < nmat - 1; m++) {
      jacobian[MatIndex(m + 1, m)] = -dpdv[m];
      jacobian[MatIndex(m + 1, m + 1)] = dpdv[m + 1];
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real ScaleDx() const {
    // Each check reduces the scale further if necessary
    Real scale = 1.0;
    // control how big of a step toward vfrac = 0 is allowed
    for (std::size_t m = 0; m < nmat; ++m) {
      if (scale * dx[m] < -params_.vfrac_safety_fac * vfrac[m]) {
        scale = robust::ratio(-params_.vfrac_safety_fac * vfrac[m], dx[m]);
      }
    }
    // control how big of a step toward rho = rho(Pmin) is allowed
    for (std::size_t m = 0; m < nmat; m++) {
      const Real rho_min = eos[m].RhoPmin(Tequil);
      const Real alpha_max = robust::ratio(rhobar[m], rho_min);
      if (alpha_max < vfrac[m]) {
        // Despite our best efforts, we're already in the unstable regime (i.e.
        // dPdV_T > 0) so we would actually want to *increase* the step instead
        // of decreasing it. As a result, this code doesn't work as intended and
        // should be skipped.
        continue;
      }
      if (scale * dx[m] > 0.5 * (alpha_max - vfrac[m])) {
        scale = robust::ratio(0.5 * (alpha_max - vfrac[m]), dx[m]);
      }
    }
    // Now apply the overall scaling
    for (std::size_t i = 0; i < neq; ++i)
      dx[i] *= scale;
    return scale;
  }

  // Update the solution and return new residual.  Possibly called repeatedly with
  // different scale factors as part of a line search
  PORTABLE_INLINE_FUNCTION
  Real TestUpdate(const Real scale, bool const cache_state = false) {
    if (cache_state) {
      // Store the current state in temp variables for first iteration of line
      // search
      for (std::size_t m = 0; m < nmat; ++m)
        vtemp[m] = vfrac[m];
    }
    for (std::size_t m = 0; m < nmat; ++m) {
      vfrac[m] = vtemp[m] + scale * dx[m];
      rho[m] = robust::ratio(rhobar[m], vfrac[m]);
      u[m] = rhobar[m] *
             eos[m].InternalEnergyFromDensityTemperature(rho[m], Tequil, GetLambda(m));
      sie[m] = robust::ratio(u[m], rhobar[m]);
      u[m] = robust::ratio(u[m], uscale);
      press[m] = robust::ratio(
          eos[m].PressureFromDensityTemperature(rho[m], Tequil, GetLambda(m)), uscale);
    }
    Residual();
    return ResidualNorm();
  }

 private:
  Real *dpdv, *vtemp;
  Real Tequil, Ttemp;
};

// ======================================================================
// fixed P solver
// ======================================================================
inline std::size_t PTESolverFixedPRequiredScratch(const std::size_t nmat,
                                                  const bool with_cache = true) {
  std::size_t neq = nmat + 1;
  return neq * neq                              // jacobian
         + 4 * neq                              // dx, residual, and sol_scratch
         + 2 * nmat                             // all the nmat sized arrays in base
         + 3 * nmat                             // all the nmat sized arrays in fixedP
         + with_cache * MAX_NUM_LAMBDAS * nmat; // the cache
}
inline size_t PTESolverFixedPRequiredScratchInBytes(const std::size_t nmat,
                                                    const bool with_cache = true) {
  return PTESolverFixedPRequiredScratch(nmat, with_cache) * sizeof(Real);
}

template <typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
class PTESolverFixedP
    : public mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer> {
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::InitBase;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::AssignIncrement;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::nmat;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::neq;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::ResidualNorm;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::vfrac_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::sie_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::eos;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::rho;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::vfrac;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::sie;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::temp;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::press;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::rho_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::uscale;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::MatIndex;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::TryIdealPTE;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::jacobian;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::residual;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::dx;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::u;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::rhobar;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::Cache;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::Tnorm;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::params_;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::lambda;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::GetLambda;

 public:
  // template the ctor to get type deduction/universal references prior to c++17
  template <typename EOS_t, typename Real_t, typename CReal_t, typename Lambda_t>
  PORTABLE_INLINE_FUNCTION
  PTESolverFixedP(const std::size_t nmat, EOS_t &&eos, const Real vfrac_tot, const Real P,
                  Real_t &&rho, Real_t &&vfrac, Real_t &&sie, Real_t &&temp,
                  CReal_t &&press, Lambda_t &lambda, Real *scratch,
                  const MixParams &params = MixParams())
      : mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>(
            nmat, nmat + 1, eos, vfrac_tot, 1.0, rho, vfrac, sie, temp, press, lambda,
            scratch, 0.0, params) {
    dpdv = AssignIncrement(scratch, nmat);
    dpdT = AssignIncrement(scratch, nmat);
    vtemp = AssignIncrement(scratch, nmat);
    Pequil = P;
  }

  PORTABLE_INLINE_FUNCTION
  Real Init() {
    // InitBase();
    // fake init base
    // rhobar is a fixed quantity: the average density of
    // material m averaged over the full PTE volume
    this->InitRhoBarandRho();

    // guess some non-zero temperature to start
    const Real Tguess = this->GetTguess();
    // set the temperature normalization
    Tnorm = Tguess;
    // calculate u normalization as internal energy guess
    uscale = 0.0;
    for (std::size_t m = 0; m < nmat; m++) {
      // scaled initial guess for temperature is just 1
      temp[m] = 1.0;
      sie[m] = eos[m].InternalEnergyFromDensityTemperature(rho[m], Tguess, GetLambda(m));
      uscale += sie[m] * rho[m];
    }

    // note the scaling of the material internal energy densities
    for (std::size_t m = 0; m < nmat; ++m) {
      u[m] = robust::ratio(sie[m] * rhobar[m], uscale);
      press[m] = robust::ratio(
          eos[m].PressureFromDensityTemperature(rho[m], Tguess, GetLambda(m)), uscale);
    }
    Residual();
    // Set the current guess for the equilibrium temperature.  Note that this is already
    // scaled.
    Tequil = temp[0];
    Ttemp = Tequil;
    return ResidualNorm();
  }

  PORTABLE_INLINE_FUNCTION
  void Residual() const {
    Real vsum = 0.0;
    for (std::size_t m = 0; m < nmat; ++m) {
      vsum += vfrac[m];
      residual[m] = robust::ratio(Pequil, uscale) - press[m];
    }
    residual[nmat] = vfrac_total - vsum;
  }

  PORTABLE_INLINE_FUNCTION
  bool CheckPTE() const {
    Real error_p = 0;
    for (std::size_t m = 0; m < nmat; ++m) {
      error_p += residual[m] * residual[m];
    }
    error_p = std::sqrt(error_p);
    error_p *= uscale;
    Real error_v = std::abs(residual[neq - 1]);
    // Check for convergence
    bool converged_p =
        (error_p < params_.pte_rel_tolerance_p || error_p < params_.pte_abs_tolerance_p);
    bool converged_v =
        (error_v < params_.pte_rel_tolerance_e || error_v < params_.pte_abs_tolerance_e);
    return converged_p && converged_v;
  }

  PORTABLE_INLINE_FUNCTION
  void Jacobian() const {
    for (std::size_t m = 0; m < nmat; m++) {
      //////////////////////////////
      // perturb volume fractions
      //////////////////////////////
      Real dv = (vfrac[m] < 0.5 ? 1.0 : -1.0) * vfrac[m] * params_.derivative_eps;
      const Real vf_pert = vfrac[m] + dv;
      const Real rho_pert = robust::ratio(rhobar[m], vf_pert);

      Real p_pert{};
      Real e_pert{};
      p_pert =
          robust::ratio(this->GetPressureFromPreferred(eos[m], rho_pert, Tnorm * Tequil,
                                                       e_pert, GetLambda(m), true),
                        uscale);
      dpdv[m] = robust::ratio((p_pert - press[m]), dv);
      //////////////////////////////
      // perturb temperature
      //////////////////////////////
      Real dT = Tequil * params_.derivative_eps;

      p_pert = robust::ratio(this->GetPressureFromPreferred(eos[m], rho[m],
                                                            Tnorm * (Tequil + dT), e_pert,
                                                            GetLambda(m), true),
                             uscale);
      dpdT[m] = robust::ratio((p_pert - press[m]), dT);
    }

    // Fill in the Jacobian
    for (std::size_t i = 0; i < neq * neq; ++i)
      jacobian[i] = 0.0;
    for (std::size_t m = 0; m < nmat; m++) {
      jacobian[MatIndex(m, m)] = dpdv[m];
      jacobian[MatIndex(m, nmat)] = dpdT[m];
    }
    for (std::size_t m = 0; m < neq - 1; ++m) {
      jacobian[MatIndex(neq - 1, m)] = 1.0;
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real ScaleDx() const {
    // Each check reduces the scale further if necessary
    Real scale = 1.0;
    // control how big of a step toward vfrac = 0 is allowed
    for (std::size_t m = 0; m < nmat; ++m) {
      if (scale * dx[m] < -params_.vfrac_safety_fac * vfrac[m]) {
        scale = -params_.vfrac_safety_fac * robust::ratio(vfrac[m], dx[m]);
      }
    }
    const Real Tnew = Tequil + scale * dx[nmat];
    // control how big of a step toward rho = rho(Pmin) is allowed
    for (std::size_t m = 0; m < nmat; m++) {
      const Real rho_min =
          std::max(eos[m].RhoPmin(Tnorm * Tequil), eos[m].RhoPmin(Tnorm * Tnew));
      const Real alpha_max = robust::ratio(rhobar[m], rho_min);
      if (alpha_max < vfrac[m]) {
        // Despite our best efforts, we're already in the unstable regime (i.e.
        // dPdV_T > 0) so we would actually want to *increase* the step instead
        // of decreasing it. As a result, this code doesn't work as intended and
        // should be skipped.
        continue;
      }
      if (scale * dx[m] > 0.5 * (alpha_max - vfrac[m])) {
        scale = 0.5 * robust::ratio(alpha_max - vfrac[m], dx[m]);
      }
    }
    // control how big of a step toward T = 0 is allowed
    if (scale * dx[nmat] < -0.95 * Tequil) {
      scale = -0.95 * robust::ratio(Tequil, dx[nmat]);
    }
    // Now apply the overall scaling
    for (std::size_t i = 0; i < neq; ++i)
      dx[i] *= scale;
    return scale;
  }

  // Update the solution and return new residual.  Possibly called repeatedly with
  // different scale factors as part of a line search
  PORTABLE_INLINE_FUNCTION
  Real TestUpdate(const Real scale, bool const cache_state = false) {
    if (cache_state) {
      // Store the current state in temp variables for first iteration of line
      // search
      Ttemp = Tequil;
      for (std::size_t m = 0; m < nmat; ++m)
        vtemp[m] = vfrac[m];
    }
    Tequil = Ttemp + scale * dx[nmat];
    for (std::size_t m = 0; m < nmat; ++m) {
      vfrac[m] = vtemp[m] + scale * dx[m];
      rho[m] = robust::ratio(rhobar[m], vfrac[m]);
      u[m] = rhobar[m] * eos[m].InternalEnergyFromDensityTemperature(
                             rho[m], Tnorm * Tequil, GetLambda(m));
      press[m] = robust::ratio(
          eos[m].PressureFromDensityTemperature(rho[m], Tnorm * Tequil, GetLambda(m)),
          uscale);
      sie[m] = robust::ratio(u[m], rhobar[m]);
      u[m] = robust::ratio(u[m], uscale);
      temp[m] = Tequil;
    }
    Residual();
    return ResidualNorm();
  }

 private:
  Real *dpdv, *dpdT, *vtemp;
  Real Tequil, Ttemp, Pequil;
};

// ======================================================================
// RhoU Solver
// ======================================================================
inline std::size_t PTESolverRhoURequiredScratch(const std::size_t nmat,
                                                const bool with_cache = true) {
  std::size_t neq = 2 * nmat;
  return neq * neq                              // jacobian
         + 4 * neq                              // dx, residual, and sol_scratch
         + 8 * nmat                             // all the nmat sized arrays
         + with_cache * MAX_NUM_LAMBDAS * nmat; // the cache
}
inline size_t PTESolverRhoURequiredScratchInBytes(const std::size_t nmat,
                                                  const bool with_cache = true) {
  return PTESolverRhoURequiredScratch(nmat, with_cache) * sizeof(Real);
}

template <typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
class PTESolverRhoU
    : public mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer> {
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::InitBase;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::AssignIncrement;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::nmat;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::neq;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::ResidualNorm;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::vfrac_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::sie_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::eos;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::rho;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::vfrac;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::sie;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::temp;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::press;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::rho_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::uscale;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::utotal_scale;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::TryIdealPTE;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::MatIndex;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::jacobian;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::residual;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::dx;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::u;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::rhobar;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::Cache;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::Tnorm;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::params_;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::lambda;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>::GetLambda;

 public:
  // template the ctor to get type deduction/universal references prior to c++17
  template <typename EOS_t, typename Real_t, typename Lambda_t>
  PORTABLE_INLINE_FUNCTION
  PTESolverRhoU(const std::size_t nmat, const EOS_t &&eos, const Real vfrac_tot,
                const Real sie_tot, Real_t &&rho, Real_t &&vfrac, Real_t &&sie,
                Real_t &&temp, Real_t &&press, Lambda_t &&lambda, Real *scratch,
                const Real Tnorm = 0.0, const MixParams &params = MixParams())
      : mix_impl::PTESolverBase<EOSIndexer, RealIndexer, LambdaIndexer>(
            nmat, 2 * nmat, eos, vfrac_tot, sie_tot, rho, vfrac, sie, temp, press, lambda,
            scratch, Tnorm, params) {
    dpdv = AssignIncrement(scratch, nmat);
    dtdv = AssignIncrement(scratch, nmat);
    dpde = AssignIncrement(scratch, nmat);
    dtde = AssignIncrement(scratch, nmat);
    vtemp = AssignIncrement(scratch, nmat);
    utemp = AssignIncrement(scratch, nmat);
  }

  PORTABLE_INLINE_FUNCTION
  Real Init() {
    InitBase();
    Residual();
    TryIdealPTE(this);
    return ResidualNorm();
  }

  PORTABLE_INLINE_FUNCTION
  void Residual() const {
    Real vsum = 0.0;
    Real esum = 0.0;
    for (std::size_t m = 0; m < nmat; ++m) {
      vsum += vfrac[m];
      esum += u[m];
    }
    residual[0] = vfrac_total - vsum;
    residual[1] = utotal_scale - esum;
    for (std::size_t m = 0; m < nmat - 1; ++m) {
      residual[2 + m] = press[m + 1] - press[m];
    }
    for (std::size_t m = nmat + 1; m < neq; m++) {
      residual[m] = temp[m - nmat] - temp[m - nmat - 1];
    }
  }

  PORTABLE_INLINE_FUNCTION
  bool CheckPTE() const {
    Real mean_p = vfrac[0] * press[0];
    Real mean_t = rhobar[0] * temp[0];
    Real error_p = 0.0;
    Real error_t = 0.0;
    for (std::size_t m = 1; m < nmat; ++m) {
      mean_p += vfrac[m] * press[m];
      mean_t += rhobar[m] * temp[m];
      error_p += residual[m + 1] * residual[m + 1];
      error_t += residual[m + nmat] * residual[m + nmat];
    }
    mean_t = robust::ratio(mean_t, rho_total);
    error_p = std::sqrt(error_p);
    error_t = std::sqrt(error_t);
    // Check for convergence
    bool converged_p = (error_p < params_.pte_rel_tolerance_p * std::abs(mean_p) ||
                        error_p < params_.pte_abs_tolerance_p);
    bool converged_t = (error_t < params_.pte_rel_tolerance_t * mean_t ||
                        error_t < params_.pte_abs_tolerance_t);
    return (converged_p && converged_t);
  }

  PORTABLE_INLINE_FUNCTION
  void Jacobian() const {
    for (std::size_t m = 0; m < nmat; m++) {
      //////////////////////////////
      // perturb volume fractions
      //////////////////////////////
      Real dv = (vfrac[m] < 0.5 ? 1.0 : -1.0) * vfrac[m] * params_.derivative_eps;
      const Real vf_pert = vfrac[m] + dv;
      const Real rho_pert = robust::ratio(rhobar[m], vf_pert);

      Real p_pert;
      Real t_pert = robust::ratio(
          eos[m].TemperatureFromDensityInternalEnergy(rho_pert, sie[m], GetLambda(m)),
          Tnorm);
      p_pert =
          robust::ratio(this->GetPressureFromPreferred(eos[m], rho_pert, Tnorm * t_pert,
                                                       sie[m], GetLambda(m), false),
                        uscale);
      dpdv[m] = robust::ratio(p_pert - press[m], dv);
      dtdv[m] = robust::ratio(t_pert - temp[m], dv);
      //////////////////////////////
      // perturb energies
      //////////////////////////////
      const Real de = std::abs(u[m]) * params_.derivative_eps;
      Real e_pert = robust::ratio(u[m] + de, rhobar[m]);

      t_pert = robust::ratio(eos[m].TemperatureFromDensityInternalEnergy(
                                 rho[m], uscale * e_pert, GetLambda(m)),
                             Tnorm);
      p_pert = robust::ratio(
          this->GetPressureFromPreferred(eos[m], rho[m], Tnorm * t_pert, uscale * e_pert,
                                         GetLambda(m), false),
          uscale);
      dpde[m] = robust::ratio(p_pert - press[m], de);
      dtde[m] = robust::ratio(t_pert - temp[m], de);
      if (std::abs(dtde[m]) < params_.min_dtde) { // must be on the cold curve
        dtde[m] = params_.derivative_eps;
      }
    }
    for (std::size_t i = 0; i < neq * neq; ++i)
      jacobian[i] = 0.0;
    // TODO(JCD): clean all this up with MatIndex
    for (std::size_t m = 0; m < nmat; ++m) {
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
    for (std::size_t m = 1; m < nmat - 1; ++m) {
      jacobian[(1 + m) * 2 * nmat + m] = -dpdv[m];
      jacobian[(2 + m) * 2 * nmat + m] = dpdv[m];
      jacobian[(nmat + m) * 2 * nmat + m] = -dtdv[m];
      jacobian[(nmat + m + 1) * 2 * nmat + m] = dtdv[m];
      jacobian[(1 + m) * 2 * nmat + nmat + m] = -dpde[m];
      jacobian[(2 + m) * 2 * nmat + nmat + m] = dpde[m];
      jacobian[(nmat + m) * 2 * nmat + nmat + m] = -dtde[m];
      jacobian[(nmat + m + 1) * 2 * nmat + nmat + m] = dtde[m];
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real ScaleDx() const {
    // Each check reduces the scale further if necessary
    Real scale = 1.0;
    for (std::size_t m = 0; m < nmat; ++m) {
      // control how big of a step toward vfrac = 0 is allowed
      if (scale * dx[m] < -0.1 * vfrac[m]) {
        scale = -0.1 * robust::ratio(vfrac[m], dx[m]);
      }
      // try to control steps toward T = 0
      const Real dt = (dtdv[m] * dx[m] + dtde[m] * dx[m + nmat]);
      if (scale * dt < -0.1 * temp[m]) {
        scale = -0.1 * robust::ratio(temp[m], dt);
      }
      const Real tt = temp[m] + scale * dt;
      const Real rho_min =
          std::max(eos[m].RhoPmin(Tnorm * temp[m]), eos[m].RhoPmin(Tnorm * tt));
      const Real alpha_max = robust::ratio(rhobar[m], rho_min);
      // control how big of a step toward rho = rho(Pmin) is allowed
      if (alpha_max < vfrac[m]) {
        // Despite our best efforts, we're already in the unstable regime (i.e.
        // dPdV_T > 0) so we would actually want to *increase* the step instead
        // of decreasing it. As a result, this code doesn't work as intended and
        // should be skipped.
        continue;
      }
      if (scale * dx[m] > 0.5 * (alpha_max - vfrac[m])) {
        scale = 0.5 * robust::ratio(alpha_max - vfrac[m], dx[m]);
      }
    }
    // Now apply the overall scaling
    for (std::size_t i = 0; i < neq; ++i)
      dx[i] *= scale;
    return scale;
  }

  // Update the solution and return new residual.  Possibly called repeatedly with
  // different scale factors as part of a line search
  PORTABLE_INLINE_FUNCTION
  Real TestUpdate(const Real scale, bool const cache_state = false) const {
    if (cache_state) {
      // Store the current state in temp variables for first iteration of line
      // search
      for (std::size_t m = 0; m < nmat; ++m) {
        vtemp[m] = vfrac[m];
        utemp[m] = u[m];
      }
    }
    for (std::size_t m = 0; m < nmat; ++m) {
      vfrac[m] = vtemp[m] + scale * dx[m];
      rho[m] = robust::ratio(rhobar[m], vfrac[m]);
      u[m] = utemp[m] + scale * dx[nmat + m];
      sie[m] = uscale * robust::ratio(u[m], rhobar[m]);
      temp[m] = robust::ratio(
          eos[m].TemperatureFromDensityInternalEnergy(rho[m], sie[m], GetLambda(m)),
          Tnorm);
      press[m] =
          robust::ratio(this->GetPressureFromPreferred(eos[m], rho[m], Tnorm * temp[m],
                                                       sie[m], GetLambda(m), false),
                        uscale);
    }
    Residual();
    return ResidualNorm();
  }

 private:
  Real *dpdv, *dtdv, *dpde, *dtde, *vtemp, *utemp;
};

// ======================================================================
// Solver loop
// ======================================================================
template <class System>
PORTABLE_INLINE_FUNCTION SolverStatus PTESolver(System &s) {
  SolverStatus status;
  Real &err = status.residual;
  bool &converged = status.converged;

  // initialize the system, fill in residual, and get its norm
  err = s.Init();

  // Pull out params
  const MixParams &params = s.GetParams();

  converged = false;
  const std::size_t pte_max_iter = s.Nmat() * params.pte_max_iter_per_mat;
  const Real residual_tol = s.Nmat() * params.pte_residual_tolerance;
  auto &niter = s.Niter();
  for (niter = 0; niter < pte_max_iter; ++niter) {
    status.max_niter = std::max(status.max_niter, niter);

    // Check for convergence
    converged = s.CheckPTE();
    if (converged) break;

    // compute the Jacobian
    s.Jacobian();

    // solve for the Newton step
    bool success = s.Solve();
    if (!success) {
      // do something to crash out?  Tell folks what happened?
      // printf("crashing out at iteration: %ld\n", niter);
      converged = false;
      break;
    }

    // possibly scale the update to stay within reasonable bounds
    Real scale = s.ScaleDx();
    PORTABLE_REQUIRE(scale <= 1.0, "PTE Solver is attempting to increase the step size");

    // Line search
    Real gradfdx = -2.0 * scale * err;
    scale = 1.0; // New scale for line search
    Real err_old = err;
    // Test the update and reset the cache the current state
    err = s.TestUpdate(scale, true /* cache_state */);
    if (err > err_old + params.line_search_alpha * gradfdx) {
      // backtrack to middle of step
      scale = 0.5;
      Real err_mid = s.TestUpdate(scale);
      if (err_mid < err && err_mid < err_old) {
        // We know the half step is better than both the full step and the
        // prior result, so try a pseudo parabolic fit to the error and find its
        // minimum. The `scale` value is bound between 0.75 and 0.25.
        scale = 0.75 + 0.5 * robust::ratio(err_mid - err, err - 2.0 * err_mid + err_old);
      }
      for (std::size_t line_iter = 0; line_iter < params.line_search_max_iter;
           line_iter++) {
        err = s.TestUpdate(scale);
        if (err < err_old + params.line_search_alpha * scale * gradfdx) break;
        // shrink the step if the error isn't reduced enough
        scale *= params.line_search_fac;
        status.max_line_niter = std::max(line_iter, status.max_line_niter);
      }
    }

    // apply fixes post update, e.g. renormalize volume fractions to deal with round-off
    s.Fixup();

    // check for the case where we have converged as much as precision allows
    if (err > 0.5 * err_old && err < residual_tol) {
      converged = true;
      break;
    }
  }
  // Call it converged even though CheckPTE never said it was because the residual is
  // small.  Helps to avoid "failures" where things have actually converged as well as
  // finite precision allows
  if (!converged && err < residual_tol) converged = true;
  // undo any scaling that was applied internally for the solver
  s.Finalize();
  return status;
}

} // namespace singularity

#endif // _SINGULARITY_EOS_CLOSURE_MIXED_CELL_MODELS_
