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
//------------------------------------------------------------------------------

#ifndef _SINGULARITY_EOS_CLOSURE_MIXED_CELL_MODELS_
#define _SINGULARITY_EOS_CLOSURE_MIXED_CELL_MODELS_

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/fast-math/logs.hpp>
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

// TODO(JCD): some of these should be exposed to consumers to allow changes to defaults
namespace mix_params {
constexpr Real derivative_eps = 3.0e-6;
constexpr Real pte_rel_tolerance_p = 1.e-6;
constexpr Real pte_rel_tolerance_e = 1.e-6;
constexpr Real pte_rel_tolerance_t = 1.e-4;
constexpr Real pte_abs_tolerance_p = 0.0;
constexpr Real pte_abs_tolerance_e = 1.e-4;
constexpr Real pte_abs_tolerance_t = 0.0;
constexpr Real pte_residual_tolerance = 1.e-8;
constexpr int pte_max_iter_per_mat = 128;
constexpr Real line_search_alpha = 1.e-2;
constexpr int line_search_max_iter = 6;
constexpr Real line_search_fac = 0.5;
constexpr Real vfrac_safety_fac = 0.95;
constexpr Real minimum_temperature = 1.e-9;
constexpr Real maximum_temperature = 1.e9;
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

PORTABLE_INLINE_FUNCTION
bool solve_Ax_b_wscr(const int n, Real *a, Real *b, Real *scr) {
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
  PORTABLE_INLINE_FUNCTION Real *operator[](const int i) { return nullptr; }
};

class CacheAccessor {
 public:
  CacheAccessor() = default;
  PORTABLE_INLINE_FUNCTION
  explicit CacheAccessor(Real *scr) : cache_(scr) {}
  PORTABLE_FORCEINLINE_FUNCTION
  Real *operator[](const int m) const { return cache_ + m * MAX_NUM_LAMBDAS; }

 private:
  Real *cache_;
};

template <typename EOSIndexer, typename RealIndexer>
class PTESolverBase {
 public:
  PTESolverBase() = delete;
  PORTABLE_INLINE_FUNCTION int Nmat() const { return nmat; }
  PORTABLE_INLINE_FUNCTION int &Niter() { return niter; }
  // Fixup is meant to be a hook for derived classes to provide arbitrary manipulations
  // after each iteration of the Newton solver.  This version just renormalizes the
  // volume fractions, which is useful to deal with roundoff error.
  PORTABLE_INLINE_FUNCTION
  virtual void Fixup() const {
    Real vsum = 0;
    for (int m = 0; m < nmat; ++m)
      vsum += vfrac[m];
    for (int m = 0; m < nmat; ++m)
      vfrac[m] *= vfrac_total / vsum;
  }
  // Finalize restores the temperatures, energies, and pressures to unscaled values from
  // the internally scaled quantities used by the solvers
  PORTABLE_INLINE_FUNCTION
  void Finalize() {
    for (int m = 0; m < nmat; m++) {
      temp[m] *= Tnorm;
      u[m] *= uscale;
      press[m] *= uscale;
    }
  }
  // Solve the linear system for the update dx
  PORTABLE_INLINE_FUNCTION
  bool Solve() const {
    for (int m = 0; m < neq; m++)
      dx[m] = residual[m];
    return solve_Ax_b_wscr(neq, jacobian, dx, sol_scratch);
  }

 protected:
  PORTABLE_INLINE_FUNCTION
  PTESolverBase(int nmats, int neqs, const EOSIndexer &eos_, const Real vfrac_tot,
                const Real sie_tot, const RealIndexer &rho_, const RealIndexer &vfrac_,
                const RealIndexer &sie_, const RealIndexer &temp_,
                const RealIndexer &press_, Real *&scratch, Real Tguess)
      : nmat(nmats), neq(neqs), niter(0), eos(eos_), vfrac_total(vfrac_tot),
        sie_total(sie_tot), rho(rho_), vfrac(vfrac_), sie(sie_), temp(temp_),
        press(press_), Tnorm(Tguess) {
    jacobian = AssignIncrement(scratch, neq * neq);
    dx = AssignIncrement(scratch, neq);
    sol_scratch = AssignIncrement(scratch, 2 * neq);
    residual = AssignIncrement(scratch, neq);
    u = AssignIncrement(scratch, nmat);
    rhobar = AssignIncrement(scratch, nmat);
    Cache = CacheAccessor(AssignIncrement(scratch, nmat * MAX_NUM_LAMBDAS));
  }

  PORTABLE_INLINE_FUNCTION
  void InitRhoBarandRho() {
    // rhobar is a fixed quantity: the average density of
    // material m averaged over the full PTE volume
    rho_total = 0.0;
    for (int m = 0; m < nmat; ++m) {
      rhobar[m] = rho[m] * vfrac[m];
      rho_total += rhobar[m];
    }
  }

  PORTABLE_INLINE_FUNCTION
  void SetVfracFromT(const Real T) {
    Real vsum = 0.0;
    // set volume fractions
    for (int m = 0; m < nmat; ++m) {
      const Real rho_min = eos[m].RhoPmin(T);
      const Real vmax = std::min(0.9 * rhobar[m] / rho_min, 1.0);
      vfrac[m] = (vfrac[m] > 0.0 ? std::min(vmax, vfrac[m]) : vmax);
      vsum += vfrac[m];
    }
    // Normalize vfrac
    for (int m = 0; m < nmat; ++m) {
      vfrac[m] *= vfrac_total / vsum;
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real GetTguess() {
    // guess some non-zero temperature to start
    // If a guess was passed in, it's saved in Tnorm
    Real Tguess{};
    if (Tnorm > 0.0) {
      Tguess = Tnorm;
    } else {
      Tguess = 300.0;
      for (int m = 0; m < nmat; ++m)
        Tguess = std::max(Tguess, temp[m]);
    }
    // check for sanity.  basically checks that the input temperatures weren't garbage
    assert(Tguess < 1.0e15);
    // iteratively increase temperature guess until all rho's are above rho_at_pmin
    const Real Tfactor = 10.0;
    bool rho_fail;
    for (int i = 0; i < 3; i++) {
      SetVfracFromT(Tguess);
      // check to make sure the normalization didn't put us below rho_at_pmin
      rho_fail = false;
      for (int m = 0; m < nmat; ++m) {
        const Real rho_min = eos[m].RhoPmin(Tguess);
        rho[m] = rhobar[m] / vfrac[m];
        if (rho[m] < rho_min) {
          rho_fail = true;
          Tguess *= Tfactor;
          break;
        }
      }
      if (!rho_fail) break;
    }

    if (rho_fail) {
      printf("rho < rho_min in PTE initialization!  Solver may not converge.");
    }
    return Tguess;
  }

  PORTABLE_INLINE_FUNCTION
  Real GetPressureFromPreferred(const EOS &eos, const Real rho, const Real T, Real sie,
                                Real *lambda, const bool do_e_lookup) const {
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
    uscale = rho_total * sie_total;

    // guess some non-zero temperature to start
    const Real Tguess = GetTguess();
    // set the temperature normalization
    Tnorm = Tguess;

    for (int m = 0; m < nmat; m++) {
      // scaled initial guess for temperature is just 1
      temp[m] = 1.0;
      sie[m] = eos[m].InternalEnergyFromDensityTemperature(rho[m], Tguess, Cache[m]);
      // note the scaling of pressure
      press[m] = this->GetPressureFromPreferred(eos[m], rho[m], Tguess, sie[m], Cache[m],
                                                false) /
                 uscale;
    }

    // note the scaling of the material internal energy densities
    for (int m = 0; m < nmat; ++m)
      u[m] = sie[m] * rhobar[m] / uscale;
  }

  PORTABLE_INLINE_FUNCTION
  Real ResidualNorm() const {
    Real norm = 0.0;
    for (int m = 0; m < neq; m++)
      norm += residual[m] * residual[m];
    return 0.5 * norm;
  }

  PORTABLE_FORCEINLINE_FUNCTION
  int MatIndex(const int &i, const int &j) const { return i * neq + j; }

  // Compute the equilibrium pressure and temperature assuming an ideal EOS for each
  // material.  Set Pideal and Tideal to the ***scaled*** solution.
  PORTABLE_INLINE_FUNCTION
  void GetIdealPTE(Real &Pideal, Real &Tideal) const {
    Real rhoBsum = 0.0;
    Real Asum = 0.0;
    for (int m = 0; m < nmat; m++) {
      Asum += vfrac[m] * press[m] / temp[m];
      rhoBsum += rho[m] * vfrac[m] * sie[m] / temp[m];
    }
    Asum *= uscale / Tnorm;
    rhoBsum /= Tnorm;
    Tideal = uscale / rhoBsum / Tnorm;
    Pideal = Tnorm * Tideal * Asum / vfrac_total / uscale;
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
    for (int m = 0; m < nmat; ++m) {
      etemp[m] = u[m];
      ptemp[m] = press[m];
      vtemp[m] = vfrac[m];
      rtemp[m] = rho[m];
    }
    Real res_norm_old = 0.0;
    for (int m = 0; m < neq; ++m) {
      res[m] = residual[m];
      res_norm_old += res[m] * res[m];
    }
    // check if the volume fractions are reasonable
    const Real alpha = Pideal / Tideal;
    for (int m = 0; m < nmat; ++m) {
      vfrac[m] *= press[m] / (temp[m] * alpha);
      if (rhobar[m] / vfrac[m] < eos[m].RhoPmin(Tnorm * Tideal)) {
        // abort because this is putting this material into a bad state
        for (int n = m; n >= 0; n--)
          vfrac[n] = vtemp[n];
        return;
      }
    }
    // fill in the rest of the state
    for (int m = 0; m < nmat; ++m) {
      rho[m] = rhobar[m] / vfrac[m];
      const Real sie_m =
          eos[m].InternalEnergyFromDensityTemperature(rho[m], Tnorm * Tideal, Cache[m]);
      u[m] = rhobar[m] * sie_m / uscale;
      press[m] = this->GetPressureFromPreferred(eos[m], rho[m], Tnorm * Tideal, sie_m,
                                                Cache[m], false) /
                 uscale;
    }
    // fill in the residual
    solver->Residual();
    Real res_norm_new = 0.0;
    for (int m = 0; m < neq; m++)
      res_norm_new += residual[m] * residual[m];

    if (res_norm_new > res_norm_old) {
      // this didn't work out so reset everything
      for (int m = 0; m < nmat; ++m) {
        vfrac[m] = vtemp[m];
        rho[m] = rtemp[m];
        u[m] = etemp[m];
        press[m] = ptemp[m];
      }
      for (int m = 0; m < neq; ++m) {
        residual[m] = res[m];
      }
    } else {
      // did work, fill in temp and energy density
      for (int m = 0; m < nmat; ++m) {
        temp[m] = Tideal;
        sie[m] = uscale * u[m] / rhobar[m];
      }
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real *AssignIncrement(Real *&scratch, const int size) const {
    Real *p = scratch;
    scratch += size;
    return p;
  }

  const int nmat, neq;
  int niter;
  const Real vfrac_total, sie_total;
  const EOSIndexer &eos;
  const RealIndexer &rho;
  const RealIndexer &vfrac;
  const RealIndexer &sie;
  const RealIndexer &temp;
  const RealIndexer &press;
  Real *jacobian, *dx, *sol_scratch, *residual, *u, *rhobar;
  CacheAccessor Cache;
  Real rho_total, uscale, Tnorm;
};

} // namespace mix_impl

template <typename EOSIndexer, typename RealIndexer>
PORTABLE_INLINE_FUNCTION Real ApproxTemperatureFromRhoMatU(
    const int nmat, EOSIndexer &&eos, const Real u_tot, RealIndexer &&rho,
    RealIndexer &&vfrac, const Real Tguess = 0.0) {
  // given material microphysical densities, volume fractions, and a total internal energy
  // density (rho e => erg/cm^3), solve for the temperature that gives the right sum
  // of material energies.  this should only be used for a rough guess since it has a
  // hard coded and fairly loose tolerance
  auto ufunc = [&](const Real T) {
    Real usum = 0.0;
    for (int m = 0; m < nmat; m++) {
      usum += rho[m] * vfrac[m] * eos[m].InternalEnergyFromDensityTemperature(rho[m], T);
    }
    return usum;
  };

  Real ulo = ufunc(mix_params::minimum_temperature);
  if (u_tot < ulo) return mix_params::minimum_temperature;
  Real uhi = ufunc(mix_params::maximum_temperature);
  if (u_tot > uhi) return mix_params::maximum_temperature;
  Real lTlo = FastMath::lg(mix_params::minimum_temperature);
  Real lThi = FastMath::lg(mix_params::maximum_temperature);
  if (Tguess > mix_params::minimum_temperature &&
      Tguess < mix_params::maximum_temperature) {
    const Real ug = ufunc(Tguess);
    if (ug < u_tot) {
      lTlo = FastMath::lg(Tguess);
      ulo = ug;
    } else {
      lThi = FastMath::lg(Tguess);
      uhi = ug;
    }
  }
  int iter = 0;
  constexpr int max_iter = 10;
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

  const Real alpha = (u_tot - ulo) / (uhi - ulo);
  return FastMath::pow2((1.0 - alpha) * lTlo + alpha * lThi);
}

inline int PTESolverRhoTRequiredScratch(const int nmat) {
  int neq = nmat + 1;
  return neq * neq                 // jacobian
         + 4 * neq                 // dx, residual, and sol_scratch
         + 6 * nmat                // all the nmat sized arrays
         + MAX_NUM_LAMBDAS * nmat; // the cache
}
inline size_t PTESolverRhoTRequiredScratchInBytes(const int nmat) {
  return PTESolverRhoTRequiredScratch(nmat) * sizeof(Real);
}

template <typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
class PTESolverRhoT : public mix_impl::PTESolverBase<EOSIndexer, RealIndexer> {
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::InitBase;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::AssignIncrement;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::nmat;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::neq;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::ResidualNorm;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::vfrac_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::sie_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::eos;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::rho;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::vfrac;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::sie;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::temp;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::press;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::rho_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::uscale;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::MatIndex;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::TryIdealPTE;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::jacobian;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::residual;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::dx;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::u;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::rhobar;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::Cache;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::Tnorm;

 public:
  // template the ctor to get type deduction/universal references prior to c++17
  template <typename EOS_t, typename Real_t, typename Lambda_t>
  PORTABLE_INLINE_FUNCTION
  PTESolverRhoT(const int nmat, EOS_t &&eos, const Real vfrac_tot, const Real sie_tot,
                Real_t &&rho, Real_t &&vfrac, Real_t &&sie, Real_t &&temp, Real_t &&press,
                Lambda_t &&lambda, Real *scratch, const Real Tguess = 0.0)
      : mix_impl::PTESolverBase<EOSIndexer, RealIndexer>(nmat, nmat + 1, eos, vfrac_tot,
                                                         sie_tot, rho, vfrac, sie, temp,
                                                         press, scratch, Tguess) {
    dpdv = AssignIncrement(scratch, nmat);
    dedv = AssignIncrement(scratch, nmat);
    dpdT = AssignIncrement(scratch, nmat);
    vtemp = AssignIncrement(scratch, nmat);
    // TODO(JCD): use whatever lambdas are passed in
    /*for (int m = 0; m < nmat; m++) {
      if (lambda[m] != nullptr) Cache[m] = lambda[m];
    }*/
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
    for (int m = 0; m < nmat; ++m) {
      vsum += vfrac[m];
      esum += u[m];
    }
    residual[0] = vfrac_total - vsum;
    // the 1 here is the scaled total internal energy density
    residual[1] = 1.0 - esum;
    for (int m = 0; m < nmat - 1; ++m) {
      residual[2 + m] = press[m + 1] - press[m];
    }
  }

  PORTABLE_INLINE_FUNCTION
  bool CheckPTE() const {
    using namespace mix_params;
    Real mean_p = vfrac[0] * press[0];
    Real error_p = 0.0;
    for (int m = 1; m < nmat; ++m) {
      mean_p += vfrac[m] * press[m];
      error_p += residual[m + 1] * residual[m + 1];
    }
    error_p = std::sqrt(error_p);
    Real error_u = std::abs(residual[1]);
    // Check for convergence
    bool converged_p = (error_p < pte_rel_tolerance_p * std::abs(mean_p) ||
                        error_p < pte_abs_tolerance_p);
    bool converged_u = (error_u < pte_rel_tolerance_e || error_u < pte_abs_tolerance_e);
    return converged_p && converged_u;
  }

  PORTABLE_INLINE_FUNCTION
  void Jacobian() const {
    using namespace mix_params;
    Real dedT_sum = 0.0;
    for (int m = 0; m < nmat; m++) {
      //////////////////////////////
      // perturb volume fractions
      //////////////////////////////
      Real dv = (vfrac[m] < 0.5 ? 1.0 : -1.0) * vfrac[m] * derivative_eps;
      const Real vf_pert = vfrac[m] + dv;
      const Real rho_pert = rhobar[m] / vf_pert;

      Real p_pert{};
      Real e_pert =
          eos[m].InternalEnergyFromDensityTemperature(rho_pert, Tnorm * Tequil, Cache[m]);
      p_pert = this->GetPressureFromPreferred(eos[m], rho_pert, Tnorm * Tequil, e_pert,
                                              Cache[m], false) /
               uscale;
      dpdv[m] = (p_pert - press[m]) / dv;
      dedv[m] = (rhobar[m] * e_pert / uscale - u[m]) / dv;
      //////////////////////////////
      // perturb temperature
      //////////////////////////////
      Real dT = Tequil * derivative_eps;
      e_pert = eos[m].InternalEnergyFromDensityTemperature(rho[m], Tnorm * (Tequil + dT),
                                                           Cache[m]);
      p_pert = this->GetPressureFromPreferred(eos[m], rho[m], Tnorm * (Tequil + dT),
                                              e_pert, Cache[m], false) /
               uscale;
      dpdT[m] = (p_pert - press[m]) / dT;
      dedT_sum += (rhobar[m] * e_pert / uscale - u[m]) / dT;
    }

    // Fill in the Jacobian
    for (int i = 0; i < neq * neq; ++i)
      jacobian[i] = 0.0;
    for (int m = 0; m < nmat; ++m) {
      jacobian[m] = 1.0;
      jacobian[neq + m] = dedv[m];
    }
    jacobian[neq + nmat] = dedT_sum;
    for (int m = 0; m < nmat - 1; m++) {
      const int ind = MatIndex(2 + m, m);
      jacobian[ind] = dpdv[m];
      jacobian[ind + 1] = -dpdv[m + 1];
      jacobian[MatIndex(2 + m, nmat)] = dpdT[m] - dpdT[m + 1];
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real ScaleDx() const {
    using namespace mix_params;
    Real scale = 1.0;
    // control how big of a step toward vfrac = 0 is allowed
    for (int m = 0; m < nmat; ++m) {
      if (scale * dx[m] < -vfrac_safety_fac * vfrac[m]) {
        scale = -vfrac_safety_fac * vfrac[m] / dx[m];
      }
    }
    const Real Tnew = Tequil + scale * dx[nmat];
    // control how big of a step toward rho = rho(Pmin) is allowed
    for (int m = 0; m < nmat; m++) {
      const Real rho_min =
          std::max(eos[m].RhoPmin(Tnorm * Tequil), eos[m].RhoPmin(Tnorm * Tnew));
      const Real alpha_max = rhobar[m] / rho_min;
      if (scale * dx[m] > 0.5 * (alpha_max - vfrac[m])) {
        scale = 0.5 * (alpha_max - vfrac[m]) / dx[m];
      }
    }
    // control how big of a step toward T = 0 is allowed
    if (scale * dx[nmat] < -0.95 * Tequil) {
      scale = -0.95 * Tequil / dx[nmat];
    }
    // Now apply the overall scaling
    for (int i = 0; i < neq; ++i)
      dx[i] *= scale;
    return scale;
  }

  // Update the solution and return new residual.  Possibly called repeatedly with
  // different scale factors as part of a line search
  PORTABLE_INLINE_FUNCTION
  Real TestUpdate(const Real scale) {
    if (scale == 1.0) {
      Ttemp = Tequil;
      for (int m = 0; m < nmat; ++m)
        vtemp[m] = vfrac[m];
    }
    Tequil = Ttemp + scale * dx[nmat];
    for (int m = 0; m < nmat; ++m) {
      vfrac[m] = vtemp[m] + scale * dx[m];
      rho[m] = rhobar[m] / vfrac[m];
      u[m] = rhobar[m] * eos[m].InternalEnergyFromDensityTemperature(
                             rho[m], Tnorm * Tequil, Cache[m]);
      sie[m] = u[m] / rhobar[m];
      u[m] /= uscale;
      temp[m] = Tequil;
      press[m] = this->GetPressureFromPreferred(eos[m], rho[m], Tnorm * Tequil, sie[m],
                                                Cache[m], false) /
                 uscale;
    }
    Residual();
    return ResidualNorm();
  }

 private:
  Real *dpdv, *dedv, *dpdT, *vtemp;
  Real Tequil, Ttemp;
};

// fixed temperature solver
inline int PTESolverFixedTRequiredScratch(const int nmat) {
  int neq = nmat;
  return neq * neq                 // jacobian
         + 4 * neq                 // dx, residual, and sol_scratch
         + 2 * nmat                // rhobar and u in base
         + 2 * nmat                // nmat sized arrays in fixed T solver
         + MAX_NUM_LAMBDAS * nmat; // the cache
}
inline size_t PTESolverFixedTRequiredScratchInBytes(const int nmat) {
  return PTESolverFixedTRequiredScratch(nmat) * sizeof(Real);
}

template <typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
class PTESolverFixedT : public mix_impl::PTESolverBase<EOSIndexer, RealIndexer> {
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::AssignIncrement;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::nmat;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::neq;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::ResidualNorm;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::vfrac_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::eos;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::rho;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::vfrac;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::sie;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::temp;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::press;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::rho_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::uscale;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::MatIndex;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::jacobian;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::residual;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::dx;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::u;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::rhobar;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::Cache;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::Tnorm;

 public:
  // template the ctor to get type deduction/universal references prior to c++17
  // allow the type of the temperature array to be different, potentially a const Real*
  template <typename EOS_t, typename Real_t, typename CReal_t, typename Lambda_t>
  PORTABLE_INLINE_FUNCTION
  PTESolverFixedT(const int nmat, EOS_t &&eos, const Real vfrac_tot, const Real T_true,
                  Real_t &&rho, Real_t &&vfrac, Real_t &&sie, CReal_t &&temp,
                  Real_t &&press, Lambda_t &&lambda, Real *scratch)
      : mix_impl::PTESolverBase<EOSIndexer, RealIndexer>(nmat, nmat, eos, vfrac_tot, 1.0,
                                                         rho, vfrac, sie, temp, press,
                                                         scratch, T_true) {
    dpdv = AssignIncrement(scratch, nmat);
    vtemp = AssignIncrement(scratch, nmat);
    Tequil = T_true;
    Ttemp = T_true;
    Tnorm = 1.0;
    // TODO(JCD): use whatever lambdas are passed in
    /*for (int m = 0; m < nmat; m++) {
      if (lambda[m] != nullptr) Cache[m] = lambda[m];
    }*/
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
    for (int m = 0; m < nmat; m++) {
      // volume fractions have been potentially reset to ensure densitites are
      // larger than rho(Pmin(Tequil)); set the physical density to reflect
      // this change in volume fraction
      rho[m] = rhobar[m] / vfrac[m];
      sie[m] = eos[m].InternalEnergyFromDensityTemperature(rho[m], Tequil, Cache[m]);
      uscale += sie[m] * rho[m];
      // note the scaling of pressure
      press[m] = eos[m].PressureFromDensityTemperature(rho[m], Tequil, Cache[m]);
    }
    for (int m = 0; m < nmat; ++m) {
      press[m] /= uscale;
      u[m] = sie[m] * rhobar[m] / uscale;
    }
    Residual();
    return ResidualNorm();
  }

  PORTABLE_INLINE_FUNCTION
  void Residual() const {
    Real vsum = 0.0;
    for (int m = 0; m < nmat; ++m) {
      vsum += vfrac[m];
    }
    // the 1 here is the scaled volume fraction
    residual[0] = vfrac_total - vsum;
    for (int m = 0; m < nmat - 1; ++m) {
      residual[1 + m] = press[m] - press[m + 1];
    }
  }

  PORTABLE_INLINE_FUNCTION
  bool CheckPTE() const {
    using namespace mix_params;
    Real mean_p = vfrac[0] * press[0];
    Real error_p = 0;
    for (int m = 1; m < nmat; ++m) {
      mean_p += vfrac[m] * press[m];
      error_p += residual[m + 1] * residual[m + 1];
    }
    error_p = std::sqrt(error_p);
    Real error_v = std::abs(residual[0]);
    // Check for convergence
    bool converged_p = (error_p < pte_rel_tolerance_p * std::abs(mean_p) ||
                        error_p < pte_abs_tolerance_p);
    bool converged_v = (error_v < pte_rel_tolerance_e || error_v < pte_abs_tolerance_e);
    return converged_p && converged_v;
  }

  PORTABLE_INLINE_FUNCTION
  void Jacobian() const {
    using namespace mix_params;
    for (int m = 0; m < nmat; m++) {
      //////////////////////////////
      // perturb volume fractions
      //////////////////////////////
      Real dv = (vfrac[m] < 0.5 ? 1.0 : -1.0) * vfrac[m] * derivative_eps;
      const Real vf_pert = vfrac[m] + dv;
      const Real rho_pert = rhobar[m] / vf_pert;

      Real p_pert =
          eos[m].PressureFromDensityTemperature(rho_pert, Tequil, Cache[m]) / uscale;
      dpdv[m] = (p_pert - press[m]) / dv;
    }

    // Fill in the Jacobian
    for (int i = 0; i < neq * neq; ++i)
      jacobian[i] = 0.0;
    for (int m = 0; m < nmat; ++m) {
      jacobian[m] = 1.0;
    }
    for (int m = 0; m < nmat - 1; m++) {
      jacobian[MatIndex(m + 1, m)] = -dpdv[m];
      jacobian[MatIndex(m + 1, m + 1)] = dpdv[m + 1];
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real ScaleDx() const {
    using namespace mix_params;
    Real scale = 1.0;
    // control how big of a step toward vfrac = 0 is allowed
    for (int m = 0; m < nmat; ++m) {
      if (scale * dx[m] < -vfrac_safety_fac * vfrac[m]) {
        scale = -vfrac_safety_fac * vfrac[m] / dx[m];
      }
    }
    // control how big of a step toward rho = rho(Pmin) is allowed
    for (int m = 0; m < nmat; m++) {
      const Real rho_min = eos[m].RhoPmin(Tequil);
      const Real alpha_max = rhobar[m] / rho_min;
      if (scale * dx[m] > 0.5 * (alpha_max - vfrac[m])) {
        scale = 0.5 * (alpha_max - vfrac[m]) / dx[m];
      }
    }
    // Now apply the overall scaling
    for (int i = 0; i < neq; ++i)
      dx[i] *= scale;
    return scale;
  }

  // Update the solution and return new residual.  Possibly called repeatedly with
  // different scale factors as part of a line search
  PORTABLE_INLINE_FUNCTION
  Real TestUpdate(const Real scale) {
    if (scale == 1.0) {
      for (int m = 0; m < nmat; ++m)
        vtemp[m] = vfrac[m];
    }
    for (int m = 0; m < nmat; ++m) {
      vfrac[m] = vtemp[m] + scale * dx[m];
      rho[m] = rhobar[m] / vfrac[m];
      u[m] = rhobar[m] *
             eos[m].InternalEnergyFromDensityTemperature(rho[m], Tequil, Cache[m]);
      sie[m] = u[m] / rhobar[m];
      u[m] /= uscale;
      press[m] = eos[m].PressureFromDensityTemperature(rho[m], Tequil, Cache[m]) / uscale;
    }
    Residual();
    return ResidualNorm();
  }

 private:
  Real *dpdv, *vtemp;
  Real Tequil, Ttemp;
};

// fixed P solver
inline int PTESolverFixedPRequiredScratch(const int nmat) {
  int neq = nmat + 1;
  return neq * neq                 // jacobian
         + 4 * neq                 // dx, residual, and sol_scratch
         + 2 * nmat                // all the nmat sized arrays in base
         + 3 * nmat                // all the nmat sized arrays in fixedP
         + MAX_NUM_LAMBDAS * nmat; // the cache
}
inline size_t PTESolverFixedPRequiredScratchInBytes(const int nmat) {
  return PTESolverFixedPRequiredScratch(nmat) * sizeof(Real);
}

template <typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
class PTESolverFixedP : public mix_impl::PTESolverBase<EOSIndexer, RealIndexer> {
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::InitBase;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::AssignIncrement;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::nmat;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::neq;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::ResidualNorm;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::vfrac_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::sie_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::eos;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::rho;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::vfrac;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::sie;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::temp;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::press;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::rho_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::uscale;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::MatIndex;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::TryIdealPTE;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::jacobian;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::residual;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::dx;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::u;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::rhobar;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::Cache;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::Tnorm;

 public:
  // template the ctor to get type deduction/universal references prior to c++17
  template <typename EOS_t, typename Real_t, typename CReal_t, typename Lambda_t>
  PORTABLE_INLINE_FUNCTION
  PTESolverFixedP(const int nmat, EOS_t &&eos, const Real vfrac_tot, const Real P,
                  Real_t &&rho, Real_t &&vfrac, Real_t &&sie, Real_t &&temp,
                  CReal_t &&press, Lambda_t &&lambda, Real *scratch)
      : mix_impl::PTESolverBase<EOSIndexer, RealIndexer>(nmat, nmat + 1, eos, vfrac_tot,
                                                         1.0, rho, vfrac, sie, temp,
                                                         press, scratch, 0.0) {
    dpdv = AssignIncrement(scratch, nmat);
    dpdT = AssignIncrement(scratch, nmat);
    vtemp = AssignIncrement(scratch, nmat);
    Pequil = P;
    // TODO(JCD): use whatever lambdas are passed in
    /*for (int m = 0; m < nmat; m++) {
      if (lambda[m] != nullptr) Cache[m] = lambda[m];
    }*/
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
    for (int m = 0; m < nmat; m++) {
      // scaled initial guess for temperature is just 1
      temp[m] = 1.0;
      sie[m] = eos[m].InternalEnergyFromDensityTemperature(rho[m], Tguess, Cache[m]);
      uscale += sie[m] * rho[m];
    }

    // note the scaling of the material internal energy densities
    for (int m = 0; m < nmat; ++m) {
      u[m] = sie[m] * rhobar[m] / uscale;
      press[m] = eos[m].PressureFromDensityTemperature(rho[m], Tguess, Cache[m]) / uscale;
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
    for (int m = 0; m < nmat; ++m) {
      vsum += vfrac[m];
      residual[m] = Pequil / uscale - press[m];
    }
    residual[nmat] = vfrac_total - vsum;
  }

  PORTABLE_INLINE_FUNCTION
  bool CheckPTE() const {
    using namespace mix_params;
    Real error_p = 0;
    for (int m = 0; m < nmat; ++m) {
      error_p += residual[m] * residual[m];
    }
    error_p = std::sqrt(error_p);
    error_p *= uscale;
    Real error_v = std::abs(residual[neq - 1]);
    // Check for convergence
    bool converged_p = (error_p < pte_rel_tolerance_p || error_p < pte_abs_tolerance_p);
    bool converged_v = (error_v < pte_rel_tolerance_e || error_v < pte_abs_tolerance_e);
    return converged_p && converged_v;
  }

  PORTABLE_INLINE_FUNCTION
  void Jacobian() const {
    using namespace mix_params;
    for (int m = 0; m < nmat; m++) {
      //////////////////////////////
      // perturb volume fractions
      //////////////////////////////
      Real dv = (vfrac[m] < 0.5 ? 1.0 : -1.0) * vfrac[m] * derivative_eps;
      const Real vf_pert = vfrac[m] + dv;
      const Real rho_pert = rhobar[m] / vf_pert;

      Real p_pert{};
      Real e_pert{};
      p_pert = this->GetPressureFromPreferred(eos[m], rho_pert, Tnorm * Tequil, e_pert,
                                              Cache[m], true) /
               uscale;
      dpdv[m] = (p_pert - press[m]) / dv;
      //////////////////////////////
      // perturb temperature
      //////////////////////////////
      Real dT = Tequil * derivative_eps;

      p_pert = this->GetPressureFromPreferred(eos[m], rho[m], Tnorm * (Tequil + dT),
                                              e_pert, Cache[m], true) /
               uscale;
      dpdT[m] = (p_pert - press[m]) / dT;
    }

    // Fill in the Jacobian
    for (int i = 0; i < neq * neq; ++i)
      jacobian[i] = 0.0;
    for (int m = 0; m < nmat; m++) {
      jacobian[MatIndex(m, m)] = dpdv[m];
      jacobian[MatIndex(m, nmat)] = dpdT[m];
    }
    for (int m = 0; m < neq - 1; ++m) {
      jacobian[MatIndex(neq - 1, m)] = 1.0;
    }
  }

  PORTABLE_INLINE_FUNCTION
  Real ScaleDx() const {
    using namespace mix_params;
    Real scale = 1.0;
    // control how big of a step toward vfrac = 0 is allowed
    for (int m = 0; m < nmat; ++m) {
      if (scale * dx[m] < -vfrac_safety_fac * vfrac[m]) {
        scale = -vfrac_safety_fac * vfrac[m] / dx[m];
      }
    }
    const Real Tnew = Tequil + scale * dx[nmat];
    // control how big of a step toward rho = rho(Pmin) is allowed
    for (int m = 0; m < nmat; m++) {
      const Real rho_min =
          std::max(eos[m].RhoPmin(Tnorm * Tequil), eos[m].RhoPmin(Tnorm * Tnew));
      const Real alpha_max = rhobar[m] / rho_min;
      if (scale * dx[m] > 0.5 * (alpha_max - vfrac[m])) {
        scale = 0.5 * (alpha_max - vfrac[m]) / dx[m];
      }
    }
    // control how big of a step toward T = 0 is allowed
    if (scale * dx[nmat] < -0.95 * Tequil) {
      scale = -0.95 * Tequil / dx[nmat];
    }
    // Now apply the overall scaling
    for (int i = 0; i < neq; ++i)
      dx[i] *= scale;
    return scale;
  }

  // Update the solution and return new residual.  Possibly called repeatedly with
  // different scale factors as part of a line search
  PORTABLE_INLINE_FUNCTION
  Real TestUpdate(const Real scale) {
    if (scale == 1.0) {
      Ttemp = Tequil;
      for (int m = 0; m < nmat; ++m)
        vtemp[m] = vfrac[m];
    }
    Tequil = Ttemp + scale * dx[nmat];
    for (int m = 0; m < nmat; ++m) {
      vfrac[m] = vtemp[m] + scale * dx[m];
      rho[m] = rhobar[m] / vfrac[m];
      u[m] = rhobar[m] * eos[m].InternalEnergyFromDensityTemperature(
                             rho[m], Tnorm * Tequil, Cache[m]);
      press[m] = eos[m].PressureFromDensityTemperature(rho[m], Tnorm * Tequil, Cache[m]) /
                 uscale;
      sie[m] = u[m] / rhobar[m];
      u[m] /= uscale;
      temp[m] = Tequil;
    }
    Residual();
    return ResidualNorm();
  }

 private:
  Real *dpdv, *dpdT, *vtemp;
  Real Tequil, Ttemp, Pequil;
};

inline int PTESolverRhoURequiredScratch(const int nmat) {
  int neq = 2 * nmat;
  return neq * neq                 // jacobian
         + 4 * neq                 // dx, residual, and sol_scratch
         + 8 * nmat                // all the nmat sized arrays
         + MAX_NUM_LAMBDAS * nmat; // the cache
}
inline size_t PTESolverRhoURequiredScratchInBytes(const int nmat) {
  return PTESolverRhoURequiredScratch(nmat) * sizeof(Real);
}

template <typename EOSIndexer, typename RealIndexer, typename LambdaIndexer>
class PTESolverRhoU : public mix_impl::PTESolverBase<EOSIndexer, RealIndexer> {
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::InitBase;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::AssignIncrement;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::nmat;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::neq;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::ResidualNorm;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::vfrac_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::sie_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::eos;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::rho;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::vfrac;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::sie;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::temp;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::press;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::rho_total;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::uscale;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::TryIdealPTE;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::MatIndex;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::jacobian;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::residual;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::dx;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::u;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::rhobar;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::Cache;
  using mix_impl::PTESolverBase<EOSIndexer, RealIndexer>::Tnorm;

 public:
  // template the ctor to get type deduction/universal references prior to c++17
  template <typename EOS_t, typename Real_t, typename Lambda_t>
  PORTABLE_INLINE_FUNCTION PTESolverRhoU(const int nmat, const EOS_t &&eos,
                                         const Real vfrac_tot, const Real sie_tot,
                                         Real_t &&rho, Real_t &&vfrac, Real_t &&sie,
                                         Real_t &&temp, Real_t &&press, Lambda_t &&lambda,
                                         Real *scratch, const Real Tguess = 0.0)
      : mix_impl::PTESolverBase<EOSIndexer, RealIndexer>(nmat, 2 * nmat, eos, vfrac_tot,
                                                         sie_tot, rho, vfrac, sie, temp,
                                                         press, scratch, Tguess) {
    dpdv = AssignIncrement(scratch, nmat);
    dtdv = AssignIncrement(scratch, nmat);
    dpde = AssignIncrement(scratch, nmat);
    dtde = AssignIncrement(scratch, nmat);
    vtemp = AssignIncrement(scratch, nmat);
    utemp = AssignIncrement(scratch, nmat);
    // TODO(JCD): use whatever lambdas are passed in
    /*for (int m = 0; m < nmat; m++) {
      if (lambda[m] != nullptr) Cache[m] = lambda[m];
    }*/
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
    for (int m = 0; m < nmat; ++m) {
      vsum += vfrac[m];
      esum += u[m];
    }
    residual[0] = vfrac_total - vsum;
    // the 1 here is the scaled total internal energy density
    residual[1] = 1.0 - esum;
    for (int m = 0; m < nmat - 1; ++m) {
      residual[2 + m] = press[m + 1] - press[m];
    }
    for (int m = nmat + 1; m < neq; m++) {
      residual[m] = temp[m - nmat] - temp[m - nmat - 1];
    }
  }

  PORTABLE_INLINE_FUNCTION
  bool CheckPTE() const {
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
    bool converged_p = (error_p < pte_rel_tolerance_p * std::abs(mean_p) ||
                        error_p < pte_abs_tolerance_p);
    bool converged_t =
        (error_t < pte_rel_tolerance_t * mean_t || error_t < pte_abs_tolerance_t);
    return (converged_p && converged_t);
  }

  PORTABLE_INLINE_FUNCTION
  void Jacobian() const {
    using namespace mix_params;
    for (int m = 0; m < nmat; m++) {
      //////////////////////////////
      // perturb volume fractions
      //////////////////////////////
      Real dv = (vfrac[m] < 0.5 ? 1.0 : -1.0) * vfrac[m] * derivative_eps;
      const Real vf_pert = vfrac[m] + dv;
      const Real rho_pert = rhobar[m] / vf_pert;

      Real p_pert;
      Real t_pert =
          eos[m].TemperatureFromDensityInternalEnergy(rho_pert, sie[m], Cache[m]) / Tnorm;
      p_pert = this->GetPressureFromPreferred(eos[m], rho_pert, Tnorm * t_pert, sie[m],
                                              Cache[m], false) /
               uscale;
      dpdv[m] = (p_pert - press[m]) / dv;
      dtdv[m] = (t_pert - temp[m]) / dv;
      //////////////////////////////
      // perturb energies
      //////////////////////////////
      const Real de = std::abs(u[m]) * derivative_eps;
      Real e_pert = (u[m] + de) / rhobar[m];

      t_pert =
          eos[m].TemperatureFromDensityInternalEnergy(rho[m], uscale * e_pert, Cache[m]) /
          Tnorm;
      p_pert = this->GetPressureFromPreferred(eos[m], rho[m], Tnorm * t_pert,
                                              uscale * e_pert, Cache[m], false) /
               uscale;
      dpde[m] = (p_pert - press[m]) / de;
      dtde[m] = (t_pert - temp[m]) / de;
      if (std::abs(dtde[m]) < 1.e-16) { // must be on the cold curve
        dtde[m] = derivative_eps;
      }
    }
    for (int i = 0; i < neq * neq; ++i)
      jacobian[i] = 0.0;
    // TODO(JCD): clean all this up with MatIndex
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
  }

  PORTABLE_INLINE_FUNCTION
  Real ScaleDx() const {
    using namespace mix_params;
    Real scale = 1.0;
    for (int m = 0; m < nmat; ++m) {
      // control how big of a step toward vfrac = 0 is allowed
      if (scale * dx[m] < -0.1 * vfrac[m]) {
        scale = -0.1 * vfrac[m] / dx[m];
      }
      // try to control steps toward T = 0
      const Real dt = (dtdv[m] * dx[m] + dtde[m] * dx[m + nmat]);
      if (scale * dt < -0.1 * temp[m]) {
        scale = -0.1 * temp[m] / dt;
      }
      const Real tt = temp[m] + scale * dt;
      const Real rho_min =
          std::max(eos[m].RhoPmin(Tnorm * temp[m]), eos[m].RhoPmin(Tnorm * tt));
      const Real alpha_max = rhobar[m] / rho_min;
      // control how big of a step toward rho = rho(Pmin) is allowed
      if (scale * dx[m] > 0.5 * (alpha_max - vfrac[m])) {
        scale = 0.5 * (alpha_max - vfrac[m]) / dx[m];
      }
    }
    // Now apply the overall scaling
    for (int i = 0; i < neq; ++i)
      dx[i] *= scale;
    return scale;
  }

  // Update the solution and return new residual.  Possibly called repeatedly with
  // different scale factors as part of a line search
  PORTABLE_INLINE_FUNCTION
  Real TestUpdate(const Real scale) const {
    if (scale == 1.0) {
      for (int m = 0; m < nmat; ++m) {
        vtemp[m] = vfrac[m];
        utemp[m] = u[m];
      }
    }
    for (int m = 0; m < nmat; ++m) {
      vfrac[m] = vtemp[m] + scale * dx[m];
      rho[m] = rhobar[m] / vfrac[m];
      u[m] = utemp[m] + scale * dx[nmat + m];
      sie[m] = uscale * u[m] / rhobar[m];
      temp[m] =
          eos[m].TemperatureFromDensityInternalEnergy(rho[m], sie[m], Cache[m]) / Tnorm;
      press[m] = this->GetPressureFromPreferred(eos[m], rho[m], Tnorm * temp[m], sie[m],
                                                Cache[m], false) /
                 uscale;
    }
    Residual();
    return ResidualNorm();
  }

 private:
  Real *dpdv, *dtdv, *dpde, *dtde, *vtemp, *utemp;
};

template <class System>
PORTABLE_INLINE_FUNCTION bool PTESolver(System &s) {
  using namespace mix_params;
  // initialize the system, fill in residual, and get its norm
  Real err = s.Init();

  bool converged = false;
  const int pte_max_iter = s.Nmat() * pte_max_iter_per_mat;
  const Real residual_tol = s.Nmat() * pte_residual_tolerance;
  auto &niter = s.Niter();
  for (niter = 0; niter < pte_max_iter; ++niter) {
    // Check for convergence
    converged = s.CheckPTE();
    if (converged) break;

    // compute the Jacobian
    s.Jacobian();

    // solve for the Newton step
    bool success = s.Solve();
    if (!success) {
      // do something to crash out?  Tell folks what happened?
      // printf("crashing out at iteration: %i\n", niter);
      converged = false;
      break;
    }

    // possibly scale the update to stay within reasonable bounds
    Real scale = s.ScaleDx();
    // const Real scale_save = scale;

    // Line search
    Real gradfdx = -2.0 * scale * err;
    scale = 1.0;
    Real err_old = err;
    err = s.TestUpdate(scale);
    if (err > err_old + line_search_alpha * gradfdx) {
      // backtrack
      Real err_mid = s.TestUpdate(0.5);
      if (err_mid < err && err_mid < err_old) {
        scale = 0.75 + 0.5 * (err_mid - err) / (err - 2.0 * err_mid + err_old);
      } else {
        scale = line_search_fac;
      }

      for (int line_iter = 0; line_iter < line_search_max_iter; line_iter++) {
        err = s.TestUpdate(scale);
        if (err < err_old + line_search_alpha * scale * gradfdx) break;
        scale *= line_search_fac;
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
  return converged;
}

} // namespace singularity

#endif // _SINGULARITY_EOS_CLOSURE_MIXED_CELL_MODELS_
