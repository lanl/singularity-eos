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

#ifndef _SINGULARITY_EOS_UTILS_ROOT_FINDING_HPP_
#define _SINGULARITY_EOS_UTILS_ROOT_FINDING_HPP_

// Implementation based on gsl root finder API
// Code originally taken from nubhlight, LA-UR-19-20336
// Miller, Ryan, Dolence, (2019). ApJ Supplements, 241(2), 30.
// TODO: this code still contains many "C"-isms.
// TOOD: if-logic could be removed for speed. Careful, though!
// ~JMM

// #include <iostream> // debug
// #include <iomanip>
#include <assert.h>
#include <math.h>
#include <ports-of-call/portability.hpp>
#include <stdio.h>
#include <tuple>

#define SINGULARITY_ROOT_DEBUG (0)
#define SINGULARITY_ROOT_VERBOSE (0)

#define SINGULARITY_MY_SIGN(x) (x > 0) - (x < 0)

namespace RootFinding1D {
constexpr const int SECANT_NITER_MAX{1000};
constexpr const int BISECT_NITER_MAX{1000};
constexpr const int BISECT_REG_MAX{1000};
constexpr const int NEWTON_RAPHSON_NITER_MAX{100};
enum class Status { SUCCESS = 0, FAIL = 1 };

/*
  // TODO: Something like this would be nice
  friend std::ostream& operator<< (std::ostream& os, const RootCounts& c) {
  Real tot = c.total();
  std::vector<Real> percents(c.nbins_,0);
  for (int i = 0; i < c.nbins_; i++) {
  percents[i] = (100.*c.counts_[i] / tot);
  }
  os << "\n********** ROOT FINDING *********\n"
  << "   ITERATIONS          PERCENTAGE\n";
  for (int i = 0; i < c.more_; i++) {
  os << "         "
  << std::right << std::setw(4)
  << i
  << "             "
  << std::right << std::setw(3)
  << percents[i] << "\n";
  }
  os << "         more             "
  << std::right << std::setw(3)
  << percents[c.more_] << "\n"
  << "*********************************\n"
  << std::endl;
  return os;
  }
*/
class RootCounts {
 private:
  static constexpr std::size_t nbins_{15};
  mutable Real counts_[nbins_];

 public:
  PORTABLE_INLINE_FUNCTION
  RootCounts() {
    for (std::size_t i{0}; i < nbins_; ++i)
      counts_[i] = 0;
  }
  PORTABLE_INLINE_FUNCTION void reset() {
    for (std::size_t i{0}; i < nbins_; ++i)
      counts_[i] = 0;
  }
  PORTABLE_INLINE_FUNCTION void increment(std::size_t i) const {
    assert(i < nbins_ && i >= 0);
#ifdef PORTABILITY_STRATEGY_NONE
    counts_[i] += 1;
#endif // PORTABILITY_STRATEGY_NONE
  }
  PORTABLE_INLINE_FUNCTION Real total() const {
    Real tot{1.e-20};
    for (std::size_t i{0}; i < nbins_; ++i)
      tot += counts_[i];
    return tot;
  }
  PORTABLE_INLINE_FUNCTION const Real &operator[](const std::size_t i) const {
    assert(i < nbins_ && i >= 0);
    return counts_[i];
  }
  PORTABLE_INLINE_FUNCTION Real &operator[](const std::size_t i) {
    assert(i < nbins_ && i >= 0);
    return counts_[i];
  }
  PORTABLE_INLINE_FUNCTION void print_counts() const {
    for (std::size_t i{0}; i < nbins_; ++i)
      printf("%e\n", counts_[i]);
  }
  PORTABLE_INLINE_FUNCTION std::size_t nBins() const { return nbins_; }
  PORTABLE_INLINE_FUNCTION std::size_t more() const { return nbins_ - 1; }
};

PORTABLE_INLINE_FUNCTION bool check_bracket(const Real ya, const Real yb) {
  return (ya * yb <= 0.0);
}

template <typename T>
PORTABLE_INLINE_FUNCTION bool set_bracket(const T &f, Real &a, const Real guess, Real &b,
                                          Real &ya, const Real yg, Real &yb,
                                          const bool &verbose = false) {
  constexpr std::size_t max_search_depth = 6;
  Real dx = b - a;
  for (std::size_t level = 0; level < max_search_depth; level++) {
    const std::size_t nlev = (1 << level);
    for (std::size_t i = 0; i < nlev; i++) {
      const Real x = a + (i + 0.5) * dx;
      const Real yx = f(x);
      if (check_bracket(yx, yg)) {
        if (x < guess) {
          a = x;
          ya = yx;
          b = guess;
          yb = yg;
        } else {
          a = guess;
          ya = yg;
          b = x;
          yb = yx;
        }
        return true;
      }
    }
  }
  // if we get here then we failed to bracket a root
  if (verbose) {
    printf("set_bracket failed to bound a root! %.14e %.14e %.14e %.14e %.14e %.14e\n", a,
           guess, b, ya, yg, yb);
  }
  return false;
}

// solves for f(x,params) - ytarget = 0
template <typename T>
PORTABLE_INLINE_FUNCTION Status regula_falsi(const T &f, const Real ytarget,
                                             const Real guess, Real a, Real b,
                                             const Real xtol, const Real ytol,
                                             Real &xroot,
                                             const RootCounts *counts = nullptr,
                                             const bool &verbose = false) {
  constexpr std::size_t max_iter = SECANT_NITER_MAX;
  auto func = [&](const Real x) { return f(x) - ytarget; };
  Real ya = func(a);
  Real yg = func(guess);
  Real yb;

  if (check_bracket(ya, yg)) {
    b = guess;
    yb = yg;
  } else {
    yb = func(b);
    if (check_bracket(yg, yb)) {
      a = guess;
      ya = yg;
    } else {
      // ya, yg, and yb have the same sign
      if (!set_bracket(func, a, guess, b, ya, yg, yb, verbose)) {
        if (verbose) {
          printf("regula_falsi failed! %.14e %.14e %.14e %.14e\n", ytarget, guess, a, b);
        }
        return Status::FAIL;
      }
    }
  }

  Real sign = (ya < 0 ? 1.0 : -1.0);
  ya *= sign;
  yb *= sign;

  std::size_t b1 = 0;
  std::size_t b2 = 0;
  std::size_t iteration_count = 0;
  while (b - a > 2.0 * xtol && (std::abs(ya) > ytol || std::abs(yb) > ytol) &&
         iteration_count < max_iter) {
    Real c = (a * yb - b * ya) / (yb - ya);
    // guard against roundoff because ya or yb is sufficiently close to zero
    if (c == a) {
      b = a;
      continue;
    } else if (c == b) {
      a = b;
      continue;
    }
    Real yc = sign * func(c);
    if (yc > 0.0) {
      b = c;
      yb = yc;
      b1++;
      ya *= (b1 > 1 ? 0.5 : 1.0);
      b2 = 0;
    } else if (yc < 0.0) {
      a = c;
      ya = yc;
      b2++;
      yb *= (b2 > 1 ? 0.5 : 1.0);
      b1 = 0;
    } else {
      a = c;
      b = c;
    }
    iteration_count++;
  }
  auto status = Status::SUCCESS;
  if (iteration_count == max_iter) {
    if (verbose) {
      printf("root finding reached the maximum number of iterations.  likely not "
             "converged\n");
    }
    status = Status::FAIL;
  }
  if (counts != nullptr) {
    if (iteration_count < counts->nBins()) {
      counts->increment(iteration_count);
    } else {
      counts->increment(counts->more());
    }
  }
  xroot = 0.5 * (a + b);
  return status;
}

// solves for f(x,params) - ytarget = 0
// WARNING: this root finding expects a different callable f than the other
// root finding methods. f should return a tuple of (f(x), f'(x)) where f'(x)
// is the derivative of f with respect to x.
template <typename T>
PORTABLE_INLINE_FUNCTION Status newton_raphson(const T &f, const Real ytarget,
                                               const Real guess, const Real a,
                                               const Real b, const Real ytol, Real &xroot,
                                               const RootCounts *counts = nullptr,
                                               const bool &verbose = false,
                                               const bool &fail_on_bound_root = true) {

  constexpr std::size_t max_iter = NEWTON_RAPHSON_NITER_MAX;
  Real _x = guess;
  Real _xold = 0.0;
  auto status = Status::SUCCESS;

  Real yg;
  Real dfunc;

  std::size_t iter;

  for (iter = 0; iter < max_iter; iter++) {
    std::tie(yg, dfunc) = f(_x); // C++11 tuple unpacking

    // check if we are converged already
    if (std::abs(yg - ytarget) < std::abs(ytol * ytarget)) break;

    // not converged; compute the next step
    _xold = _x;
    _x = _x - (yg - ytarget) / dfunc;

    // check if we are out of bounds
    // CAUTION: we do not set the root to the boundary value in this case
    // because one might want to handle this on a case by case basis
    // (e.g. if the boundary is a physical boundary, then one might want to
    // set the root to the boundary value).
    // Per default, we fail if the root is out of bounds controlled by
    // fail_on_bound_root.
    if ((_x <= a && _xold <= a) || (_x >= b && _xold >= b)) {
      if (verbose) {
        printf("newton_raphson out of bounds! %.14e %.14e %.14e %.14e\n", ytarget, guess,
               a, b);
      }
      if (fail_on_bound_root) {
        status = Status::FAIL;
      }
      break;
    }
    _x = std::max(std::min(_x, b), a);
  }
  if (iter >= max_iter) {
    if (verbose) {
      printf("root finding reached the maximum number of iterations.  likely not "
             "converged\n");
    }
    status = Status::FAIL;
  }

  if (counts != nullptr) {
    if (iter < counts->nBins()) {
      counts->increment(iter);
    } else {
      counts->increment(counts->more());
    }
  }
  xroot = _x;
  return status;
}

template <typename T>
PORTABLE_INLINE_FUNCTION Status findRoot(const T &f, const Real ytarget, Real xguess,
                                         const Real xmin, const Real xmax,
                                         const Real xtol, const Real ytol, Real &xroot,
                                         const RootCounts *counts = nullptr) {
  Status status;

  // first check if we're at the max or min values
  const Real fmax = f(xmax);
  const Real errmax = fabs(fmax - ytarget) / (fabs(fmax) + ytol);
  if (errmax < ytol) {
    xroot = xmax;
    if (counts != nullptr) {
      counts->increment(0);
    }
    return Status::SUCCESS;
  }
  const Real fmin = f(xmin);
  const Real errmin = fabs(fmin - ytarget) / (fabs(fmin) + ytol);
  if (errmin < ytol) {
    xroot = xmin;
    if (counts != nullptr) {
      counts->increment(0);
    }
    return Status::SUCCESS;
  }

  if (xguess >= xmax) xguess = xmax - xtol;
  if (xguess <= xmin) xguess = xmin + xtol;

  // Next try Secant
  status = secant(f, ytarget, xguess, xmin, xmax, xtol, ytol, xroot, counts);
  if (status == Status::SUCCESS) return status;

#if SINGULARITY_ROOT_DEBUG
  if (isnan(xroot)) {
    fprintf(stderr, "xroot is nan after secant\n");
  }
#endif

// Secant failed. Try bisection.
#if SINGULARITY_ROOT_VERBOSE
  fprintf(stderr,
          "\n\nRoot finding. Secant failed. Trying bisection.\n"
          "\txguess  = %.10g\n"
          "\tytarget = %.10g\n"
          "\txmin    = %.10g\n"
          "\txmax    = %.10g\n",
          xguess, ytarget, xmin, xmax);
#endif
  status = bisect(f, ytarget, xguess, xmin, xmax, xtol, ytol, xroot, counts);

  // Check for something horrible happening
  if (isnan(xroot) || isinf(xroot)) {
#if SINGULARITY_ROOT_DEBUG
    fprintf(stderr, "xroot is nan after bisection\n");
#endif
    return Status::FAIL;
  }
  if (xroot < xmin) return Status::FAIL;
  if (xroot > xmax) return Status::FAIL;

  return status;
}

template <typename T>
PORTABLE_INLINE_FUNCTION Status secant(const T &f, const Real ytarget, const Real xguess,
                                       const Real xmin, const Real xmax, const Real xtol,
                                       const Real ytol, Real &xroot,
                                       const RootCounts *counts = nullptr) {
  Real dx;
  Real x_last, y, yp, ym, dyNum, dyDen, dy;

  Real x = xguess;
  std::size_t iter{0};
  for (iter = 0; iter < SECANT_NITER_MAX; ++iter) {
    x_last = x;
    dx = fabs(1.e-7 * x) + xtol;
    y = f(x) - ytarget;
    yp = f(x + dx);
    ym = f(x - dx);
    dyNum = yp - ym;
    dyDen = (2. * dx);
    dy = dyNum / dyDen;
    x -= y / dy;
    if (x < xmin) x = xmin;
    if (x > xmax) x = xmax;
    if (isnan(x) || isinf(x)) {
      // can't recover from this
#if SINGULARITY_ROOT_DEBUG
      fprintf(stderr,
              "\n\n[secant]: NAN or out-of-bounds detected!\n"
              "\txguess  = %.10e\n"
              "\tytarget = %.10e\n"
              "\tx       = %.10e\n"
              "\tx_last  = %.10e\n"
              "\tx_min   = %.10e\n"
              "\tx_max   = %.10e\n"
              "\ty       = %.10e\n"
              "\tdx      = %.10e\n"
              "\typ      = %.10e\n"
              "\tym      = %.10e\n"
              "\tdyNum   = %.10e\n"
              "\tdyDen   = %.10e\n"
              "\tdy      = %.10e\n"
              "\titer    = %d\n"
              "\tsign x  = %d\n",
              xguess, ytarget, x, x_last, xmin, xmax, y, dx, yp, ym, dyNum, dyDen, dy,
              iter, (int)SINGULARITY_MY_SIGN(x));
#endif
      if (counts != nullptr) {
        counts->increment(counts->more());
      }
      return Status::FAIL;
    }
    if (fabs(x - x_last) / (fabs(x) + xtol) < xtol) break;
  }
  // increment iter
  ++iter;
  if (counts != nullptr) {
    if (iter < counts->nBins()) {
      counts->increment(iter);
    } else {
      counts->increment(counts->more());
    }
  }
  xroot = x;

  y = f(x);
  const Real frac_error = fabs(y - ytarget) / (fabs(y) + ytol);
#if SINGULARITY_ROOT_DEBUG
  if (frac_error > ytol) {
    fprintf(stderr,
            "\n\n[secant]: Failed via too large yerror.\n"
            "\tfractional error = %.10e\n"
            "\tx                = %.10e\n"
            "\ty                = %.10e\n"
            "\tytarget          = %.10e\n"
            "\typ               = %.10e\n"
            "\tym               = %.10e\n"
            "\tdy               = %.10e\n"
            "\tdx               = %.10e\n"
            "\titer             = %d\n",
            frac_error, x, y, ytarget, yp, ym, dy, dx, iter);
  }
  if (fabs(x - x_last) > xtol) {
    fprintf(stderr,
            "\n\n[secant]: failed via dx too big.\n"
            "\tfractional error = %.10e\n"
            "\tx                = %.10e\n"
            "\tx_last           = %.10e\n"
            "\tdx               = %.10e\n"
            "\ty                = %.10e\n"
            "\tytarget          = %.10e\n"
            "\typ               = %.10e\n"
            "\tym               = %.10e\n"
            "\tdy               = %.10e\n"
            "\tdx               = %.10e\n"
            "\titer             = %d\n",
            frac_error, x, x_last, fabs(x - x_last), y, ytarget, yp, ym, dy, dx, iter);
  }
#endif

  const int secant_failed =
      (fabs(x - x_last) > xtol || fabs(frac_error) > ytol || isnan(x) || isinf(x));
  return secant_failed ? Status::FAIL : Status::SUCCESS;
}

template <typename T>
PORTABLE_INLINE_FUNCTION Status bisect(const T &f, const Real ytarget, const Real xguess,
                                       const Real xmin, const Real xmax, const Real xtol,
                                       const Real ytol, Real &xroot,
                                       const RootCounts *counts = nullptr) {
  Real xl, xr, fl, fr, dx;

  Real grow = 0.01;
  Real x = xguess;
  if (fabs(x) < xtol) {
    x += 2. * xtol;
  }
  // do { // Try to find reasonable region for bisection
  for (std::size_t i{0}; i < BISECT_REG_MAX; ++i) {
    dx = fabs(grow * x);
    xl = x - dx;
    xr = x + dx;
    fl = f(xl) - ytarget;
    fr = f(xr) - ytarget;
    grow *= 1.1;
    if (fl * fr < 0.0 || xl < xmin || xr > xmax) break;
    if (i > BISECT_REG_MAX - 2) {
#if SINGULARITY_ROOT_DEBUG
      fprintf(stderr,
              "\n\n[Bisect]: expanding region failed\n"
              "\txl   = %.10g\n"
              "\txr   = %.10g\n"
              "\tfl   = %.10g\n"
              "\tfr   = %.10g\n"
              "\tdx   = %.10g\n"
              "\tgrow = %.10g\n",
              xl, xr, fl, fr, dx, grow);
#endif
      return Status::FAIL;
    }
  }

  // force back onto the bisection region
  if (xr > xmax) {
    xr = xmax;
    fr = f(xr) - ytarget;
  }
  if (xl < xmin) {
    xl = xmin;
    fl = f(xl) - ytarget;
  }

  // if they have the same sign, change that.
  // if we can't fix it, fail.
  if (fl * fr > 0) {
    xl = xmin;
    fl = f(xl) - ytarget;
    if (fl * fr > 0) {
      xr = xmax;
      fr = f(xr) - ytarget;
      if (fl * fr > 0) {
#if SINGULARITY_ROOT_DEBUG
        Real il = f(xl);
        Real ir = f(xr);
        fprintf(stderr,
                "\n\n[bisect]: fl*fr > 0!\n"
                "\txguess  = %.10e\n"
                "\tytarget = %.10e\n"
                "\txl      = %.10e\n"
                "\txr      = %.10e\n"
                "\tfl      = %.10e\n"
                "\tfr      = %.10e\n"
                "\til      = %.10e\n"
                "\tir      = %.10e\n",
                xguess, ytarget, xl, xr, fl, fr, il, ir);
        std::size_t nx = 300;
        Real dx = (xmax - xmin) / (nx - 1);
        fprintf(stderr, "Area map:\nx\ty\n");
        for (std::size_t i = 0; i < nx; i++) {
          fprintf(stderr, "%.4f\t%.4e\n", x + i * dx, f(x + i * dx));
        }
#endif
        return Status::FAIL;
      }
    }
  }

  for (std::size_t i{0}; i < BISECT_NITER_MAX; ++i) {
    Real xm = 0.5 * (xl + xr);
    Real fm = f(xm) - ytarget;
    if (fl * fm <= 0) {
      xr = xm;
      fr = fm;
    } else {
      xl = xm;
      fl = fm;
    }
    if (xr - xl < xtol) {
      xroot = 0.5 * (xl + xr);
      return Status::SUCCESS;
    }
  }

  xroot = 0.5 * (xl + xr);

  if (isnan(xroot)) {
#if SINGULARITY_ROOT_DEBUG
    Real il = f(xl);
    Real ir = f(xr);
    fprintf(stderr,
            "\n\n[bisect]: NAN DETECTED!\n"
            "\txguess  = %.10e\n"
            "\tytarget = %.10e\n"
            "\txl      = %.10e\n"
            "\txr      = %.10e\n"
            "\tdx      = %.10e\n"
            "\tgrow    = %.10e\n"
            "\txtol    = %.10e\n"
            "\tfl      = %.10e\n"
            "\tfr      = %.10e\n"
            "\til      = %.10e\n"
            "\tir      = %.10e\n"
            "\txmin    = %.10e\n"
            "\txmax    = %.10e\n",
            xguess, ytarget, xl, xr, dx, grow, xtol, fl, fr, il, ir, xmin, xmax);
#endif
  }
  return Status::FAIL;
}
// ----------------------------------------------------------------------
} // namespace RootFinding1D

#undef SINGULARITY_ROOT_DEBUG
#undef SINGULARITY_ROOT_VERBOSE
#undef SINGULARITY_MY_SIGN

#endif // _SINGULARITY_EOS_UTILS_ROOT_FINDING_HPP_
