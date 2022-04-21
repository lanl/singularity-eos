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
#include "../ports-of-call/portability.hpp"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define SINGULARITY_ROOT_DEBUG (0)
#define SINGULARITY_ROOT_VERBOSE (0)
// TODO: mauneyc: this isn't used here, and I can't
// find it used elsewhere. it also doesn't get unset.
// should it be ditched?
#define SINGULARITY_ROOT_NAN_OK (0)

#define SINGULARITY_MY_SIGN(x) (x > 0) - (x < 0)

namespace RootFinding1D {
constexpr const int SECANT_NITER_MAX{1000};
constexpr const int BISECT_NITER_MAX{1000};
constexpr const int BISECT_REG_MAX{1000};
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
  static constexpr int nbins_{15};
  mutable Real counts_[nbins_];

 public:
  PORTABLE_INLINE_FUNCTION
  RootCounts() {
    for (int i{0}; i < nbins_; ++i)
      counts_[i] = 0;
  }
  PORTABLE_INLINE_FUNCTION void reset() {
    for (int i{0}; i < nbins_; ++i)
      counts_[i] = 0;
  }
  PORTABLE_INLINE_FUNCTION void increment(int i) const {
    assert(i < nbins_ && i >= 0);
    counts_[i] += 1;
  }
  PORTABLE_INLINE_FUNCTION Real total() const {
    Real tot{1.e-20};
    for (int i{0}; i < nbins_; ++i)
      tot += counts_[i];
    return tot;
  }
  PORTABLE_INLINE_FUNCTION const Real &operator[](const int i) const {
    assert(i < nbins_ && i >= 0);
    return counts_[i];
  }
  PORTABLE_INLINE_FUNCTION Real &operator[](const int i) {
    assert(i < nbins_ && i >= 0);
    return counts_[i];
  }
  PORTABLE_INLINE_FUNCTION void print_counts() const {
    for (int i{0}; i < nbins_; ++i)
      printf("%e\n", counts_[i]);
  }
  PORTABLE_INLINE_FUNCTION int nBins() const { return nbins_; }
  PORTABLE_INLINE_FUNCTION int more() const { return nbins_ - 1; }
};

// solves for f(x,params) - ytarget = 0
template <typename T>
PORTABLE_INLINE_FUNCTION Status findRoot(const T &f, const Real ytarget, Real xguess,
                                         const Real xmin, const Real xmax,
                                         const Real xtol, const Real ytol, Real &xroot,
                                         const RootCounts &counts) {
  Status status;

  // first check if we're at the max or min values
  // TODO: these short-circuits almost never happen. Should they be removed?
  // ~JMM
  const Real fmax = f(xmax);
  const Real errmax = fabs(fmax - ytarget) / (fabs(fmax) + ytol);
  if (errmax < ytol) {
    xroot = xmax;
    counts.increment(0);
    return Status::SUCCESS;
  }
  const Real fmin = f(xmin);
  const Real errmin = fabs(fmin - ytarget) / (fabs(fmin) + ytol);
  if (errmin < ytol) {
    xroot = xmin;
    counts.increment(0);
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
                                       const RootCounts &counts) {
  Real dx;
  Real x_last, y, yp, ym, dyNum, dyDen, dy;

  Real x = xguess;
  unsigned int iter{0};
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
      counts.increment(counts.more());
      return Status::FAIL;
    }
    if (fabs(x - x_last) / (fabs(x) + xtol) < xtol) break;
  }
  // increment iter
  ++iter;
  if (iter < counts.nBins()) {
    counts.increment(iter);
  } else {
    counts.increment(counts.more());
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
                                       const RootCounts &counts) {
  Real xl, xr, fl, fr, dx;

  Real grow = 0.01;
  Real x = xguess;
  if (fabs(x) < xtol) {
    x += 2. * xtol;
  }
  // do { // Try to find reasonable region for bisection
  for (int i{0}; i < BISECT_REG_MAX; ++i) {
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
        int nx = 300;
        Real dx = (xmax - xmin) / (nx - 1);
        fprintf(stderr, "Area map:\nx\ty\n");
        for (int i = 0; i < nx; i++) {
          fprintf(stderr, "%.4f\t%.4e\n", x + i * dx, f(x + i * dx));
        }
#endif
        return Status::FAIL;
      }
    }
  }

  for (int i{0}; i < BISECT_NITER_MAX; ++i) {
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
