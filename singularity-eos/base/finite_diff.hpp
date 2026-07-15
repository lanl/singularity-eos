
//------------------------------------------------------------------------------
// © 2026. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_BASE_FINITE_DIFF_HPP_
#define _SINGULARITY_EOS_BASE_FINITE_DIFF_HPP_
#include <optional>
#include <singularity-utils/robust_utils.hpp>

namespace singularity {

namespace finite_diff {

template <typename Func>
Real centralDifference(Func &&f, Real x, std::optional<Real> pert = std::nullopt) {
  //  f(x) is a function that calls some g(x1,x2) with either x1 or x2 fixed
  Real h = pert ? *pert : std::max(std::abs(x) * 1e-6, 1e-12);
  Real f_plus = f(x + h);
  Real f_minus = f(x - h);
  return robust::ratio(f_plus - f_minus, 2 * h);
}

template <typename Func>
Real forwardDifference(Func &&f, Real x, std::optional<Real> pert = std::nullopt) {
  //  f(x) is a function that calls some g(x1,x2) with either x1 or x2 fixed
  Real h = pert ? *pert : std::max(std::abs(x) * 1e-6, 1e-12);
  Real f0 = f(x);
  Real fp = f(x + h);
  Real fpp = f(x + 2 * h);
  return robust::ratio(-3 * f0 + 4 * fp - fpp, 2 * h);
}

template <typename Func>
Real forwardDifferenceOrder1(Func &&f, Real x, std::optional<Real> pert = std::nullopt) {
  //  f(x) is a function that calls some g(x1,x2) with either x1 or x2 fixed
  Real h = pert ? *pert : std::max(std::abs(x) * 1e-6, 1e-12);
  Real f0 = f(x);
  Real fp = f(x + h);
  return robust::ratio(fp - f0, h);
}

template <typename Func>
Real backwardDifference(Func &&f, Real x, std::optional<Real> pert = std::nullopt) {
  //  f(x) is a function that calls some g(x1,x2) with either x1 or x2 fixed
  Real h = pert ? *pert : std::max(std::abs(x) * 1e-6, 1e-12);
  Real f0 = f(x);
  Real fm = f(x - h);
  Real fmm = f(x - 2 * h);
  return robust::ratio(3 * f0 - 4 * fm + fmm, 2 * h);
}

template <typename Func>
Real backwardDifferenceOrder1(Func &&f, Real x, std::optional<Real> pert = std::nullopt) {
  //  f(x) is a function that calls some g(x1,x2) with either x1 or x2 fixed
  Real h = pert ? *pert : std::max(std::abs(x) * 1e-6, 1e-12);
  Real f0 = f(x);
  Real fm = f(x - h);
  return robust::ratio(f0 - fm, h);
}

template <typename Func>
Real finiteDifference(Func &&f, Real x, Real xmin, Real xmax,
                      std::optional<Real> pert = std::nullopt) {
  Real h = pert ? *pert : std::max(std::abs(x) * 1e-6, 1e-12);
  const Real left = x - xmin;
  const Real right = xmax - x;

  if (left >= h && right >= h) {
    return centralDifference(f, x, h);
  } else if (right >= 2 * h) {
    return forwardDifferenceOrder1(f, x, h);
  } else {
    return backwardDifferenceOrder1(f, x, h);
  }
}

} // end namespace finite_diff
} // end namespace singularity
#endif // SINGULARITY_EOS_BASE_FINITE_DIFF_HPP_
