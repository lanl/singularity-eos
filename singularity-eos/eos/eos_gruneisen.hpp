//------------------------------------------------------------------------------
// © 2021-2024. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_EOS_GRUNEISEN_HPP_
#define _SINGULARITY_EOS_EOS_GRUNEISEN_HPP_

#include <cmath>
#include <cstdio>

#include <limits>

#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/math_utils.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {

using namespace eos_base;

// COMMENT: This is meant to be an implementation of the Steinberg version of
// the Gruneisen EOS which should correspond to eostype(3) in xRAGE and
// /[...]/eos/gruneisen in FLAG
class Gruneisen : public EosBase<Gruneisen> {
 public:
  Gruneisen() = default;
  PORTABLE_INLINE_FUNCTION
  Gruneisen(const Real C0, const Real s1, const Real s2, const Real s3, const Real G0,
            const Real b, const Real rho0, const Real T0, const Real P0, const Real Cv,
            const Real rho_max)
      : _C0(C0), _s1(s1), _s2(s2), _s3(s3), _G0(G0), _b(b), _rho0(rho0), _T0(T0), _P0(P0),
        _Cv(Cv), _rho_max(rho_max) {
    // Warn user when provided rho_max is greater than the computed rho_max
#ifndef NDEBUG
    const Real computed_rho_max = ComputeRhoMax(s1, s2, s3, rho0);
    if (rho_max > RHOMAX_SAFETY * computed_rho_max) {
      printf(
          "WARNING: Provided rho_max, %g, is greater than the computed rho_max of %g.\n",
          rho_max, computed_rho_max);
      printf("         States beyond %g g/cm^3 are unphysical (i.e. imaginary sound "
             "speeds).\n",
             computed_rho_max);
    }
#endif
  }
  // Constructor when rho_max isn't specified automatically determines _rho_max
  PORTABLE_INLINE_FUNCTION
  Gruneisen(const Real C0, const Real s1, const Real s2, const Real s3, const Real G0,
            const Real b, const Real rho0, const Real T0, const Real P0, const Real Cv)
      : _C0(C0), _s1(s1), _s2(s2), _s3(s3), _G0(G0), _b(b), _rho0(rho0), _T0(T0), _P0(P0),
        _Cv(Cv), _rho_max(RHOMAX_SAFETY * ComputeRhoMax(s1, s2, s3, rho0)) {}
  static PORTABLE_INLINE_FUNCTION Real ComputeRhoMax(const Real s1, const Real s2,
                                                     const Real s3, const Real rho0);
  PORTABLE_INLINE_FUNCTION Real
  MaxStableDensityAtTemperature(const Real temperature) const;
  Gruneisen GetOnDevice() { return *this; }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _T0 + sie / _Cv;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temp,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _Cv * (temp - _T0);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  PressureFromDensityTemperature(const Real rho, const Real temperature,
                                 Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MinInternalEnergyFromDensity(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperatummmmmmre,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _Cv;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _Cv;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return Gamma(std::min(rho, _rho_max));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return Gamma(std::min(rho, _rho_max));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
          const unsigned long output,
          Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                         Real &bmod, Real &dpde, Real &dvdt,
                         Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  // Generic functions provided by the base class. These contain e.g. the vector
  // overloads that use the scalar versions declared here
  SG_ADD_BASE_CLASS_USINGS(Gruneisen)
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return 0;
  }
  static inline unsigned long max_scratch_size(unsigned int nelements) { return 0; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    static constexpr char s1[]{"Gruneisen Params: "};
    printf("%s C0:%e s1:%e s2:%e s3:%e\n  G0:%e b:%e rho0:%e T0:%e\n  P0:%eCv:%e "
           "rho_max:%e\n",
           s1, _C0, _s1, _s2, _s3, _G0, _b, _rho0, _T0, _P0, _Cv, _rho_max);
  }
  template <typename Indexer_t>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const;
  inline void Finalize() {}
  static std::string EosType() { return std::string("Gruneisen"); }
  static std::string EosPyType() { return EosType(); }

 private:
  Real _C0, _s1, _s2, _s3, _G0, _b, _rho0, _T0, _P0, _Cv, _rho_max;
  // static constexpr const char _eos_type[] = {"Gruneisen"};
  PORTABLE_INLINE_FUNCTION
  Real Gamma(const Real rho_in) const {
    const Real rho = std::min(rho_in, _rho_max);
    return rho < _rho0 ? _G0 : _G0 * _rho0 / rho + _b * (1 - _rho0 / rho);
  }
  PORTABLE_INLINE_FUNCTION
  Real dPres_drho_e(const Real rho, const Real sie) const;
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
  // Scaling factor for density singularity (reference pressure blows up). Consistent with
  // Pagosa implementation
  static constexpr Real RHOMAX_SAFETY = 0.99;
};

namespace gruneisen_utils {
PORTABLE_INLINE_FUNCTION Real find_min_bounded_val(const Real val1, const Real val2,
                                                   const Real min, const Real max) {
  // Try to find the minimum bounded value. If none exists, returns negative infinity
  const Real minval = std::min(val1, val2);
  const Real maxval = std::max(val1, val2);
  if (minval > min && minval < max) {
    // Minimum of two values is bounded so return it
    return minval;
  } else if (maxval > min && maxval < max) {
    // Minimum isn't bounded but maximum is so return that
    return maxval;
  } else {
    // Neither value is bounded so just return negative infinity to clearly show that the
    // value is out of bounds
    return -std::numeric_limits<Real>::infinity();
  }
}
PORTABLE_INLINE_FUNCTION bool is_near_zero(const Real val, const Real tol) {
  return std::abs(val) < tol;
}
} // namespace gruneisen_utils
/*
The Gruneisen EOS diverges at a specific compression. Ensure that the maximum density is
below the smallest singularity in the reference pressure curve
*/
PORTABLE_INLINE_FUNCTION Real Gruneisen::ComputeRhoMax(const Real s1, const Real s2,
                                                       const Real s3, const Real rho0) {
  using namespace gruneisen_utils;
  using namespace math_utils;

  // Polynomial from the denominator of the reference pressure curve:
  auto poly = [=](Real eta) {
    return 1 - s1 * eta - s2 * pow<2>(eta) - s3 * pow<3>(eta);
  };

  // Comparisons to zero are actually made to
  static constexpr Real EPS = robust::EPS();

  // First find the eta root. A negative root indicates that there is no maximum density.
  Real root = -1.; // Non-sensical root means none exists
  Real discriminant, root1, root2;
  if (is_near_zero(s2, EPS) && is_near_zero(s3, EPS) && s1 > 1.) {
    // Linear Us-up analytic root
    root = 1. / s1;
  } else if (is_near_zero(s3, EPS) && !is_near_zero(s2, EPS)) {
    // Quadratic Us-up
    discriminant = pow<2>(s1) + 4. * s2;
    if (discriminant < 0.) {
      root = -1.; // Imaginary roots, so no limit
    } else {
      root1 = (s1 + std::sqrt(discriminant)) / (-2. * s2);
      root2 = (s1 - std::sqrt(discriminant)) / (-2. * s2);
      root = find_min_bounded_val(root1, root2, 0., 1.);
    }
  } else if (!is_near_zero(s3, EPS)) {
    // Cubic Us-up: we'll use an iterative method to search for the minimum bounded root
    // Note: when the discriminant of a cubic is less than zero, only one real root exists
    //       so we don't need anything fancy to find the root.
    discriminant = -18. * s3 * s2 * s1 - 4. * pow<3>(-s2) + pow<2>(s2) * pow<2>(s1) -
                   4. * s3 * pow<3>(s1) - 27. * pow<2>(s3);
    Real minbound = 0.;
    Real maxbound = 1.;
    if (discriminant >= 0) {
      // Three real roots (roots may have multiplicity). We need to use the extrema to
      // ensure we have a proper bracket in which to search for the root (if they exist).
      // Note: If the discriminant is positive, then `pow<2>(s2) - 3 * s3 * s1` must
      //       necessarily be positive (the reverse does not hold). Thus the extrema below
      //       are not imaginary.
      const Real extremum1 = (s2 + std::sqrt(pow<2>(s2) - 3. * s3 * s1)) / (-3. * s3);
      const Real extremum2 = (s2 - std::sqrt(pow<2>(s2) - 3. * s3 * s1)) / (-3. * s3);
      const Real min_extremum = std::min(extremum1, extremum2);
      const Real max_extremum = std::max(extremum1, extremum2);
      // Because poly(eta = 0) = 1, the only possible root will lie in an area of the
      // cubic that has negative slope. That area is determined by the extrema.
      if (s3 < 0.) {
        // Cubic is *increasing*
        // The only place with negative slope is between the extrema so that's where the
        // relevant root will be.
        minbound = std::max(min_extremum, minbound);
        maxbound = std::min(std::max(max_extremum, minbound), maxbound);
      } else {
        // Cubic is *decreasing*
        // The only place with negative slope is outside the extrema. Further, the only
        // possibility for multiple bound roots will occur when the extremum are both
        // positive, this is the only case we need to focus on.
        if (min_extremum > 0) {
          maxbound = std::min(min_extremum, maxbound);
        }
      } // s3 > 0 ; else
    }   // discriminant >= 0
    if (poly(minbound) * poly(maxbound) < 0.) {
      // Root is appropriately bounded
      using RootFinding1D::regula_falsi;
      using RootFinding1D::RootCounts;
      using RootFinding1D::Status;
      static constexpr Real factor = 100;
      static constexpr Real xtol = factor * EPS;
      static constexpr Real ytol = factor / 10 * EPS;
      static constexpr Real eta_guess = 0.001;
      auto status =
          regula_falsi(poly, 0., eta_guess, minbound, maxbound, xtol, ytol, root);
      if (status != Status::SUCCESS) {
        // Root finder failed even though the solution was bracketed... this is an error
        EOS_ERROR("Gruneisen initialization: Cubic root find failed. Maximum "
                  "density cannot be"
                  " automatically determined and must be manually specified\n.");
      }
    } else if (is_near_zero(poly(minbound), EPS)) {
      root = minbound;
    } else if (is_near_zero(poly(maxbound), EPS)) {
      root = maxbound;
    } else {
      // Root doesn't lie within physical bounds for eta so no maximum density exists
      root = -1.;
    }
  } // Linear/Quadratic/Cubic
  // `root` is a value of eta so it should be greater than zero and less than 1
  if (root > 0. && root < 1.) {
    return rho0 / ((1 - root) + EPS); // convert from eta to compression
  } else {
    // No bounded root exists so there is no upper-limit on the compression
    return std::numeric_limits<Real>::infinity();
  }
}

PORTABLE_INLINE_FUNCTION Real Gruneisen::dPres_drho_e(const Real rho_in,
                                                      const Real sie) const {
  using namespace math_utils;
  const Real rho = std::min(rho_in, _rho_max);
  if (rho < _rho0) {
    return pow<2>(_C0) + Gamma(rho) * sie;
  } else {
    const Real eta = 1 - _rho0 / rho;
    const Real s = _s1 * eta + _s2 * pow<2>(eta) + _s3 * pow<3>(eta);
    const Real ds = _s1 + 2 * _s2 * eta + 3 * _s3 * pow<2>(eta);
    const Real deta = _rho0 / pow<2>(rho);
    const Real dGam = (_b - _G0) * deta;
    const Real P_H = _P0 + pow<2>(_C0) * _rho0 * eta / pow<2>(1 - s);
    const Real dP_H = math_utils::pow<2>(_C0) * _rho0 / pow<2>(1 - s) * deta *
                      (1 + 2 * eta * ds / (1 - s));
    const Real E_H = (P_H + _P0) * eta / _rho0 / 2.;
    const Real dE_H = deta * (P_H + _P0) / _rho0 / 2. + eta / _rho0 / 2 * dP_H;
    return dP_H + Gamma(rho) * (sie - E_H) + rho * dGam * (sie - E_H) -
           rho * Gamma(rho) * dE_H;
  }
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real Gruneisen::PressureFromDensityInternalEnergy(
    const Real rho_in, const Real sie, Indexer_t &&lambda) const {
  using namespace math_utils;
  const Real rho = std::min(rho_in, _rho_max);
  Real P_H;
  Real E_H;
  if (rho >= _rho0) {
    const Real eta = 1 - _rho0 / rho;
    const Real s = _s1 * eta + _s2 * pow<2>(eta) + _s3 * pow<3>(eta);
    P_H = _P0 + pow<2>(_C0) * _rho0 * eta / pow<2>(1 - s);
    E_H = (P_H + _P0) * eta / _rho0 / 2.;
  } else {
    // This isn't thermodynamically consistent but it's widely used for expansion
    P_H = _P0 + pow<2>(_C0) * (rho - _rho0);
    E_H = 0.;
  }
  return P_H + Gamma(rho) * rho * (sie - E_H);
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real
Gruneisen::MinInternalEnergyFromDensity(const Real rho_in, Indexer_t &&lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  MinInternalEnergyIsNotEnabled("Gruneisen");
  return 0.0;
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real Gruneisen::EntropyFromDensityInternalEnergy(
    const Real rho_in, const Real sie, Indexer_t &&lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  EntropyIsNotEnabled("Gruneisen");
  return 1.0;
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real Gruneisen::BulkModulusFromDensityInternalEnergy(
    const Real rho_in, const Real sie, Indexer_t &&lambda) const {
  using namespace gruneisen_utils;
  const Real rho = std::min(rho_in, _rho_max);
  // The if statement exists here to avoid the divide by zero
  if (rho < _rho0) {
    return rho * math_utils::pow<2>(_C0) +
           _G0 * (rho * sie + PressureFromDensityInternalEnergy(rho, sie));
  } else {
    const Real dPdr_e = dPres_drho_e(rho, sie);
    const Real dPde_r = rho * Gamma(rho);
    // Thermodynamic identity
    return rho * dPdr_e + PressureFromDensityInternalEnergy(rho, sie) / rho * dPde_r;
  }
}
// Below are "unimplemented" routines
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real Gruneisen::PressureFromDensityTemperature(
    const Real rho_in, const Real temp, Indexer_t &&lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  return PressureFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real Gruneisen::EntropyFromDensityTemperature(
    const Real rho_in, const Real temp, Indexer_t &&lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  const Real sie = InternalEnergyFromDensityTemperature(rho, temp);
  return EntropyFromDensityInternalEnergy(rho, sie);
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real Gruneisen::BulkModulusFromDensityTemperature(
    const Real rho_in, const Real temp, Indexer_t &&lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  return BulkModulusFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_INLINE_FUNCTION Real
Gruneisen::MaxStableDensityAtTemperature(const Real temperature) const {
  // Because of the constant cold curve assumption, this EOS is thermodynamically
  // inconsistent, and leads to states that are effectively off the EOS surface even
  // though the temperature is positive. Mathematically, this means that the \Gamma\rho(e
  // - e_H) term becomes highly negative and dominates the positive P_H term at some
  // density. This results in a local maximum in the pressure as a function of density for
  // a given temperature. Beyond this point, the pressure decreases with increasing
  // density, which is thermodynamcially unstable. In a thermodynamically consistent EOS,
  // the cold curve energy would increase with density, leading an appropriate bounding at
  // T=0.

  // Since E_H and P_H are monotonic up to the singularity given by _rho_max, if the
  // derivative of the pressure is negative at _rho_max, a maximum exists and this should
  // be the highest density for the isotherm. If this maximum doesn't exist, then the EOS
  // is thermodynamically consistent up to the maximum density for this isotherm.
  Real slope_at_max_density =
      dPres_drho_e(_rho_max, InternalEnergyFromDensityTemperature(_rho_max, temperature));
  if (slope_at_max_density >= 0) {
    // No maximum pressure before _rho_max
    return _rho_max;
  }

  // Maximum pressure should exist... do a root find to locate where the derivative is
  // zero
  auto dPdrho_T = PORTABLE_LAMBDA(const Real r) {
    return dPres_drho_e(r, InternalEnergyFromDensityTemperature(r, temperature));
  };
  Real rho_lower = 0.9 * _rho0;
  Real rho_upper = std::min(_rho_max, 1.0e4);
  Real rho_at_max_P;
  using RootFinding1D::regula_falsi;
  using RootFinding1D::Status;
  auto status = regula_falsi(dPdrho_T, 0., _rho0, rho_lower, rho_upper, 1.0e-8, 1.0e-8,
                             rho_at_max_P);
  if (status != Status::SUCCESS) {
    // Root finder failed even though the solution should be bracketed
    EOS_ERROR("Gruneisen::MaxStableDensityAtTemperature: "
              "Root find failed to find maximum P at T despite aparent bracket\n");
  }
  return rho_at_max_P;
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void Gruneisen::DensityEnergyFromPressureTemperature(
    const Real press, const Real temp, Indexer_t &&lambda, Real &rho, Real &sie) const {
  // We have a branch at rho0, so we need to decide, based on our pressure, whether we
  // should be above or below rho0
  Real Pref = PressureFromDensityTemperature(_rho0, temp);
  Real rho_lower;
  Real rho_upper;
  // Pick bounds appropriate depending on whether in compression or expansion
  if (press < Pref) {
    rho_lower = 0.;
    rho_upper = 1.1 * _rho0;
  } else {
    rho_lower = 0.9 * _rho0;
    // Find maximum thermodynamically _stable_ density at this temperature. Use this to
    // check if we're actually on the EOS surface or not
    rho_upper = MaxStableDensityAtTemperature(temp);
    auto pres_max = PressureFromDensityTemperature(rho_upper, temp);
    if (press > pres_max) {
      // We're off the EOS surface
      using PortsOfCall::printf;
      printf("Requested pressure, %.15g, exceeds maximum, %.15g, for temperature, %.15g",
             press, pres_max, temp);
      PORTABLE_ALWAYS_THROW_OR_ABORT("Input pressure is off EOS surface");
    }
  }
  auto PofRatT = PORTABLE_LAMBDA(const Real r) {
    return PressureFromDensityTemperature(r, temp);
  };
  using RootFinding1D::regula_falsi;
  using RootFinding1D::Status;
  auto status =
      regula_falsi(PofRatT, press, _rho0, rho_lower, rho_upper, 1.0e-8, 1.0e-8, rho);
  if (status != Status::SUCCESS) {
    // Root finder failed even though the solution was bracketed... this is an error
    EOS_ERROR("Gruneisen::DensityEnergyFromPressureTemperature: "
              "Root find failed to find a solution given P, T\n");
  }
  sie = InternalEnergyFromDensityTemperature(rho, temp);
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void
Gruneisen::FillEos(Real &rho_in, Real &temp, Real &sie, Real &press, Real &cv, Real &bmod,
                   const unsigned long output, Indexer_t &&lambda) const {
  // The following could be sped up with work!
  const unsigned long input = ~output;
  if (thermalqs::temperature & input && thermalqs::pressure & input) {
    DensityEnergyFromPressureTemperature(press, temp, lambda, rho_in, sie);
  } else if (thermalqs::density & output ||
             thermalqs::specific_internal_energy & output) {
    // Error out on density or energy output because they're currently required as inputs
    EOS_ERROR("Gruneisen FillEos: Density and energy are currently required inputs "
              "except when pressure and temperature are inputs.\n");
  }
  const Real rho = std::min(rho_in, _rho_max);
  if (output & thermalqs::pressure) press = PressureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::temperature)
    temp = TemperatureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::bulk_modulus)
    bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::specific_heat)
    cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
}

// TODO(JMM): pre-cache these rather than recomputing them each time
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void
Gruneisen::ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                                  Real &bmod, Real &dpde, Real &dvdt,
                                  Indexer_t &&lambda) const {
  rho = _rho0;
  temp = _T0;
  sie = 0;
  press = PressureFromDensityInternalEnergy(rho, sie, lambda);
  cv = SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
  bmod = BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
  dpde = PressureFromDensityInternalEnergy(rho, sie, lambda) / sie;
  // TODO: chad please fix this one for me. This is wrong.
  Real gm1 = GruneisenParamFromDensityInternalEnergy(rho, sie, lambda) * rho;
  dvdt = gm1 * cv / bmod;
}
} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_GRUNEISEN_HPP_
