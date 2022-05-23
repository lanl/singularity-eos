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

#include <sstream>
#include <limits>

#include <root-finding-1d/root_finding.hpp>
#include <singularity-eos/eos/eos.hpp>

namespace singularity {

PORTABLE_INLINE_FUNCTION Real square(const Real x) { return x * x; }
PORTABLE_INLINE_FUNCTION Real cube(const Real x) { return x * x * x; }
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
/*
The Gruneisen EOS diverges at a specific compression. Ensure that the maximum density is
below the smallest singularity in the reference pressure curve
*/
PORTABLE_FUNCTION Real Gruneisen::ComputeRhoMax(const Real s1, const Real s2, const Real s3,
                                            const Real rho0) {
  // Polynomial from the denominator of the reference pressure curve:
  auto poly = [=](Real eta) { return 1 - s1 * eta - s2 * square(eta) - s3 * cube(eta); };

  // First find the eta root. A negative root indicates that there is no maximum density.
  Real root = 0; // Non-sensical root means a root hasn't been found yet
  Real discriminant, root1, root2;
  if (s2 == 0 && s3 == 0 && s1 > 0) {
    // Linear Us-up analytic root
    root = 1 / s1;
  } else if (s3 == 0 && s2 != 0) {
    // Quadratic Us-up
    discriminant = square(s1) + 4 * s2;
    if (discriminant < 0) {
      root = -1; // Imaginary roots, so no limit
    } else {
      root1 = (s1 + std::sqrt(square(s1) + 4 * s2)) / (-2 * s2);
      root2 = (s1 - std::sqrt(square(s1) + 4 * s2)) / (-2 * s2);
      root = find_min_bounded_val(root1, root2, 0, 1);
    }
  } else if (s3 != 0) {
    // Cubic Us-up
    discriminant = -18 * s3 * s2 * s1 - 4 * cube(-s2) + square(s2) * square(s1) -
                   4 * s3 * cube(s1) - 27 * square(s3);
    // Use discriminant to help find roots
    if (discriminant == 0) {
      // Easy analytical roots (probably not the case)
      if (square(s2) == 3 * s3 * s1) {
        // Single root with multiplicity of 3 (inflection point is root)
        // Also guards against divide by zero
        root = -s2 / 3 / s3;
      } else {
        // One root has multiplicity of 2 (local extremum is a root)
        root1 = (-9 * s3 - s2 * s1) / 2 / (square(s2) - 3 * s3 * s1);
        root2 = (-4 * s3 * s2 * s1 - 9 * square(s3) + cube(s2)) /
                (-s3 * (square(s2) - 3 * s3 * s1));
        root = find_min_bounded_val(root1, root2, 0, 1);
      }
    } else {
      // We'll need an iterative method to search for the minimum bounded root
      Real minbound = 0;
      Real maxbound = 1;
      if (discriminant > 0) {
        // Three real roots. We need to use the extrema to ensure we have a proper bracket
        // in which to search for the root.
        const Real extremum1 =
            (2 * s2 + std::sqrt(square(-2 * s2) - 4 * 3 * s3 * s1)) / (-2 * 3 * s3);
        const Real extremum2 =
            (2 * s2 - std::sqrt(square(-2 * s2) - 4 * 3 * s3 * s1)) / (-2 * 3 * s3);
        const Real min_extremum = std::min(extremum1, extremum2);
        const Real max_extremum = std::max(extremum1, extremum2);
        if (s3 < 0) {
          // Cubic is *increasing*
          // Because poly(eta = 0) = 1, the only possible root for an increasing function
          // will lie between the etrema.
          minbound = std::max(min_extremum, minbound);
          maxbound = std::min(std::max(max_extremum, minbound), maxbound);
        } else {
          // Cubic is *decreasing*
          // Because poly(eta = 0) = 1, the only possible root for a decreasing function
          // will lie outside of the extrema. Further, the only possibility for multiple
          // bound roots will occur when the extremum are both positive
          if (min_extremum > 0) {
            maxbound = std::min(min_extremum, maxbound);
          }
        }
      }
      if (poly(minbound) * poly(maxbound) <= 0) {
        // Root is appropriately bounded
        using RootFinding1D::bisect;
        using RootFinding1D::RootCounts;
        using RootFinding1D::Status;
        RootCounts counts;
        const Real xtol = 1.e-08;
        const Real ytol = 1.e-08;
        const Real eta_guess = 0.001;
        auto status =
            bisect(poly, 0., eta_guess, minbound, maxbound, xtol, ytol, root, counts);
        if (status != Status::SUCCESS) {
          // Root finder failed even though the solution was bracketed... this is an error
          std::stringstream errorMessage;
          errorMessage << "Gruneisen initialization: Cubic root find failed. Maximum "
                          "density cannot be"
                       << " automatically determined and must be manually specified."
                       << std::endl;
          EOS_ERROR(errorMessage.str().c_str());
        }
      } else {
        // Root doesn't lie within physical bounds for eta so no maximum density exists
        root = -1.;
      }
    } // Cubic discriminant >= 0
  }   // Linear/Quadratic/Cubic
  // `root` is a value of eta so it should be greater than zero
  if (root > 0) {
    return rho0 / (1 - root);
  } else {
    // No bounded root exists so there is no upper-limit on the compression
    return std::numeric_limits<Real>::infinity();
  }
}
PORTABLE_FUNCTION Real Gruneisen::Gamma(const Real rho_in) const {
  const Real rho = std::min(rho_in, _rho_max);
  return rho < _rho0 ? _G0 : _G0 * _rho0 / rho + _b * (1 - _rho0 / rho);
}
PORTABLE_FUNCTION Real Gruneisen::dPres_drho_e(const Real rho_in, const Real sie) const {
  const Real rho = std::min(rho_in, _rho_max);
  if (rho < _rho0) {
    return square(_C0) + Gamma(rho) * sie;
  } else {
    const Real eta = 1 - _rho0 / rho;
    const Real s = _s1 * eta + _s2 * square(eta) + _s3 * cube(eta);
    const Real ds = _s1 + 2 * _s2 * eta + 3 * _s3 * square(eta);
    const Real deta = _rho0 / square(rho);
    const Real dGam = (_b - _G0) * deta;
    const Real P_H = _P0 + square(_C0) * _rho0 * eta / square(1 - s);
    const Real dP_H =
        square(_C0) * _rho0 / square(1 - s) * deta * (1 + 2 * eta * ds / (1 - s));
    const Real E_H = (P_H + _P0) * eta / _rho0 / 2.;
    const Real dE_H = deta * (P_H + _P0) / _rho0 / 2. + eta / _rho0 / 2 * dP_H;
    return dP_H + Gamma(rho) * (sie - E_H) + rho * dGam * (sie - E_H) -
           rho * Gamma(rho) * dE_H;
  }
}
PORTABLE_FUNCTION Real Gruneisen::InternalEnergyFromDensityTemperature(
    const Real rho_in, const Real temp, Real *lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  return _Cv * (temp - _T0);
}
PORTABLE_FUNCTION Real Gruneisen::TemperatureFromDensityInternalEnergy(
    const Real rho_in, const Real sie, Real *lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  return _T0 + sie / _Cv;
}
PORTABLE_FUNCTION Real Gruneisen::PressureFromDensityInternalEnergy(const Real rho_in,
                                                                    const Real sie,
                                                                    Real *lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  Real P_H;
  Real E_H;
  if (rho >= _rho0) {
    const Real eta = 1 - _rho0 / rho;
    const Real s = _s1 * eta + _s2 * square(eta) + _s3 * cube(eta);
    P_H = _P0 + square(_C0) * _rho0 * eta / square(1 - s);
    E_H = (P_H + _P0) * eta / _rho0 / 2.;
  } else {
    // This isn't thermodynamically consistent but it's widely used for expansion
    P_H = _P0 + square(_C0) * (rho - _rho0);
    E_H = 0.;
  }
  return P_H + Gamma(rho) * rho * (sie - E_H);
}
PORTABLE_FUNCTION Real Gruneisen::SpecificHeatFromDensityInternalEnergy(
    const Real rho_in, const Real sie, Real *lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  return _Cv;
}
PORTABLE_FUNCTION Real Gruneisen::BulkModulusFromDensityInternalEnergy(
    const Real rho_in, const Real sie, Real *lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  // The if statement exists here to avoid the divide by zero
  if (rho < _rho0) {
    return rho * square(_C0) +
           _G0 * (rho * sie + PressureFromDensityInternalEnergy(rho, sie));
  } else {
    const Real dPdr_e = dPres_drho_e(rho, sie);
    const Real dPde_r = rho * Gamma(rho);
    // Thermodynamic identity
    return rho * dPdr_e + PressureFromDensityInternalEnergy(rho, sie) / rho * dPde_r;
  }
}
PORTABLE_FUNCTION
Real Gruneisen::GruneisenParamFromDensityInternalEnergy(const Real rho_in, const Real sie,
                                                        Real *lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  return Gamma(rho);
}
// Below are "unimplemented" routines
PORTABLE_FUNCTION Real Gruneisen::PressureFromDensityTemperature(const Real rho_in,
                                                                 const Real temp,
                                                                 Real *lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  return PressureFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_FUNCTION Real Gruneisen::SpecificHeatFromDensityTemperature(const Real rho_in,
                                                                     const Real temp,
                                                                     Real *lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  return SpecificHeatFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_FUNCTION Real Gruneisen::BulkModulusFromDensityTemperature(const Real rho_in,
                                                                    const Real temp,
                                                                    Real *lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  return BulkModulusFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_FUNCTION
Real Gruneisen::GruneisenParamFromDensityTemperature(const Real rho_in, const Real temp,
                                                     Real *lambda) const {
  const Real rho = std::min(rho_in, _rho_max);
  return Gamma(rho);
}
PORTABLE_FUNCTION void Gruneisen::DensityEnergyFromPressureTemperature(
    const Real press, const Real temp, Real *lambda, Real &rho, Real &sie) const {
  sie = _Cv * (temp - _T0);
  // We have a branch at rho0, so we need to decide, based on our pressure, whether we
  // should be above or below rho0
  Real Pref = PressureFromDensityInternalEnergy(_rho0, sie);
  if (press < Pref) {
    rho = (press - _P0 + _C0 * _C0 * _rho0) / (_C0 * _C0 + _G0 * sie);
  } else { // We are in compression; iterate
    auto residual = [&](const Real r) {
      return press - PressureFromDensityInternalEnergy(r, sie);
    };
    Real rho1 = _rho0, res1 = residual(rho1), slope = _G0 * sie + _C0 * _C0, rho2, res2,
         rhom, resm;
    rho2 = (rho > rho1 + 1e-3) ? rho : rho1 + res1 / slope;
    res2 = residual(rho2);
    for (int i = 0; i < 20; ++i) {
      slope = (rho2 - rho1) / (res2 - res1 + 1.0e-10);
      rhom = rho2 - res2 * slope;
      resm = residual(rhom);
      if (resm / press < 1e-8) break;
      rho1 = rho2;
      res1 = res2;
      rho2 = rhom;
      res2 = resm;
    }
    rho = rhom;
  }
}
PORTABLE_FUNCTION void Gruneisen::FillEos(Real &rho_in, Real &temp, Real &sie, Real &press,
                                          Real &cv, Real &bmod,
                                          const unsigned long output,
                                          Real *lambda) const {
  // The following could be sped up with work!
  const unsigned long input = ~output;
  if (thermalqs::temperature & input && thermalqs::pressure & input) {
    DensityEnergyFromPressureTemperature(press, temp, lambda, rho_in, sie);
  }
  else if (thermalqs::density & output || thermalqs::specific_internal_energy & output) {
    // Error out on density or energy output because they're currently required as inputs
    std::stringstream errorMessage;
    errorMessage << "Gruneisen FillEos: Density and energy are currently required inputs "
                 << "except when pressure and temperature are inputs"
                 << std::endl;
    EOS_ERROR(errorMessage.str().c_str());
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
PORTABLE_FUNCTION
void Gruneisen::ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press,
                                       Real &cv, Real &bmod, Real &dpde, Real &dvdt,
                                       Real *lambda) const {
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
