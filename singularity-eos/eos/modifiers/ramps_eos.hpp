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

#ifndef _SINGULARITY_EOS_EOS_RAMPS_EOS_
#define _SINGULARITY_EOS_EOS_RAMPS_EOS_

#include "stdio.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <utility>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/eos_error.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {

using namespace eos_base;

using singularity::robust::ratio;

template <typename T>
void pAlpha2BilinearRampParams(const T &eos, const Real alpha0, const Real Pe,
                               const Real Pc, Real &r0, Real &a, Real &b, Real &c) {
  // get reference conditions
  Real rho0, T0, sie0, P0, cv0, bmod0, dpde0, dvdt0;
  Real rmid, r1;
  eos.ValuesAtReferenceState(rho0, T0, sie0, P0, cv0, bmod0, dpde0, dvdt0);
  // calculate r0, ensure alpha0 > 1
  assert(alpha0 > 1.0);
  r0 = ratio(rho0, alpha0);
  // calculate rmid
  auto rmid_func = PORTABLE_LAMBDA(const Real x) {
    return eos.PressureFromDensityTemperature(alpha0 * x, T0);
  };
  RootFinding1D::RootCounts co{};
  // get upper bound to density informed by the reference
  // bulk modulus
  const Real max_exp_arg = std::log(std::numeric_limits<Real>::max() * 0.99);
  const Real exp_arg = std::min(max_exp_arg, ratio((2.0 * Pc - P0), bmod0));
  const Real rho_ub = rho0 * std::exp(exp_arg);
  // finds where rmid_func = Pe
  RootFinding1D::findRoot(rmid_func, Pe, rho0, r0, rho_ub, 1.e-12, 1.e-12, rmid, co);
  // calculate r1
  auto r1_func = PORTABLE_LAMBDA(const Real x) {
    return eos.PressureFromDensityTemperature(x, T0);
  };
  // finds where r1_func = Pc
  RootFinding1D::findRoot(r1_func, Pc, rmid, r0, rho_ub, 1.e-12, 1.e-12, r1, co);
  // a
  a = ratio(r0 * Pe, rmid - r0);
  // b
  b = ratio(r0 * (Pc - Pe), r1 - rmid);
  // c
  c = ratio(Pc * rmid - Pe * r1, r0 * (Pc - Pe));
  return;
}

template <typename T>
class BilinearRampEOS : public EosBase<BilinearRampEOS<T>> {
 public:
  // move semantics ensures dynamic memory comes along for the ride
  BilinearRampEOS(T &&t, const Real r0, const Real a, const Real b, const Real c)
      : t_(std::forward<T>(t)), r0_(r0), a_(a), b_(b), c_(c),
        rmid_(r0 * (a - b * c) / (a - b)), Pmid_(a * (rmid_ / r0 - 1.0)) {
    // add input parameter checks to ensure validity of the ramp
    assert(r0 > 0.0);
    assert(a > 0.0);
    assert(b >= 0);
    assert(a != b);
  }
  BilinearRampEOS() = default;

  using BaseType = T;

  // give me std::format or fmt::format...
  static std::string EosType() {
    return std::string("BilinearRampEOS<") + T::EosType() + std::string(">");
  }

  static std::string EosPyType() { return std::string("BilinearRamp") + T::EosPyType(); }

  auto GetOnDevice() { return BilinearRampEOS<T>(t_.GetOnDevice(), r0_, a_, b_, c_); }
  inline void Finalize() { t_.Finalize(); }

  PORTABLE_INLINE_FUNCTION
  Real get_ramp_pressure(Real rho) const {
    const Real r_r0{ratio(rho, r0_)};
    const Real p_ramp{rho < r0_     ? 0.0
                      : rho < rmid_ ? a_ * (r_r0 - 1.0)
                                    : b_ * (r_r0 - c_)};
    return p_ramp;
  }

  PORTABLE_INLINE_FUNCTION
  Real get_ramp_density(Real P) const {
    const Real rho_ramp{P < Pmid_ ? r0_ * (ratio(P, a_) + 1.0)
                                  : r0_ * (ratio(P, b_) + c_)};
    return rho_ramp;
  }

  PORTABLE_INLINE_FUNCTION
  Real get_ramp_dpdrho(Real rho) const {
    const Real dpdr{rho < r0_ ? 0.0 : rho < rmid_ ? ratio(a_, r0_) : ratio(b_, r0_)};
    return dpdr;
  }

  PORTABLE_FUNCTION
  Real TemperatureFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const {
    return t_.TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  }
  PORTABLE_FUNCTION
  Real InternalEnergyFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const {
    return t_.InternalEnergyFromDensityTemperature(rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                         Real *lambda = nullptr) const {
    // ramp pressure
    const Real p_ramp{get_ramp_pressure(rho)};
    // pressure from eos
    const Real p_eos{t_.PressureFromDensityInternalEnergy(rho, sie, lambda)};
    // return max(p_ramp, p_eos)
    return p_eos < p_ramp ? p_ramp : p_eos;
  }
  PORTABLE_FUNCTION
  Real EntropyFromDensityInternalEnergy(const Real rho, const Real sie,
                                         Real *lambda = nullptr) const {
    return t_.EntropyFromDensityInternalEnergy(rho, sie, lambda);
  }
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityInternalEnergy(const Real rho, const Real sie,
                                             Real *lambda = nullptr) const {
    return t_.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
  }
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const {
    const Real p_ramp{get_ramp_pressure(rho)};
    const Real p_eos{t_.PressureFromDensityInternalEnergy(rho, sie, lambda)};

    return p_eos < p_ramp ? rho * get_ramp_dpdrho(rho)
                          : t_.BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
  }
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityInternalEnergy(const Real rho, const Real sie,
                                               Real *lambda = nullptr) const {
    return t_.GruneisenParamFromDensityInternalEnergy(rho, sie, lambda);
  }
  PORTABLE_FUNCTION
  Real PressureFromDensityTemperature(const Real rho, const Real temperature,
                                      Real *lambda = nullptr) const {
    const Real p_ramp{get_ramp_pressure(rho)};
    const Real p_eos{t_.PressureFromDensityTemperature(rho, temperature, lambda)};
    return p_eos < p_ramp ? p_ramp : p_eos;
  }
  PORTABLE_FUNCTION
  Real EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                         Real *lambda = nullptr) const {
    return t_.EntropyFromDensityTemperature(rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityTemperature(const Real rho, const Real temperature,
                                          Real *lambda = nullptr) const {
    return t_.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityTemperature(const Real rho, const Real temperature,
                                         Real *lambda = nullptr) const {
    const Real p_ramp{get_ramp_pressure(rho)};
    const Real p_eos{t_.PressureFromDensityTemperature(rho, temperature, lambda)};
    return p_eos < p_ramp
               ? rho * get_ramp_dpdrho(rho)
               : t_.BulkModulusFromDensityTemperature(rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const {
    return t_.GruneisenParamFromDensityTemperature(rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  void FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
               const unsigned long output, Real *lambda = nullptr) const {
    // output passed into internal filleos can't include pressure
    const unsigned long ramp_out = output & ~thermalqs::pressure;
    // if pressure is output, calculate it first
    if (output & thermalqs::pressure) {
      // maybe switch case on preferred input and check for not output of the input
      // for now sie lookup is prioritized
      if (!(output & thermalqs::specific_internal_energy)) {
        press = PressureFromDensityInternalEnergy(rho, energy, lambda);
      } else if (!(output & thermalqs::temperature)) {
        press = PressureFromDensityTemperature(rho, temp, lambda);
      }
    }
    // call internal filleos
    t_.FillEos(rho, temp, energy, press, cv, bmod, ramp_out, lambda);
    // fill ramp density
    Real rho_ramp{rho};
    if (output & thermalqs::density) {
      assert(!(output & thermalqs::pressure));
      rho_ramp = get_ramp_density(press);
    }
    rho = rho_ramp < rho ? rho_ramp : rho;
    // bulk modulus
    bmod = BulkModulusFromDensityInternalEnergy(rho, energy, lambda);
    return;
  }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return t_.nlambda(); }

  static constexpr unsigned long PreferredInput() { return T::PreferredInput(); }

  PORTABLE_FUNCTION void PrintParams() const {
    t_.PrintParams();
    printf("Ramp Params:\nr0=%.14e\na=%.14e\nb=%.14e\nc=%.14e\n", r0_, a_, b_, c_);
  }
  PORTABLE_FUNCTION
  void DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                            Real *lambda, Real &rho, Real &sie) const {
    t_.DensityEnergyFromPressureTemperature(press, temp, lambda, rho, sie);
  }

  PORTABLE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const {
    t_.ValuesAtReferenceState(rho, temp, sie, press, cv, bmod, dpde, dvdt, lambda);
  }

  // Vector functions that overload the scalar versions declared here.
  SG_ADD_BASE_CLASS_USINGS(BilinearRampEOS<T>)

 private:
  T t_;
  Real r0_;
  Real a_;
  Real b_;
  Real c_;
  Real rmid_;
  Real Pmid_;
};

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_RAMPS_EOS_
