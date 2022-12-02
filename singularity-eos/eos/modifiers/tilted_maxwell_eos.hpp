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

#ifndef _SINGULARITY_EOS_EOS_TILTED_MAXWELL_EOS_
#define _SINGULARITY_EOS_EOS_TILTED_MAXWELL_EOS_

#include "stdio.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <tuple>
#include <utility>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/eos_error.hpp>
#include <singularity-eos/base/math_utils.hpp>
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {

using namespace eos_base;
using namespace RootFinding1D;

template <typename T>
class TiltedMaxwellEOS : public EosBase<TiltedMaxwellEOS<T>> {
 public:
  // Generic functions provided by the base class. These contain
  // e.g. the vector overloads that use the scalar versions declared
  // here We explicitly list, rather than using the macro because we
  // overload some methods.

  // TODO(JMM): The modifier EOS's should probably call the specific
  // sub-functions of the class they modify so that they can leverage,
  // e.g., an especially performant or special version of these
  using EosBase<TiltedMaxwellEOS<T>>::TemperatureFromDensityInternalEnergy;
  using EosBase<TiltedMaxwellEOS<T>>::InternalEnergyFromDensityTemperature;
  using EosBase<TiltedMaxwellEOS<T>>::PressureFromDensityTemperature;
  using EosBase<TiltedMaxwellEOS<T>>::PressureFromDensityInternalEnergy;
  using EosBase<TiltedMaxwellEOS<T>>::SpecificHeatFromDensityTemperature;
  using EosBase<TiltedMaxwellEOS<T>>::SpecificHeatFromDensityInternalEnergy;
  using EosBase<TiltedMaxwellEOS<T>>::BulkModulusFromDensityTemperature;
  using EosBase<TiltedMaxwellEOS<T>>::BulkModulusFromDensityInternalEnergy;
  using EosBase<TiltedMaxwellEOS<T>>::GruneisenParamFromDensityTemperature;
  using EosBase<TiltedMaxwellEOS<T>>::GruneisenParamFromDensityInternalEnergy;
  using EosBase<TiltedMaxwellEOS<T>>::PTofRE;
  using EosBase<TiltedMaxwellEOS<T>>::FillEos;

  using BaseType = T;

  // give me std::format or fmt::format...
  static std::string EosType() {
    return std::string("TiltedMaxwellEOS<") + T::EosType() + std::string(">");
  }

  static std::string EosPyType() { return std::string("Shifted") + T::EosPyType(); }

  // move semantics ensures dynamic memory comes along for the ride
  TiltedMaxwellEOS(T &&t, Real tilt_fac) : t_(std::forward<T>(t)), tilt_(tilt_fac) {}
  TiltedMaxwellEOS() = default;

  Real g(Real u) const { return (1.0 - tilt_) * u + tilt_ * u * u;} 
  Real gp(Real u) const { return 1.0 - tilt_ + 2 * tilt_ * u;} 
  
  const Real tilt_;
  static constexpr Real tol = 1.e-10;
  static constexpr Real eos_rho_min = 1.e-8;
  static constexpr Real eos_rho_max = 1.e2;

  auto TestMaxwellConditions(Real pin, Real temp, Real rho_p_min, Real rho_p_max) const { 
    // Find densities above and below spinodal region with p(rho_l) = pin and p(rho_h) = pin  
    Real rho_l, rho_h;
    RootCounts counts_l, counts_h;
    auto stat_l = regula_falsi([this, temp](Real rho){return this->t_.PressureFromDensityTemperature(rho, temp);}, 
                               pin, 0.8 * rho_p_max, 
                               eos_rho_min, rho_p_max,
                               tol, tol,
                               rho_l, counts_l);
    auto stat_h = regula_falsi([this, temp](Real rho){return this->t_.PressureFromDensityTemperature(rho, temp);}, 
                               pin, 1.1 * rho_p_min, 
                               rho_p_min, eos_rho_max,
                               tol, tol,
                               rho_h, counts_h);
    // Integrate Maxwell construction integral (i.e. the chemical potential difference) 
    Real maxwell_int = math_utils::GLQuad32([this, temp, pin](Real rho){return (this->t_.PressureFromDensityTemperature(rho, temp) - pin)/(rho*rho);}, rho_l, rho_h);
    return std::make_tuple(maxwell_int, rho_l, rho_h);
  }
  
  auto GetMaxwellConstruction(Real temp) const { 
    constexpr double eps = std::sqrt(std::numeric_limits<Real>::epsilon());
    
    // Minimize dPdrho to split isotherm into two regions, if minimum is greater than zero 
    // there is no spinodal region
    auto dpdrho = [this, eps, temp](double rho) { 
      return (this->t_.PressureFromDensityTemperature(rho * (1.0 + eps), temp) 
      - this->t_.PressureFromDensityTemperature(rho * (1.0 - eps), temp)) / (2 * eps * rho);
    };
    auto d2pdrho2 = [this, eps, temp](double rho) { 
      return (this->t_.PressureFromDensityTemperature(rho * (1.0 + eps), temp) 
            - 2*this->t_.PressureFromDensityTemperature(rho, temp)
            + this->t_.PressureFromDensityTemperature(rho * (1.0 - eps), temp)) / (eps * eps * rho * rho);
    };
    Real rho_dpdrho_min = 0.0; 
    RootCounts counts;
    auto stat = regula_falsi(d2pdrho2, 
                             0.0, 1.0, 
                             0.1, 3.0,
                             tol, tol,
                             rho_dpdrho_min, counts);
    if (dpdrho(rho_dpdrho_min) >= 0.0) return std::make_tuple(0.0, 0.0, 0.0, false); 

    // Find spinodal region by finding maximum pressure below and minimum pressure above rho_dpdrho_min
    Real rho_p_max = 0.0; 
    Real rho_p_min = 0.0;
    stat = regula_falsi(dpdrho, 
                        0.0, 0.8 * rho_dpdrho_min, 
                        eos_rho_min, rho_dpdrho_min, 
                        tol, tol,
                        rho_p_max, counts);
    stat = regula_falsi(dpdrho, 
                        0.0, 1.2 * rho_dpdrho_min, 
                        rho_dpdrho_min, eos_rho_max, 
                        tol, tol,
                        rho_p_min, counts);
    Real p_max = std::min(t_.PressureFromDensityTemperature(rho_p_max, temp), 
                          t_.PressureFromDensityTemperature(eos_rho_max, temp));
    Real p_min = std::max(t_.PressureFromDensityTemperature(rho_p_min, temp),
                          t_.PressureFromDensityTemperature(eos_rho_min, temp));

    // Now need to search in between p_min and p_max for the pressure obeying the Maxwell construction conditions 
    // This requires a function that takes a p_in, finds a rho_low with p_in = p(rho_low, temp) and a rho_hi 
    // with p_in = p(rho_hi, temp), then returns the Maxwell construction integral between those two densities.
    Real p_maxwell;
    stat = regula_falsi([this, temp, rho_p_max, rho_p_min](double p) -> double {return std::get<0>(this->TestMaxwellConditions(p, temp, rho_p_min, rho_p_max));},
                        0.0, 0.5 * (p_max + p_min), 
                        p_min, p_max,
                        tol, tol, 
                        p_maxwell, counts);
    const auto [f, rho_maxwell_low, rho_maxwell_hi] = TestMaxwellConditions(p_maxwell, temp, rho_p_min, rho_p_max); 

    return std::make_tuple(rho_maxwell_low, rho_maxwell_hi, p_maxwell, true);
  }

  auto GetTiltedMaxwellConstruction(Real temp) const {
    const auto [rho_maxwell_low, rho_maxwell_hi, p_maxwell, has_spinodal] = this->GetMaxwellConstruction(temp);
    if (!has_spinodal) return std::make_tuple(0.0, 0.0, 0.0, false); 
    
    auto func = [this, temp, p_maxwell](std::array<Real, 2> x){
      const Real Pl = this->t_.PressureFromDensityTemperature(x[0], temp); 
      const Real Ph = this->t_.PressureFromDensityTemperature(x[1], temp);
      const Real PD = math_utils::GLQuad32([this, temp](Real rho){return this->t_.PressureFromDensityTemperature(rho, temp)/(rho * rho);}, x[0], x[1])/(1/x[0] - 1/x[1]);
      return std::array<Real, 2>{(Pl - PD * this->gp(1.0)) / p_maxwell, (Ph - PD * this->gp(0.0)) / p_maxwell};
    };
    
    // NR Iteration
    constexpr double eps = std::sqrt(std::numeric_limits<Real>::epsilon());
    constexpr int max_iters = 50;
    std::array<Real, 2> x{1.e-1 * rho_maxwell_low, 1.5 * rho_maxwell_hi};
    int iter;
    for (iter = 0; iter < max_iters; ++iter) { 
      // Build the Jacobian 
      Real J[2][2]; 
      for (int iv = 0; iv < 2; ++iv) {
        std::array<Real, 2>xp{x};
        std::array<Real, 2>xm{x};
        xp[iv] *= (1.0 + eps); 
        xm[iv] *= (1.0 - eps);  
        const auto dfp = func(xp);
        const auto dfm = func(xm); 
        for (int fi = 0; fi < 2; ++fi) {
          J[iv][fi] = (dfp[fi] - dfm[fi]) / (xp[iv] - xm[iv]); 
        }
      }

      // Invert the jacobian and update 
      auto fc = func(x); 
      const Real det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
      const std::array<Real, 2> dx{-(fc[0]*J[1][1] - fc[1]*J[1][0]) / det,
                                   -(fc[1]*J[0][0] - fc[0]*J[0][1]) / det};

      const double alpha = std::min(std::min(0.1 * std::abs(x[0] / dx[0]), 1.0), 0.1 * std::abs(x[1]/dx[1]));
      x[0] += alpha * dx[0];
      x[1] += alpha * dx[1];
      // Check convergence 
      if (std::abs(fc[0]) < tol && std::abs(fc[1]) < tol) {
        printf("iter = %i fc[0] = %e fc[1] = %e\n", iter, fc[0], fc[1]);
        break;
      }
    } 
    const Real PD = math_utils::GLQuad32([this, temp](Real rho) {return this->t_.PressureFromDensityTemperature(rho, temp)/(rho * rho);}, x[0], x[1])/(1/x[0] - 1/x[1]);
    return std::make_tuple(x[0], x[1], PD, true); 
  } 

  auto GetOnDevice() { return TiltedMaxwellEOS<T>(t_.GetOnDevice()); }
  inline void Finalize() { t_.Finalize(); }

  PORTABLE_FUNCTION
  Real TemperatureFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const {
    return t_.TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  }
  PORTABLE_FUNCTION
  Real InternalEnergyFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const {
    Real energy = t_.InternalEnergyFromDensityTemperature(rho, temperature, lambda);
    return energy;
  }
  PORTABLE_FUNCTION
  Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                         Real *lambda = nullptr) const {
    return t_.PressureFromDensityInternalEnergy(rho, sie, lambda);
  }
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityInternalEnergy(const Real rho, const Real sie,
                                             Real *lambda = nullptr) const {
    return t_.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
  }
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const {
    return t_.BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
  }
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityInternalEnergy(const Real rho, const Real sie,
                                               Real *lambda = nullptr) const {
    return t_.GruneisenParamFromDensityInternalEnergy(rho, sie, lambda);
  }
  PORTABLE_FUNCTION
  Real PressureFromDensityTemperature(const Real rho, const Real temperature,
                                      Real *lambda = nullptr) const {
    return t_.PressureFromDensityTemperature(rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityTemperature(const Real rho, const Real temperature,
                                          Real *lambda = nullptr) const {
    return t_.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityTemperature(const Real rho, const Real temperature,
                                         Real *lambda = nullptr) const {
    return t_.BulkModulusFromDensityTemperature(rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const {
    return t_.GruneisenParamFromDensityTemperature(rho, temperature, lambda);
  }
  PORTABLE_FUNCTION
  void FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
               const unsigned long output, Real *lambda = nullptr) const {
    Real senergy;
    switch (t_.PreferredInput()) {
    case thermalqs::density | thermalqs::temperature:
      t_.FillEos(rho, temp, energy, press, cv, bmod, output, lambda);
      break;
    case thermalqs::density | thermalqs::specific_internal_energy:
      senergy = energy;
      t_.FillEos(rho, temp, senergy, press, cv, bmod, output, lambda);
      break;
    default:
      EOS_ERROR("Didn't find a valid input for TiltedMaxwellEOS::FillEOS\n");
    }
  }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return t_.nlambda(); }

  static constexpr unsigned long PreferredInput() { return T::PreferredInput(); }

  PORTABLE_FUNCTION void PrintParams() const {
    t_.PrintParams();
  }
  PORTABLE_FUNCTION
  void DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                            Real *lambda, Real &rho, Real &sie) const {
    t_.DensityEnergyFromPressureTemperature(press, temp, lambda, rho, sie);
    sie = sie;
  }

  PORTABLE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const {
    t_.ValuesAtReferenceState(rho, temp, sie, press, cv, bmod, dpde, dvdt, lambda);
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumDensity() const {
    return t_.MinimumDensity();
  }
  PORTABLE_FORCEINLINE_FUNCTION Real MinimumTemperature() const {
    return t_.MinimumTemperature();
  }

 private:
  T t_;
};

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_SHIFTED_EOS_
