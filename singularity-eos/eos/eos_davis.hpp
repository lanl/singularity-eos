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

#ifndef _SINGULARITY_EOS_EOS_EOS_DAVIS_HPP_
#define _SINGULARITY_EOS_EOS_EOS_DAVIS_HPP_

#include <cmath>
#include <cstdio>

#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/math_utils.hpp>
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {

using namespace eos_base;

class DavisReactants : public EosBase<DavisReactants> {
 public:
  DavisReactants() = default;
  PORTABLE_INLINE_FUNCTION
  DavisReactants(const Real rho0, const Real e0, const Real P0, const Real T0,
                 const Real A, const Real B, const Real C, const Real G0, const Real Z,
                 const Real alpha, const Real Cv0)
      : _rho0(rho0), _e0(e0), _P0(P0), _T0(T0), _A(A), _B(B), _C(C), _G0(G0), _Z(Z),
        _alpha(alpha), _Cv0(Cv0) {}
  DavisReactants GetOnDevice() { return *this; }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    const Real es = Es(rho);
    const Real tmp = std::pow((1.0 + _alpha) / (Ts(rho) * _Cv0) * (sie - es) + 1.0,
                              1.0 / (1.0 + _alpha));
    if (tmp > 0) return Ts(rho) * tmp;
    return Ts(rho) + (sie - es) / _Cv0; // This branch is a negative temperature
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const {
    const Real t_s = Ts(rho);
    return Es(rho) +
           _Cv0 * t_s / (1.0 + _alpha) * (std::pow(temp / t_s, 1.0 + _alpha) - 1.0);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const {
    return PressureFromDensityInternalEnergy(
        rho, InternalEnergyFromDensityTemperature(rho, temp));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return Ps(rho) + Gamma(rho) * rho * (sie - Es(rho));
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  MinInternalEnergyFromDensity(const Real rho, Indexer_t &&lambda = nullptr) const {
    MinInternalEnergyIsNotEnabled("DavisReactants");
    return 0.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    EntropyIsNotEnabled("DavisReactants");
    return 1.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    EntropyIsNotEnabled("DavisReactants");
    return 1.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const {
    return SpecificHeatFromDensityInternalEnergy(
        rho, InternalEnergyFromDensityTemperature(rho, temp));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return _Cv0 / std::pow((1 + _alpha) / (Ts(rho) * _Cv0) * (sie - Es(rho)) + 1,
                           -_alpha / (1 + _alpha));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const {
    return BulkModulusFromDensityInternalEnergy(
        rho, InternalEnergyFromDensityTemperature(rho, temp));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return Gamma(rho);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return Gamma(rho);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                        Real &cv, Real &bmod, const unsigned long output,
                                        Indexer_t &&lambda = nullptr) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                         Real &bmod, Real &dpde, Real &dvdt,
                         Indexer_t &&lambda = nullptr) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const;
  // Generic functions provided by the base class. These contain e.g. the vector
  // overloads that use the scalar versions declared here
  SG_ADD_BASE_CLASS_USINGS(DavisReactants)
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return 0;
  }
  static inline unsigned long max_scratch_size(unsigned int nelements) { return 0; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    static constexpr char s1[]{"DavisReactants Params: "};
    printf("%srho0:%e e0:%e P0:%e\nT0:%e A:%e B:%e\nC:%e G0:%e Z:%e\nalpha:%e "
           "Cv0:%e\n",
           s1, _rho0, _e0, _P0, _T0, _A, _B, _C, _G0, _Z, _alpha, _Cv0);
  }
  void inline Finalize() {}
  static std::string EosType() { return std::string("DavisReactants"); }
  static std::string EosPyType() { return EosType(); }

 private:
  static constexpr Real onethird = 1.0 / 3.0;
  Real _rho0, _e0, _P0, _T0, _A, _B, _C, _G0, _Z, _alpha, _Cv0;
  // static constexpr const char _eos_type[] = "DavisReactants";
  static constexpr unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
  PORTABLE_INLINE_FUNCTION Real Ps(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real Es(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real Ts(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real Gamma(const Real rho) const;
};

class DavisProducts : public EosBase<DavisProducts> {
 public:
  DavisProducts() = default;
  PORTABLE_INLINE_FUNCTION
  DavisProducts(const Real a, const Real b, const Real k, const Real n, const Real vc,
                const Real pc, const Real Cv, const Real E0)
      : _a(a), _b(b), _k(k), _n(n), _vc(vc), _pc(pc), _Cv(Cv), _E0(E0) {}
  DavisProducts GetOnDevice() { return *this; }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return Ts(rho) + (sie - Es(rho)) / _Cv;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const {
    return _Cv * (temp - Ts(rho)) + Es(rho);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const {
    return PressureFromDensityInternalEnergy(
        rho, InternalEnergyFromDensityTemperature(rho, temp));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return Ps(rho) + rho * Gamma(rho) * (sie - Es(rho));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  MinInternalEnergyFromDensity(const Real rho, Indexer_t &&lambda = nullptr) const {
    MinInternalEnergyIsNotEnabled("DavisProducts");
    return 0.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    EntropyIsNotEnabled("DavisProducts");
    return 1.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    EntropyIsNotEnabled("DavisProducts");
    return 1.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return _Cv;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return _Cv;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temp, Indexer_t &&lambda = nullptr) const {
    return BulkModulusFromDensityInternalEnergy(
        rho, InternalEnergyFromDensityTemperature(rho, temp));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    return Gamma(rho);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    return Gamma(rho);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                        Real &cv, Real &bmod, const unsigned long output,
                                        Indexer_t &&lambda = nullptr) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                         Real &bmod, Real &dpde, Real &dvdt,
                         Indexer_t &&lambda = nullptr) const;
  // Generic functions provided by the base class. These contain e.g. the vector
  // overloads that use the scalar versions declared here
  SG_ADD_BASE_CLASS_USINGS(DavisProducts)
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return 0;
  }
  static inline unsigned long max_scratch_size(unsigned int nelements) { return 0; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    static constexpr char s1[]{"DavisProducts Params: "};
    printf("%sa:%e b:%e k:%e\nn:%e vc:%e pc:%e\nCv:%e E0:%e\n", s1, _a, _b, _k, _n, _vc,
           _pc, _Cv, _E0);
  }
  inline void Finalize() {}
  static std::string EosType() { return std::string("DavisProducts"); }
  static std::string EosPyType() { return EosType(); }

 private:
  static constexpr Real onethird = 1.0 / 3.0;
  Real _a, _b, _k, _n, _vc, _pc, _Cv, _E0;
  // static constexpr const char _eos_type[] = "DavisProducts";
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
  PORTABLE_INLINE_FUNCTION Real F(const Real rho) const {
    const Real vvc = 1.0 / (rho * _vc);
    return 2.0 * _a / (std::pow(vvc, 2 * _n) + 1.0);
  }
  PORTABLE_INLINE_FUNCTION Real Ps(const Real rho) const {
    const Real vvc = 1 / (rho * _vc);
    return _pc * std::pow(0.5 * (std::pow(vvc, _n) + std::pow(vvc, -_n)), _a / _n) /
           std::pow(vvc, _k + _a) * (_k - 1.0 + F(rho)) / (_k - 1.0 + _a);
  }
  PORTABLE_INLINE_FUNCTION Real Es(const Real rho) const {
    const Real vvc = 1 / (rho * _vc);
    const Real ec = _pc * _vc / (_k - 1.0 + _a);
    // const Real de = ecj-(Es(rho0)-_E0);
    return ec * std::pow(0.5 * (std::pow(vvc, _n) + std::pow(vvc, -_n)), _a / _n) /
           std::pow(vvc, _k - 1.0 + _a);
  }
  PORTABLE_INLINE_FUNCTION Real Ts(const Real rho) const {
    const Real vvc = 1 / (rho * _vc);
    return std::pow(2.0, -_a * _b / _n) * _pc * _vc / (_Cv * (_k - 1 + _a)) *
           std::pow(0.5 * (std::pow(vvc, _n) + std::pow(vvc, -_n)), _a / _n * (1 - _b)) /
           std::pow(vvc, _k - 1.0 + _a * (1 - _b));
  }
  PORTABLE_INLINE_FUNCTION Real Gamma(const Real rho) const {
    return _k - 1.0 + (1.0 - _b) * F(rho);
  }
};

PORTABLE_INLINE_FUNCTION Real DavisReactants::Ps(const Real rho) const {
  using namespace math_utils;
  const Real y = 1.0 - _rho0 / rho;
  const Real phat = 0.25 * _A * _A / _B * _rho0;
  const Real b4y = 4.0 * _B * y;

  if (rho >= _rho0) {
    return phat *
           (b4y +
            0.5 * (pow<2>(b4y) + onethird * (pow<3>(b4y) + _C * pow<4>(b4y) * 0.25)) +
            pow<2>(y) / pow<4>(1 - y));
  } else {
    return phat * (std::exp(b4y) - 1.0);
  }
}
PORTABLE_INLINE_FUNCTION Real DavisReactants::Es(const Real rho) const {
  const Real y = 1 - _rho0 / rho;
  const Real phat = 0.25 * _A * _A / _B * _rho0;
  const Real b4y = 4 * _B * y;
  Real e_s;
  if (y > 0.0) {
    const Real z = rho / _rho0 - 1;
    e_s = 0.5 * y * b4y *
              (1.0 + onethird * b4y * (1.0 + 0.25 * b4y * (1.0 + _C * 0.2 * b4y))) +
          onethird * math_utils::pow<3>(z);
  } else {
    e_s = -y - (1.0 - std::exp(b4y)) / (4.0 * _B);
  }
  return _e0 + _P0 * (1.0 / _rho0 - 1.0 / rho) + phat / _rho0 * e_s;
}
PORTABLE_INLINE_FUNCTION Real DavisReactants::Ts(const Real rho) const {
  if (rho >= _rho0) {
    const Real y = 1 - _rho0 / rho;
    return _T0 * std::exp(-_Z * y) * std::pow(_rho0 / rho, -_G0 - _Z);
  } else {
    return _T0 * std::pow(_rho0 / rho, -_G0);
  }
}
PORTABLE_INLINE_FUNCTION Real DavisReactants::Gamma(const Real rho) const {
  if (rho >= _rho0) {
    const Real y = 1 - _rho0 / rho;
    return _G0 + _Z * y;
  } else {
    return _G0;
  }
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real DavisReactants::BulkModulusFromDensityInternalEnergy(
    const Real rho, const Real sie, Indexer_t &&lambda) const {
  using namespace math_utils;
  const Real y = 1 - _rho0 / rho;
  const Real phat = 0.25 * _A * _A / _B * _rho0;
  const Real b4y = 4 * _B * y;
  const Real gamma = Gamma(rho);
  const Real esv = -Ps(rho);
  const Real psv =
      (rho >= _rho0)
          ? -phat * _rho0 *
                (4 * _B * (1 + b4y + 0.5 * (pow<2>(b4y) + _C / 3 * pow<3>(b4y))) +
                 3 * y / pow<4>(1 - y) + 4 * pow<2>(y) / pow<5>(1 - y))
          : -phat * 4 * _B * _rho0 * std::exp(b4y);
  const Real gammav = (rho >= _rho0) ? _Z * _rho0 : 0.0;
  return -(psv + (sie - Es(rho)) * rho * (gammav - gamma * rho) - gamma * rho * esv) /
         rho;
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void DavisReactants::DensityEnergyFromPressureTemperature(
    const Real press, const Real temp, Indexer_t &&lambda, Real &rho, Real &sie) const {
  // First, solve P=P(rho,T) for rho.  Note P(rho,e) has an sie-es term, which is only a
  // function of T
  auto PofRatT = [&](const Real r) {
    return (Ps(r) + Gamma(r) * r * _Cv0 * Ts(r) / (1 + _alpha) *
                        (std::pow(temp / Ts(r), 1 + _alpha) - 1.0));
  };
  using RootFinding1D::regula_falsi;
  using RootFinding1D::Status;
  auto status = regula_falsi(PofRatT, press, _rho0, 1.0e-5, 1.0e3, 1.0e-8, 1.0e-8, rho);
  if (status != Status::SUCCESS) {
    // Root finder failed even though the solution was bracketed... this is an error
    EOS_ERROR("DavisReactants::DensityEnergyFromPressureTemperature: "
              "Root find failed to find a solution given P, T\n");
  }
  sie = InternalEnergyFromDensityTemperature(rho, temp);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void DavisReactants::FillEos(Real &rho, Real &temp, Real &sie,
                                                      Real &press, Real &cv, Real &bmod,
                                                      const unsigned long output,
                                                      Indexer_t &&lambda) const {
  if (output & thermalqs::pressure) press = PressureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::temperature)
    temp = TemperatureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::bulk_modulus)
    bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::specific_heat)
    cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
}

// TODO: Chad please decide if this is sane
// TODO(JMM): Pre-cache values instead of computing inline
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void
DavisReactants::ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press,
                                       Real &cv, Real &bmod, Real &dpde, Real &dvdt,
                                       Indexer_t &&lambda) const {
  rho = _rho0;
  temp = _T0;
  sie = InternalEnergyFromDensityTemperature(_rho0, _T0);
  press = PressureFromDensityInternalEnergy(_rho0, sie);
  cv = SpecificHeatFromDensityInternalEnergy(_rho0, sie);
  bmod = BulkModulusFromDensityInternalEnergy(_rho0, sie);
  dpde = Gamma(_rho0) * _rho0;
  // chad please fix this. it's wrong
  Real gm1 = GruneisenParamFromDensityTemperature(_rho0, _T0) * _rho0;
  dvdt = gm1 * cv / bmod;
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real DavisProducts::BulkModulusFromDensityInternalEnergy(
    const Real rho, const Real sie, Indexer_t &&lambda) const {
  using namespace math_utils;
  const Real vvc = 1 / (rho * _vc);
  const Real Fx = -4 * _a * std::pow(vvc, 2 * _n - 1) / pow<2>(1 + std::pow(vvc, 2 * _n));
  const Real tmp = std::pow(0.5 * (std::pow(vvc, _n) + std::pow(vvc, -_n)), _a / _n) /
                   std::pow(vvc, _k + _a);
  const Real tmp_x =
      0.5 * _a * (std::pow(vvc, _n - 1) - std::pow(vvc, -_n - 1)) *
          std::pow(0.5 * (std::pow(vvc, _n) + std::pow(vvc, -_n)), _a / _n - 1) /
          std::pow(vvc, _k + _a) -
      (_k + _a) * tmp / vvc;
  const Real psv = _pc / (_k - 1 + _a) * (tmp * Fx + (_k - 1 + F(rho)) * tmp_x) / _vc;
  // const Real esv = _pc*_vc/(_k-1+_a)*(tmp+vvc*tmp_x)/_vc;
  const Real esv = _pc / (_k - 1 + _a) * (tmp + vvc * tmp_x);
  const Real gamma = Gamma(rho);
  const Real gammav = (1 - _b) * Fx * _vc;
  return -(psv + (sie - Es(rho)) * rho * (gammav - gamma * rho) - gamma * rho * esv) /
         rho;
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void DavisProducts::DensityEnergyFromPressureTemperature(
    const Real press, const Real temp, Indexer_t &&lambda, Real &rho, Real &sie) const {
  auto PofRatT = [&](const Real r) {
    return (Ps(r) + Gamma(r) * r * _Cv * (temp - Ts(r)));
  };
  using RootFinding1D::regula_falsi;
  using RootFinding1D::Status;
  const Real rho0 = 1.0 / _vc;
  auto status = regula_falsi(PofRatT, press, rho0, 1.0e-5, 1.0e3, 1.0e-8, 1.0e-8, rho);
  if (status != Status::SUCCESS) {
    // Root finder failed even though the solution was bracketed... this is an error
    EOS_ERROR("DavisProducts::DensityEnergyFromPressureTemperature: "
              "Root find failed to find a solution given P, T\n");
  }
  sie = InternalEnergyFromDensityTemperature(rho, temp);
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void
DavisProducts::FillEos(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                       Real &bmod, const unsigned long output, Indexer_t &&lambda) const {
  if (output & thermalqs::pressure) press = PressureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::temperature)
    temp = TemperatureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::bulk_modulus)
    bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::specific_heat)
    cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
}
// TODO: pre-cache values instead of computing them
// TODO: chad please decide if these choices are sane
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void
DavisProducts::ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press,
                                      Real &cv, Real &bmod, Real &dpde, Real &dvdt,
                                      Indexer_t &&lambda) const {
  rho = 1.0 / _vc;
  sie = 2. * Es(rho);
  temp = TemperatureFromDensityInternalEnergy(rho, sie);
  press = PressureFromDensityInternalEnergy(rho, sie);
  cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
  bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
  dpde = Gamma(rho) * rho;
  // chad please fix this. it's wrong
  Real gm1 = GruneisenParamFromDensityTemperature(rho, temp) * rho;
  dvdt = gm1 * cv / bmod;
}
} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_EOS_DAVIS_HPP_
