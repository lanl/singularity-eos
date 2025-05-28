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

#ifndef _SINGULARITY_EOS_EOS_EOS_DAVIS_HPP_
#define _SINGULARITY_EOS_EOS_EOS_DAVIS_HPP_

#include <cmath>
#include <cstdio>

#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/math_utils.hpp>
#include <singularity-eos/base/robust_utils.hpp>
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
                 const Real alpha, const Real Cv0,
                 const MeanAtomicProperties &AZbar = MeanAtomicProperties())
      : _rho0(rho0), _spvol0(1.0 / rho0), _e0(e0), _P0(P0), _T0(T0), _A(A), _B(B), _C(C),
        _G0(G0), _Z(Z), _alpha(alpha), _i1pa(1.0 / (1.0 + _alpha)), _Cv0(Cv0),
        _A2oB(_A * _A / _B), _AZbar(AZbar) {
    CheckParams();
  }
  DavisReactants GetOnDevice() { return *this; }
  PORTABLE_INLINE_FUNCTION
  void CheckParams() const {
    PORTABLE_REQUIRE(_rho0 > robust::SMALL(),
                     "Reference density must be strictly positive");
    PORTABLE_REQUIRE(_T0 >= 0, "Reference temperature must be positive");
    PORTABLE_REQUIRE(_B > robust::SMALL(), "B must be strictly positive");
    PORTABLE_REQUIRE(_alpha >= 0, "alpha must be positive");
    _AZbar.CheckParams();
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real power_base = DimlessEdiff(rho, sie);
    if (power_base <= 0) {
      // This case would result in an imaginary temperature (i.e. negative), but we won't
      // allow that so return zero
      return 0.;
    }
    const Real tmp = std::pow(power_base, _i1pa);
    return Ts(rho) * tmp;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temp,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real t_s = Ts(rho);
    PORTABLE_REQUIRE(temp >= 0, "Negative temperature provided");
    return Es(rho) +
           _Cv0 * t_s * _i1pa * (std::pow(robust::ratio(temp, t_s), 1.0 + _alpha) - 1.0);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temp,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PressureFromDensityInternalEnergy(
        rho, InternalEnergyFromDensityTemperature(rho, temp));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return Ps(rho) + Gamma(rho) * rho * (sie - Es(rho));
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MinInternalEnergyFromDensity(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    // Minimum enegy is when the returned temperature is zero. This only happens
    // when the base to the exponent is zero (see T(rho, e) equation)
    const Real es = Es(rho);
    const Real ts = Ts(rho);
    return es - (_Cv0 * ts) * _i1pa;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    EntropyIsNotEnabled("DavisReactants");
    return 1.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    EntropyIsNotEnabled("DavisReactants");
    return 1.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temp,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return SpecificHeatFromDensityInternalEnergy(
        rho, InternalEnergyFromDensityTemperature(rho, temp));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const Real power_base = DimlessEdiff(rho, sie);
    if (power_base <= 0) {
      // Return zero heat capacity instead of an imaginary value
      return 0.;
    }
    return robust::ratio(_Cv0, std::pow(power_base, -_alpha * _i1pa));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temp,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return BulkModulusFromDensityInternalEnergy(
        rho, InternalEnergyFromDensityTemperature(rho, temp));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return Gamma(rho);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return Gamma(rho);
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
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const;
  // Default accessors for mean atomic mass/number
  SG_ADD_DEFAULT_MEAN_ATOMIC_FUNCTIONS(_AZbar)
  // Generic functions provided by the base class. These contain e.g. the vector
  // overloads that use the scalar versions declared here
  SG_ADD_BASE_CLASS_USINGS(DavisReactants)

  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    static constexpr char s1[]{"DavisReactants Params: "};
    printf("%srho0:%e e0:%e P0:%e\nT0:%e A:%e B:%e\nC:%e G0:%e Z:%e\nalpha:%e "
           "Cv0:%e\n",
           s1, _rho0, _e0, _P0, _T0, _A, _B, _C, _G0, _Z, _alpha, _Cv0);
    _AZbar.PrintParams();
  }
  void inline Finalize() {}
  static std::string EosType() { return std::string("DavisReactants"); }
  static std::string EosPyType() { return EosType(); }

 private:
  static constexpr Real onethird = 1.0 / 3.0;
  Real _rho0, _spvol0, _e0, _P0, _T0, _A, _B, _C, _G0, _Z, _alpha, _i1pa, _Cv0, _A2oB;
  MeanAtomicProperties _AZbar;
  // static constexpr const char _eos_type[] = "DavisReactants";
  static constexpr unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
  PORTABLE_FORCEINLINE_FUNCTION Real DimlessEdiff(const Real rho, const Real sie) const;
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
                const Real pc, const Real Cv,
                const MeanAtomicProperties &AZbar = MeanAtomicProperties())
      : _a(a), _b(b), _k(k), _n(n), _aon(_a / _n), _vc(vc), _pc(pc), _Cv(Cv),
        _AZbar(AZbar) {
    CheckParams();
  }
  PORTABLE_INLINE_FUNCTION
  void CheckParams() const {
    PORTABLE_REQUIRE(_n > robust::SMALL(), "n must be strictly positive");
    PORTABLE_REQUIRE(_Cv > robust::SMALL(), "Cv must be strictly positive");
    _AZbar.CheckParams();
  }
  DavisProducts GetOnDevice() { return *this; }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return Ts(rho) + robust::ratio(sie - Es(rho), _Cv);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temp,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return _Cv * (temp - Ts(rho)) + Es(rho);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temp,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PressureFromDensityInternalEnergy(
        rho, InternalEnergyFromDensityTemperature(rho, temp));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return Ps(rho) + rho * Gamma(rho) * (sie - Es(rho));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MinInternalEnergyFromDensity(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    MinInternalEnergyIsNotEnabled("DavisProducts");
    return 0.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    EntropyIsNotEnabled("DavisProducts");
    return 1.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    EntropyIsNotEnabled("DavisProducts");
    return 1.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature,
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
      const Real rho, const Real temp,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return BulkModulusFromDensityInternalEnergy(
        rho, InternalEnergyFromDensityTemperature(rho, temp));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return Gamma(rho);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return Gamma(rho);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const;
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
  // Default accessors for mean atomic mass/number
  SG_ADD_DEFAULT_MEAN_ATOMIC_FUNCTIONS(_AZbar)
  // Generic functions provided by the base class. These contain e.g. the vector
  // overloads that use the scalar versions declared here
  SG_ADD_BASE_CLASS_USINGS(DavisProducts)

  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return 0;
  }
  static inline unsigned long max_scratch_size(unsigned int nelements) { return 0; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    static constexpr char s1[]{"DavisProducts Params: "};
    printf("%sa:%e b:%e k:%e\nn:%e vc:%e pc:%e\nCv:%e \n", s1, _a, _b, _k, _n, _vc, _pc,
           _Cv);
    _AZbar.PrintParams();
  }
  inline void Finalize() {}
  static std::string EosType() { return std::string("DavisProducts"); }
  static std::string EosPyType() { return EosType(); }

 private:
  static constexpr Real onethird = 1.0 / 3.0;
  Real _a, _b, _k, _n, _aon, _vc, _pc, _Cv;
  MeanAtomicProperties _AZbar;
  // static constexpr const char _eos_type[] = "DavisProducts";
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
  PORTABLE_INLINE_FUNCTION Real F(const Real rho) const {
    if (rho <= 0) {
      return 0.;
    }
    const Real vvc = robust::ratio(1.0, rho * _vc);
    return 2.0 * robust::ratio(_a, std::pow(vvc, 2 * _n) + 1.0);
  }
  PORTABLE_INLINE_FUNCTION Real Ps(const Real rho) const {
    if (rho <= 0) {
      return 0.;
    }
    const Real vvc = robust::ratio(1, rho * _vc);
    return _pc *
           robust::ratio(std::pow(0.5 * (std::pow(vvc, _n) + std::pow(vvc, -_n)), _aon),
                         std::pow(vvc, _k + _a)) *
           robust::ratio(_k - 1.0 + F(rho), _k - 1.0 + _a);
  }
  PORTABLE_INLINE_FUNCTION Real Es(const Real rho) const {
    if (rho <= 0) {
      return 0.;
    }
    const Real vvc = robust::ratio(1.0, rho * _vc);
    const Real ec = _pc * robust::ratio(_vc, _k - 1.0 + _a);
    // const Real de = ecj-(Es(rho0)-_E0);
    return robust::ratio(
        ec * std::pow(0.5 * (std::pow(vvc, _n) + std::pow(vvc, -_n)), _aon),
        std::pow(vvc, _k - 1.0 + _a));
  }
  PORTABLE_INLINE_FUNCTION Real Ts(const Real rho) const {
    if (rho <= 0) {
      return 0.;
    }
    const Real vvc = robust::ratio(1.0, rho * _vc);
    return std::pow(2.0, -_a * robust::ratio(_b, _n)) * _pc *
           robust::ratio(_vc, _Cv * (_k - 1 + _a)) *
           robust::ratio(
               std::pow(0.5 * (std::pow(vvc, _n) + std::pow(vvc, -_n)), _aon * (1 - _b)),
               std::pow(vvc, _k - 1.0 + _a * (1 - _b)));
  }
  PORTABLE_INLINE_FUNCTION Real Gamma(const Real rho) const {
    return _k - 1.0 + (1.0 - _b) * F(rho);
  }
};

PORTABLE_FORCEINLINE_FUNCTION Real DavisReactants::DimlessEdiff(const Real rho,
                                                                const Real sie) const {
  return robust::ratio(1.0 + _alpha, Ts(rho) * _Cv0) * (sie - Es(rho)) + 1.0;
}

PORTABLE_INLINE_FUNCTION Real DavisReactants::Ps(const Real rho) const {
  using namespace math_utils;
  const Real y = 1.0 - robust::ratio(_rho0, std::max(rho, 0.));
  const Real phat = 0.25 * _A2oB * _rho0;
  const Real b4y = 4.0 * _B * y;

  if (rho >= _rho0) {
    return phat *
           (b4y +
            0.5 * (pow<2>(b4y) + onethird * (pow<3>(b4y) + _C * pow<4>(b4y) * 0.25)) +
            robust::ratio(pow<2>(y), pow<4>(1 - y)));
  } else {
    return phat * (std::exp(b4y) - 1.0);
  }
}
PORTABLE_INLINE_FUNCTION Real DavisReactants::Es(const Real rho) const {
  const Real y = 1 - robust::ratio(_rho0, std::max(rho, 0.));
  const Real phat = 0.25 * _A2oB * _rho0;
  const Real b4y = 4 * _B * y;
  Real e_s;
  if (y > 0.0) {
    const Real z = rho * _spvol0 - 1;
    e_s = 0.5 * y * b4y *
              (1.0 + onethird * b4y * (1.0 + 0.25 * b4y * (1.0 + _C * 0.2 * b4y))) +
          onethird * math_utils::pow<3>(z);
  } else {
    e_s = -y - robust::ratio(1.0 - std::exp(b4y), 4.0 * _B);
  }
  return _e0 + _P0 * (_spvol0 - robust::ratio(1.0, std::max(rho, 0.))) +
         phat * _spvol0 * e_s;
}
PORTABLE_INLINE_FUNCTION Real DavisReactants::Ts(const Real rho) const {
  const Real rho0overrho = robust::ratio(_rho0, std::max(rho, 0.));
  const Real y = 1 - rho0overrho;
  return _T0 * std::exp(-_Z * y) * std::pow(rho0overrho, -_G0 - _Z);
}
PORTABLE_INLINE_FUNCTION Real DavisReactants::Gamma(const Real rho) const {
  if (rho >= _rho0) {
    const Real y = 1 - robust::ratio(_rho0, std::max(rho, 0.));
    return _G0 + _Z * y;
  } else {
    return _G0;
  }
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real DavisReactants::BulkModulusFromDensityInternalEnergy(
    const Real rho, const Real sie, Indexer_t &&lambda) const {
  using namespace math_utils;
  const Real y = 1 - robust::ratio(_rho0, std::max(rho, 0.));
  const Real phat = 0.25 * _A2oB * _rho0;
  const Real b4y = 4 * _B * y;
  const Real gamma = Gamma(std::max(rho, 0.));
  const Real esv = -Ps(std::max(rho, 0.));
  const Real psv =
      (rho >= _rho0)
          ? -phat * _rho0 *
                (4 * _B * (1 + b4y + 0.5 * (pow<2>(b4y) + onethird * _C * pow<3>(b4y))) +
                 3 * robust::ratio(y, pow<4>(1 - y)) +
                 4 * robust::ratio(pow<2>(y), pow<5>(1 - y)))
          : -phat * 4 * _B * _rho0 * std::exp(b4y);
  const Real gammav = (rho >= _rho0) ? _Z * _rho0 : 0.0;
  const Real numerator =
      -(psv + (sie - Es(rho)) * std::max(rho, 0.) * (gammav - gamma * std::max(rho, 0.)) -
        gamma * std::max(rho, 0.) * esv);
  return robust::ratio(numerator, std::max(rho, 0.));
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void DavisReactants::DensityEnergyFromPressureTemperature(
    const Real press, const Real temp, Indexer_t &&lambda, Real &rho, Real &sie) const {
  // First, solve P=P(rho,T) for rho.  Note P(rho,e) has an sie-es term, which is only a
  // function of T
  PORTABLE_REQUIRE(temp >= 0, "Negative temperature provided");
  auto PofRatT = [&](const Real r) {
    return (Ps(r) + Gamma(r) * r * _Cv0 * Ts(r) * _i1pa *
                        (std::pow(robust::ratio(temp, Ts(r)), 1 + _alpha) - 1.0));
  };
  using RootFinding1D::regula_falsi;
  using RootFinding1D::Status;
  auto status = regula_falsi(PofRatT, press, _rho0, 1.0e-5, 1.0e3, 1.0e-8, 1.0e-8, rho);
  if (status != Status::SUCCESS) {
    // Root finder failed even though the solution was bracketed... this is an error
    EOS_ERROR("DavisReactants::DensityEnergyFromPressureTemperature: "
              "Root find failed to find a solution given P, T\n");
  }
  rho = std::max(0., rho);
  sie = InternalEnergyFromDensityTemperature(rho, temp);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void DavisReactants::FillEos(Real &rho, Real &temp, Real &sie,
                                                      Real &press, Real &cv, Real &bmod,
                                                      const unsigned long output,
                                                      Indexer_t &&lambda) const {
  // BROKEN: This can only handle density-energy inputs! MUST BE CHANGED
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
  dvdt = gm1 * robust::ratio(cv, bmod);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real DavisProducts::BulkModulusFromDensityInternalEnergy(
    const Real rho, const Real sie, Indexer_t &&lambda) const {
  using namespace math_utils;
  if (rho <= 0) {
    return 0.;
  }
  const Real vvc = robust::ratio(1.0, rho * _vc);
  const Real Fx =
      -4 * _a *
      robust::ratio(std::pow(vvc, 2 * _n - 1), pow<2>(1 + std::pow(vvc, 2 * _n)));
  const Real tmp =
      robust::ratio(std::pow(0.5 * (std::pow(vvc, _n) + std::pow(vvc, -_n)), _aon),
                    std::pow(vvc, _k + _a));
  const Real tmp_x =
      0.5 * _a * (std::pow(vvc, _n - 1) - std::pow(vvc, -_n - 1)) *
          robust::ratio(
              std::pow(0.5 * (std::pow(vvc, _n) + std::pow(vvc, -_n)), _aon - 1),
              std::pow(vvc, _k + _a)) -
      (_k + _a) * robust::ratio(tmp, vvc);
  const Real psv = robust::ratio(_pc, _k - 1 + _a) *
                   robust::ratio(tmp * Fx + (_k - 1 + F(rho)) * tmp_x, _vc);
  // const Real esv = _pc*_vc/(_k-1+_a)*(tmp+vvc*tmp_x)/_vc;
  const Real esv = robust::ratio(_pc, _k - 1 + _a) * (tmp + vvc * tmp_x);
  const Real gamma = Gamma(rho);
  const Real gammav = (1 - _b) * Fx * _vc;
  return -robust::ratio(
      psv + (sie - Es(rho)) * rho * (gammav - gamma * rho) - gamma * rho * esv, rho);
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void DavisProducts::DensityEnergyFromPressureTemperature(
    const Real press, const Real temp, Indexer_t &&lambda, Real &rho, Real &sie) const {
  PORTABLE_REQUIRE(temp >= 0, "Negative temperature provided");
  auto PofRatT = [&](const Real r) {
    return (Ps(r) + Gamma(r) * r * _Cv * (temp - Ts(r)));
  };
  using RootFinding1D::regula_falsi;
  using RootFinding1D::Status;
  const Real rho0 = robust::ratio(1.0, _vc);
  auto status = regula_falsi(PofRatT, press, rho0, 1.0e-5, 1.0e3, 1.0e-8, 1.0e-8, rho);
  if (status != Status::SUCCESS) {
    // Root finder failed even though the solution was bracketed... this is an error
    EOS_ERROR("DavisProducts::DensityEnergyFromPressureTemperature: "
              "Root find failed to find a solution given P, T\n");
  }
  rho = std::max(rho, 0.);
  sie = InternalEnergyFromDensityTemperature(rho, temp);
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void
DavisProducts::FillEos(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                       Real &bmod, const unsigned long output, Indexer_t &&lambda) const {
  // BROKEN: This can only handle density-energy inputs! MUST BE CHANGED
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
  rho = robust::ratio(1.0, _vc);
  sie = 2. * Es(rho);
  temp = TemperatureFromDensityInternalEnergy(rho, sie);
  press = PressureFromDensityInternalEnergy(rho, sie);
  cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
  bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
  dpde = Gamma(rho) * rho;
  // chad please fix this. it's wrong
  Real gm1 = GruneisenParamFromDensityTemperature(rho, temp) * rho;
  dvdt = gm1 * robust::ratio(cv, bmod);
}
} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_EOS_DAVIS_HPP_
