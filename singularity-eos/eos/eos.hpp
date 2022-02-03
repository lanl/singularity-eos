//------------------------------------------------------------------------------
// © 2021. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_EOS_EOS_HPP_
#define _SINGULARITY_EOS_EOS_EOS_HPP_

#include "stdio.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <utility>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/eos/eos_variant.hpp>
#include <singularity-eos/eos/eos_base.hpp>

#ifdef SPINER_USE_HDF
#include <fast-math/logs.hpp>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <root-finding-1d/root_finding.hpp>
#include <spiner/databox.hpp>
#include <spiner/spiner_types.hpp>
#endif // SPINER_USE_HDF

#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/eos_error.hpp>
#include <singularity-eos/eos/modifiers/eos_unitsystem.hpp>
#include <singularity-eos/eos/modifiers/relativistic_eos.hpp>
#include <singularity-eos/eos/modifiers/scaled_eos.hpp>
#include <singularity-eos/eos/modifiers/shifted_eos.hpp>

#ifdef SINGULARITY_USE_EOSPAC
#include <eos_Interface.h>
#endif

namespace singularity {

using namespace eos_base;

class IdealGas {
 public:
  IdealGas() = default;
  PORTABLE_INLINE_FUNCTION IdealGas(Real gm1, Real Cv)
      : _Cv(Cv), _gm1(gm1), _rho0(_P0 / (_gm1 * _Cv * _T0)), _sie0(_Cv * _T0),
        _bmod0((_gm1 + 1) * _gm1 * _rho0 * _Cv * _T0), _dpde0(_gm1 * _rho0),
        _dvdt0(1. / (_rho0 * _T0)) {}

  IdealGas GetOnDevice() { return *this; }
  PORTABLE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real PressureFromDensityTemperature(const Real rho,
                                                        const Real temperature,
                                                        Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                                           Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real SpecificHeatFromDensityTemperature(const Real rho,
                                                            const Real temperature,
                                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real BulkModulusFromDensityTemperature(const Real rho,
                                                           const Real temperature,
                                                           Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                 Real &cv, Real &bmod, const unsigned long output,
                                 Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&temperatures,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&sies,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityTemperature(ConstRealIndexer &&rhos,
                                      ConstRealIndexer &&temperatures,
                                      RealIndexer &&pressures,
                                      const int num,
                                      LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&sies,
                                         RealIndexer &&pressures,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                          ConstRealIndexer &&temperatures,
                                          RealIndexer &&cvs,
                                          const int num,
                                          LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                             ConstRealIndexer &&sies,
                                             RealIndexer &&cvs,
                                             const int num,
                                             LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&temperatures,
                                         RealIndexer &&bmods,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&bmods,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&gm1s,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&sies,
                                               RealIndexer &&gm1s,
                                               const int num,
                                               LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename LambdaIndexer>
  inline
  void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
               RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
               const int num, const unsigned long output,
               LambdaIndexer &&lambdas) const;
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  PORTABLE_FUNCTION void PrintParams() const {
    printf("Ideal Gas Parameters:\nGamma = %g\nCv    = %g\n", _gm1 + 1.0, _Cv);
  }
  PORTABLE_FUNCTION void DensityEnergyFromPressureTemperature(const Real press,
                                                              const Real temp,
                                                              Real *lambda, Real &rho,
                                                              Real &sie) const;
  inline void Finalize() {}
  static std::string EosType() { return std::string("IdealGas"); }

 private:
  Real _Cv, _gm1;
  // reference values
  Real _rho0, _sie0, _bmod0, _dpde0, _dvdt0;
  static constexpr const Real _T0 = ROOM_TEMPERATURE;
  static constexpr const Real _P0 = ATMOSPHERIC_PRESSURE;
  // static constexpr const char _eos_type[] = {"IdealGas"};
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
};

// COMMENT: This is meant to be an implementation of the Steinberg version of
// the Gruneisen EOS which should correspond to eostype(3) in xRAGE and
// /[...]/eos/gruneisen in FLAG
class Gruneisen {
 public:
  Gruneisen() = default;
  PORTABLE_INLINE_FUNCTION
  Gruneisen(const Real C0, const Real s1, const Real s2, const Real s3, const Real G0,
            const Real b, const Real rho0, const Real T0, const Real P0, const Real Cv)
      : _C0(C0), _s1(s1), _s2(s2), _s3(s3), _G0(G0), _b(b), _rho0(rho0), _T0(T0), _P0(P0),
        _Cv(Cv) {}
  Gruneisen GetOnDevice() { return *this; }
  PORTABLE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real PressureFromDensityTemperature(const Real rho,
                                                        const Real temperature,
                                                        Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                                           Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real SpecificHeatFromDensityTemperature(const Real rho,
                                                            const Real temperatummmmmmre,
                                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real BulkModulusFromDensityTemperature(const Real rho,
                                                           const Real temperature,
                                                           Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                 Real &cv, Real &bmod, const unsigned long output,
                                 Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&temperatures,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&sies,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityTemperature(ConstRealIndexer &&rhos,
                                      ConstRealIndexer &&temperatures,
                                      RealIndexer &&pressures,
                                      const int num,
                                      LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&sies,
                                         RealIndexer &&pressures,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                          ConstRealIndexer &&temperatures,
                                          RealIndexer &&cvs,
                                          const int num,
                                          LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                             ConstRealIndexer &&sies,
                                             RealIndexer &&cvs,
                                             const int num,
                                             LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&temperatures,
                                         RealIndexer &&bmods,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&bmods,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&gm1s,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&sies,
                                               RealIndexer &&gm1s,
                                               const int num,
                                               LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename LambdaIndexer>
  inline
  void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
               RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
               const int num, const unsigned long output,
               LambdaIndexer &&lambdas) const;
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  PORTABLE_FUNCTION void PrintParams() const {
    static constexpr char s1[]{"Gruneisen Params: "};
    printf("%sC0:%e s1:%e s2:%e\ns3:%e G0:%e b:%e\nrho0:%e T0:%e P0:%e\nCv:%E\n", s1, _C0,
           _s1, _s2, _s3, _G0, _b, _rho0, _T0, _P0, _Cv);
  }
  PORTABLE_FUNCTION void DensityEnergyFromPressureTemperature(const Real press,
                                                              const Real temp,
                                                              Real *lambda, Real &rho,
                                                              Real &sie) const;
  inline void Finalize() {}
  static std::string EosType() { return std::string("Gruneisen"); }

 private:
  Real _C0, _s1, _s2, _s3, _G0, _b, _rho0, _T0, _P0, _Cv;
  // static constexpr const char _eos_type[] = {"Gruneisen"};
  PORTABLE_INLINE_FUNCTION
  Real Gamma(const Real rho) const {
    return rho < _rho0 ? _G0 : _G0 * _rho0 / rho + _b * (1 - _rho0 / rho);
  }
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
};

// COMMENT: This is meant to be an implementation of the "standard" JWL as
// implemented in xRAGE for eostype(1).  It does not include any energy shifting
class JWL {
 public:
  JWL() = default;
  PORTABLE_INLINE_FUNCTION JWL(const Real A, const Real B, const Real R1, const Real R2,
                               const Real w, const Real rho0, const Real Cv)
      : _A(A), _B(B), _R1(R1), _R2(R2), _w(w), _rho0(rho0), _Cv(Cv) {}
  JWL GetOnDevice() { return *this; }
  PORTABLE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real PressureFromDensityTemperature(const Real rho,
                                                        const Real temperature,
                                                        Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                                           Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real SpecificHeatFromDensityTemperature(const Real rho,
                                                            const Real temperature,
                                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real BulkModulusFromDensityTemperature(const Real rho,
                                                           const Real temperature,
                                                           Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                 Real &cv, Real &bmod, const unsigned long output,
                                 Real *lambda = nullptr) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&temperatures,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&sies,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityTemperature(ConstRealIndexer &&rhos,
                                      ConstRealIndexer &&temperatures,
                                      RealIndexer &&pressures,
                                      const int num,
                                      LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&sies,
                                         RealIndexer &&pressures,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                          ConstRealIndexer &&temperatures,
                                          RealIndexer &&cvs,
                                          const int num,
                                          LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                             ConstRealIndexer &&sies,
                                             RealIndexer &&cvs,
                                             const int num,
                                             LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&temperatures,
                                         RealIndexer &&bmods,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&bmods,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&gm1s,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&sies,
                                               RealIndexer &&gm1s,
                                               const int num,
                                               LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename LambdaIndexer>
  inline
  void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
               RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
               const int num, const unsigned long output,
               LambdaIndexer &&lambdas) const;
  PORTABLE_FUNCTION
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  PORTABLE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const;
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  PORTABLE_FUNCTION void PrintParams() const {
    static constexpr char s1[]{"JWL Params: "};
    printf("%sA:%e B:%e R1: %e\nR2:%e w:%e rho0:%e\nCv:%e\n", s1, _A, _B, _R1, _R2, _w,
           _rho0, _Cv);
  }
  PORTABLE_FUNCTION void DensityEnergyFromPressureTemperature(const Real press,
                                                              const Real temp,
                                                              Real *lambda, Real &rho,
                                                              Real &sie) const;
  inline void Finalize() {}
  static std::string EosType() { return std::string("JWL"); }

 private:
  Real _A, _B, _R1, _R2, _w, _rho0, _Cv;
  PORTABLE_INLINE_FUNCTION Real ReferenceEnergy(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real ReferencePressure(const Real rho) const;
  // static constexpr const char _eos_type[] = "JWL";
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
};

class DavisReactants {
 public:
  DavisReactants() = default;
  PORTABLE_INLINE_FUNCTION
  DavisReactants(const Real rho0, const Real e0, const Real P0, const Real T0,
                 const Real A, const Real B, const Real C, const Real G0, const Real Z,
                 const Real alpha, const Real Cv0)
      : _rho0(rho0), _e0(e0), _P0(P0), _T0(T0), _A(A), _B(B), _C(C), _G0(G0), _Z(Z),
        _alpha(alpha), _Cv0(Cv0) {}
  DavisReactants GetOnDevice() { return *this; }
  PORTABLE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real PressureFromDensityTemperature(const Real rho,
                                                        const Real temperature,
                                                        Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                                           Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real SpecificHeatFromDensityTemperature(const Real rho,
                                                            const Real temperature,
                                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real BulkModulusFromDensityTemperature(const Real rho,
                                                           const Real temperature,
                                                           Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                 Real &cv, Real &bmod, const unsigned long output,
                                 Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const;
  PORTABLE_FUNCTION void DensityEnergyFromPressureTemperature(const Real press,
                                                              const Real temp,
                                                              Real *lambda, Real &rho,
                                                              Real &sie) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&temperatures,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&sies,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityTemperature(ConstRealIndexer &&rhos,
                                      ConstRealIndexer &&temperatures,
                                      RealIndexer &&pressures,
                                      const int num,
                                      LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&sies,
                                         RealIndexer &&pressures,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                          ConstRealIndexer &&temperatures,
                                          RealIndexer &&cvs,
                                          const int num,
                                          LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                             ConstRealIndexer &&sies,
                                             RealIndexer &&cvs,
                                             const int num,
                                             LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&temperatures,
                                         RealIndexer &&bmods,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&bmods,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&gm1s,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&sies,
                                               RealIndexer &&gm1s,
                                               const int num,
                                               LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename LambdaIndexer>
  inline
  void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
               RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
               const int num, const unsigned long output,
               LambdaIndexer &&lambdas) const;
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  PORTABLE_FUNCTION void PrintParams() const {
    static constexpr char s1[]{"DavisReactants Params: "};
    printf("%srho0:%e e0:%e P0:%e\nT0:%e A:%e B:%e\nC:%e G0:%e Z:%e\nalpha:%e "
           "Cv0:%e\n",
           s1, _rho0, _e0, _P0, _T0, _A, _B, _C, _G0, _Z, _alpha, _Cv0);
  }
  void inline Finalize() {}
  static std::string EosType() { return std::string("DavisReactants"); }

 private:
  Real _rho0, _e0, _P0, _T0, _A, _B, _C, _G0, _Z, _alpha, _Cv0;
  // static constexpr const char _eos_type[] = "DavisReactants";
  static constexpr unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
  PORTABLE_INLINE_FUNCTION Real Ps(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real Es(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real Ts(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real Gamma(const Real rho) const;
};

class DavisProducts {
 public:
  DavisProducts() = default;
  PORTABLE_INLINE_FUNCTION
  DavisProducts(const Real a, const Real b, const Real k, const Real n, const Real vc,
                const Real pc, const Real Cv, const Real E0)
      : _a(a), _b(b), _k(k), _n(n), _vc(vc), _pc(pc), _Cv(Cv), _E0(E0) {}
  DavisProducts GetOnDevice() { return *this; }
  PORTABLE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real PressureFromDensityTemperature(const Real rho,
                                                        const Real temperature,
                                                        Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                                           Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real SpecificHeatFromDensityTemperature(const Real rho,
                                                            const Real temperature,
                                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real BulkModulusFromDensityTemperature(const Real rho,
                                                           const Real temperature,
                                                           Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION void DensityEnergyFromPressureTemperature(const Real press,
                                                              const Real temp,
                                                              Real *lambda, Real &rho,
                                                              Real &sie) const;
  PORTABLE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                 Real &cv, Real &bmod, const unsigned long output,
                                 Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&temperatures,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&sies,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityTemperature(ConstRealIndexer &&rhos,
                                      ConstRealIndexer &&temperatures,
                                      RealIndexer &&pressures,
                                      const int num,
                                      LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&sies,
                                         RealIndexer &&pressures,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                          ConstRealIndexer &&temperatures,
                                          RealIndexer &&cvs,
                                          const int num,
                                          LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                             ConstRealIndexer &&sies,
                                             RealIndexer &&cvs,
                                             const int num,
                                             LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&temperatures,
                                         RealIndexer &&bmods,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&bmods,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&gm1s,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&sies,
                                               RealIndexer &&gm1s,
                                               const int num,
                                               LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename LambdaIndexer>
  inline
  void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
               RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
               const int num, const unsigned long output,
               LambdaIndexer &&lambdas) const;
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return 0; }
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  PORTABLE_FUNCTION void PrintParams() const {
    static constexpr char s1[]{"DavisProducts Params: "};
    printf("%sa:%e b:%e k:%e\nn:%e vc:%e pc:%e\nCv:%e E0:%e\n", s1, _a, _b, _k, _n, _vc,
           _pc, _Cv, _E0);
  }
  inline void Finalize() {}
  static std::string EosType() { return std::string("DavisProducts"); }

 private:
  Real _a, _b, _k, _n, _vc, _pc, _Cv, _E0;
  // static constexpr const char _eos_type[] = "DavisProducts";
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::specific_internal_energy;
  PORTABLE_INLINE_FUNCTION Real F(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real Ps(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real Es(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real Ts(const Real rho) const;
  PORTABLE_INLINE_FUNCTION Real Gamma(const Real rho) const;
};

#ifdef SPINER_USE_HDF
/*
  Tables all have indep. variables log10(rho), log10(T)

  Extrapolation strategy:
  ------------------------------------------------------------
  We use cold curve for low temperatures/energies
  and assume ideal gas for large temperatures/energies.

  For low densities, we floor the density. For high densities, we
  we use log-linear extrapolation.
*/
class SpinerEOSDependsRhoT {
 public:
  // A weakly typed index map for lambdas
  struct Lambda {
    enum Index { lRho = 0, lT = 1 };
  };

  SpinerEOSDependsRhoT(const std::string &filename, int matid,
                       bool reproduciblity_mode = false);
  SpinerEOSDependsRhoT(const std::string &filename, const std::string &materialName,
                       bool reproducibility_mode = false);
  PORTABLE_INLINE_FUNCTION
  SpinerEOSDependsRhoT() : memoryStatus_(DataStatus::Deallocated) {}

  SpinerEOSDependsRhoT GetOnDevice();

  PORTABLE_FUNCTION
  Real TemperatureFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real InternalEnergyFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real PressureFromDensityTemperature(const Real rho, const Real temperature,
                                      Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                         Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityTemperature(const Real rho, const Real temperature,
                                          Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityInternalEnergy(const Real rho, const Real sie,
                                             Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityTemperature(const Real rho, const Real temperature,
                                         Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityInternalEnergy(const Real rho, const Real sie,
                                               Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  void DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                            Real *lambda, Real &rho, Real &sie) const;
  PORTABLE_FUNCTION
  void PTofRE(const Real rho, const Real sie, Real *lambda, Real &press, Real &temp,
              Real &dpdr, Real &dpde, Real &dtdr, Real &dtde) const;
  PORTABLE_FUNCTION
  void FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
               const unsigned long output, Real *lambda = nullptr) const;

  PORTABLE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&temperatures,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&sies,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityTemperature(ConstRealIndexer &&rhos,
                                      ConstRealIndexer &&temperatures,
                                      RealIndexer &&pressures,
                                      const int num,
                                      LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&sies,
                                         RealIndexer &&pressures,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                          ConstRealIndexer &&temperatures,
                                          RealIndexer &&cvs,
                                          const int num,
                                          LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                             ConstRealIndexer &&sies,
                                             RealIndexer &&cvs,
                                             const int num,
                                             LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&temperatures,
                                         RealIndexer &&bmods,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&bmods,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&gm1s,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&sies,
                                               RealIndexer &&gm1s,
                                               const int num,
                                               LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename LambdaIndexer>
  inline
  void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
               RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
               const int num, const unsigned long output,
               LambdaIndexer &&lambdas) const;
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  std::string filename() const { return std::string(filename_); }
  std::string materialName() const { return std::string(materialName_); }
  int matid() const { return matid_; }
  Real lRhoOffset() const { return lRhoOffset_; }
  Real lTOffset() const { return lTOffset_; }
  Real rhoMin() const { return rho_(lRhoMin_); }
  Real rhoMax() const { return rho_(lRhoMax_); }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    static constexpr char s1[]{"SpinerEOS Parameters:"};
    static constexpr char s2[]{"depends on log_10(rho) and log_10(sie)"};
    static constexpr char s3[]{"EOS file   = "};
    static constexpr char s4[]{"EOS mat ID = "};
    static constexpr char s5[]{"EOS name   = "};
    printf("%s\n\t%s\n\t%s%s\n\t%s%i\n\t%s%s\n", s1, s2, s3, filename_, s4, matid_, s5,
           materialName_);
    return;
  }
  static PORTABLE_FORCEINLINE_FUNCTION int nlambda() { return _n_lambda; }
  inline RootFinding1D::Status rootStatus() const { return status_; }
  inline TableStatus tableStatus() const { return whereAmI_; }
  RootFinding1D::RootCounts counts;
  void Finalize();
  static std::string EosType() { return std::string("SpinerEOSDependsRhoT"); }

 private:
  herr_t loadDataboxes_(const std::string &matid_str, hid_t file, hid_t lTGroup,
                        hid_t coldGroup);
  void fixBulkModulus_();
  void setlTColdCrit_();

  static PORTABLE_FORCEINLINE_FUNCTION Real toLog_(const Real x, const Real offset) {
    // return std::log10(x + offset + EPS);
    // return std::log10(std::abs(std::max(x,-offset) + offset)+EPS);
    return Math::log10(std::abs(std::max(x, -offset) + offset) + EPS);
  }
  static PORTABLE_FORCEINLINE_FUNCTION Real fromLog_(const Real lx, const Real offset) {
    return std::pow(10., lx) - offset;
  }
  PORTABLE_INLINE_FUNCTION Real __attribute__((always_inline))
  lRho_(const Real rho) const noexcept {
    Real out = toLog_(rho, lRhoOffset_);
    return out;
    // return out < lRhoMin_ ? lRhoMin_ : out;
  }
  PORTABLE_INLINE_FUNCTION Real __attribute__((always_inline))
  lT_(const Real T) const noexcept {
    return toLog_(T, lTOffset_);
  }
  PORTABLE_INLINE_FUNCTION Real __attribute__((always_inline))
  rho_(const Real lRho) const noexcept {
    Real rho = fromLog_(lRho, lRhoOffset_);
    return rho < 0 ? 0 : rho;
  }
  PORTABLE_INLINE_FUNCTION Real __attribute__((always_inline))
  T_(const Real lT) const noexcept {
    return fromLog_(lT, lTOffset_);
  }

  PORTABLE_FUNCTION
  Real lTFromlRhoSie_(const Real lRho, const Real sie, TableStatus &whereAmI,
                      Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real lTFromlRhoP_(const Real lRho, const Real press, TableStatus &whereAmI,
                    Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real lRhoFromPlT_(const Real P, const Real lT, TableStatus &whereAmI,
                    Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  void getLogsRhoT_(const Real rho, const Real temperature, Real &lRho, Real &lT,
                    Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real sieFromlRhoTlT_(const Real lRho, const Real T, const Real lT,
                       const TableStatus &whereAmI) const;
  PORTABLE_FUNCTION
  Real PFromRholRhoTlT_(const Real rho, const Real lRho, const Real T, const Real lT,
                        const TableStatus &whereAmI) const;
  PORTABLE_FUNCTION
  Real CvFromlRholT_(const Real lRho, const Real lT, const TableStatus &whereAmI) const;
  PORTABLE_FUNCTION
  Real bModFromRholRhoTlT_(const Real rho, const Real lRho, const Real T, const Real lT,
                           const TableStatus &whereAmI) const;
  PORTABLE_FUNCTION
  TableStatus getLocDependsRhoSie_(const Real lRho, const Real sie) const;
  PORTABLE_FUNCTION
  TableStatus getLocDependsRhoT_(const Real lRho, const Real lT) const;
  // PORTABLE_FUNCTION
  // Real sieToColdInterval_(const Real lRho, const Real sie) const;

  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::temperature;
  // static constexpr const char _eos_type[] {"SpinerEOSDependsRhoT"};
  static constexpr const int numDataBoxes_ = 12;
  Spiner::DataBox P_, sie_, bMod_, dPdRho_, dPdE_, dTdRho_, dTdE_, dEdRho_, dEdT_;
  Spiner::DataBox PMax_, sielTMax_, dEdTMax_, gm1Max_;
  Spiner::DataBox lTColdCrit_;
  Spiner::DataBox PCold_, sieCold_, bModCold_;
  Spiner::DataBox dPdRhoCold_, dPdECold_, dTdRhoCold_, dTdECold_, dEdTCold_;
  int numRho_, numT_;
  Real lRhoMin_, lRhoMax_, rhoMax_;
  Real lRhoMinSearch_;
  Real lTMin_, lTMax_, TMax_;
  Real rhoNormal_, TNormal_, sieNormal_, PNormal_;
  Real CvNormal_, bModNormal_, dPdENormal_, dVdTNormal_;
  Real lRhoOffset_, lTOffset_; // offsets must be non-negative
  const char *filename_;
  const char *materialName_;
  int matid_;
  bool reproducible_;
  // whereAmI_ and status_ used only for reporting. They are not thread-safe.
  mutable TableStatus whereAmI_ = TableStatus::OnTable;
  mutable RootFinding1D::Status status_ = RootFinding1D::Status::SUCCESS;
  static constexpr const Real ROOT_THRESH = 1e-14; // TODO: experiment
  DataStatus memoryStatus_ = DataStatus::Deallocated;
  static constexpr const int _n_lambda = 2;
  static constexpr const char *_lambda_names[2] = {"log(rho)", "log(T)"};
};

/*
  TODO(JMM): Extrapolation Strategy
  ----------------------------------
  Currently the bottom of the table is the bound.
  We do constant extrapolation off the bottom of the table
  and linear extrapolation off the top.

  Since the table is in log-log space, this means extrapolation off
  the top is a power law. Extrapolation off the bottom is constant.
  This worked for nubhlight. But it may or may not work here. We will
  potentially need to revisit this.

  The best solution might be a three-part EOS matched to the bottom of
  the table, containing:
  - A polytropic term
  - A photon pressure term
  - An ideal gas term
  mitigated by Ye and (1-Ye) to control how important each term is.
 */
class SpinerEOSDependsRhoSie {
 public:
  struct SP5Tables {
    Spiner::DataBox P, bMod, dPdRho, dPdE, dTdRho, dTdE, dEdRho;
  };
  PORTABLE_INLINE_FUNCTION SpinerEOSDependsRhoSie()
      : memoryStatus_(DataStatus::Deallocated) {}
  SpinerEOSDependsRhoSie(const std::string &filename, int matid,
                         bool reproducibility_mode = false);
  SpinerEOSDependsRhoSie(const std::string &filename, const std::string &materialName,
                         bool reproducibility_mode = false);
  SpinerEOSDependsRhoSie GetOnDevice();

  PORTABLE_FUNCTION
  Real TemperatureFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real InternalEnergyFromDensityTemperature(const Real rho, const Real T,
                                            Real *lambda = nullptr) const;

  PORTABLE_FUNCTION
  Real PressureFromDensityTemperature(const Real rho, const Real T,
                                      Real *lambda = nullptr) const;

  PORTABLE_FUNCTION
  Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                         Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityTemperature(const Real rho, const Real T,
                                          Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityInternalEnergy(const Real rho, const Real sie,
                                             Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityTemperature(const Real rho, const Real T,
                                         Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityTemperature(const Real rho, const Real T,
                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityInternalEnergy(const Real rho, const Real sie,
                                               Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  void DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                            Real *lambda, Real &rho, Real &sie) const;
  PORTABLE_FUNCTION
  void PTofRE(const Real rho, const Real sie, Real *lambda, Real &press, Real &temp,
              Real &dpdr, Real &dpde, Real &dtdr, Real &dtde) const;
  PORTABLE_FUNCTION
  void FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
               const unsigned long output, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&temperatures,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&sies,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityTemperature(ConstRealIndexer &&rhos,
                                      ConstRealIndexer &&temperatures,
                                      RealIndexer &&pressures,
                                      const int num,
                                      LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&sies,
                                         RealIndexer &&pressures,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                          ConstRealIndexer &&temperatures,
                                          RealIndexer &&cvs,
                                          const int num,
                                          LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                             ConstRealIndexer &&sies,
                                             RealIndexer &&cvs,
                                             const int num,
                                             LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&temperatures,
                                         RealIndexer &&bmods,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&bmods,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&gm1s,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&sies,
                                               RealIndexer &&gm1s,
                                               const int num,
                                               LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename LambdaIndexer>
  inline
  void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
               RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
               const int num, const unsigned long output,
               LambdaIndexer &&lambdas) const;
  PORTABLE_FUNCTION
  unsigned long PreferredInput() const { return _preferred_input; }
  std::string filename() const { return std::string(filename_); }
  std::string materialName() const { return std::string(materialName_); }
  int matid() const { return matid_; }
  Real lRhoOffset() const { return lRhoOffset_; }
  Real lTOffset() const { return lTOffset_; }
  Real lEOffset() const { return lEOffset_; }
  Real rhoMin() const { return fromLog_(lRhoMin_, lRhoOffset_); }
  Real rhoMax() const { return fromLog_(lRhoMax_, lRhoOffset_); }
  Real TMin() const { return fromLog_(sie_.range(0).min(), lTOffset_); }
  Real TMax() const { return fromLog_(sie_.range(0).max(), lTOffset_); }
  Real sieMin() const { return fromLog_(T_.range(0).min(), lEOffset_); }
  Real sieMax() const { return fromLog_(T_.range(0).max(), lEOffset_); }

  static PORTABLE_FORCEINLINE_FUNCTION int nlambda() { return _n_lambda; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    static constexpr char s1[]{"SpinerEOS Parameters:"};
    static constexpr char s2[]{"depends on log_10(rho) and log_10(sie)"};
    static constexpr char s3[]{"EOS file   = "};
    static constexpr char s4[]{"EOS mat ID = "};
    static constexpr char s5[]{"EOS name   = "};
    printf("%s\n\t%s\n\t%s%s\n\t%s%i\n\t%s%s\n", s1, s2, s3, filename_, s4, matid_, s5,
           materialName_);
    return;
  }
  RootFinding1D::Status rootStatus() const { return status_; }
  RootFinding1D::RootCounts counts;
  static std::string EosType() { return std::string("SpinerEOSDependsRhoSie"); }
  void Finalize();

 private:
  herr_t loadDataboxes_(const std::string &matid_str, hid_t file, hid_t lTGroup,
                        hid_t lEGroup);
  void calcBMod_(SP5Tables &tables);

  static PORTABLE_FORCEINLINE_FUNCTION Real toLog_(const Real x, const Real offset) {
    // return std::log10(std::abs(std::max(x,-offset) + offset)+EPS);
    return Math::log10(std::abs(std::max(x, -offset) + offset) + EPS);
  }
  static PORTABLE_FORCEINLINE_FUNCTION Real fromLog_(const Real lx, const Real offset) {
    return std::pow(10., lx) - offset;
  }
  PORTABLE_FUNCTION
  Real interpRhoT_(const Real rho, const Real T, const Spiner::DataBox &db,
                   Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real interpRhoSie_(const Real rho, const Real sie, const Spiner::DataBox &db,
                     Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real lRhoFromPlT_(const Real P, const Real lT, Real *lambda) const;

  Spiner::DataBox sie_; // depends on (rho,T)
  Spiner::DataBox T_;   // depends on (rho, sie)
  SP5Tables dependsRhoT_;
  SP5Tables dependsRhoSie_;
  int numRho_;
  Real rhoNormal_, TNormal_, sieNormal_, PNormal_;
  Real CvNormal_, bModNormal_, dPdENormal_, dVdTNormal_;
  Real lRhoMin_, lRhoMax_, rhoMax_;
  Spiner::DataBox PlRhoMax_, dPdRhoMax_;

  Real lRhoOffset_, lTOffset_, lEOffset_; // offsets must be non-negative

  static constexpr unsigned long _preferred_input =
      thermalqs::density | thermalqs::temperature;
  // static constexpr const char _eos_type[] = "SpinerEOSDependsRhoSie";
  const char *filename_;
  const char *materialName_;
  int matid_;
  bool reproducible_;
  mutable RootFinding1D::Status status_;
  static constexpr const int _n_lambda = 1;
  static constexpr const char *_lambda_names[1] = {"log(rho)"};
  DataStatus memoryStatus_ = DataStatus::Deallocated;
};

// Note the Stellar Collapse tables have units of:
// 1. Ye (unitless)
// 2. log(MeV) for temperature
// 3. log(g/cm^3) for density.
//
// TODO(JMM): For now the bottom of the table is a floor and the top
// is linear extrapolation in log-log space. We should reconsider this
// and introduce extrapolation as needed.
class StellarCollapse {
 public:
  // A weakly typed index map for lambdas
  struct Lambda {
    enum Index { Ye = 0, lT = 1 };
  };

  StellarCollapse(const std::string &filename, bool use_sp5 = false,
                  bool filter_bmod = true);

  // Saves to an SP5 file
  void Save(const std::string &filename);

  PORTABLE_INLINE_FUNCTION
  StellarCollapse() : memoryStatus_(DataStatus::Deallocated) {}

  StellarCollapse GetOnDevice();

  PORTABLE_FUNCTION
  Real TemperatureFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real InternalEnergyFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real PressureFromDensityTemperature(const Real rho, const Real temperature,
                                      Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                         Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityTemperature(const Real rho, const Real temperature,
                                          Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real SpecificHeatFromDensityInternalEnergy(const Real rho, const Real sie,
                                             Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityTemperature(const Real rho, const Real temperature,
                                         Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real BulkModulusFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  Real GruneisenParamFromDensityInternalEnergy(const Real rho, const Real sie,
                                               Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  void DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                            Real *lambda, Real &rho, Real &sie) const;
  PORTABLE_FUNCTION
  void FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
               const unsigned long output, Real *lambda = nullptr) const;

  PORTABLE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&temperatures,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&sies,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityTemperature(ConstRealIndexer &&rhos,
                                      ConstRealIndexer &&temperatures,
                                      RealIndexer &&pressures,
                                      const int num,
                                      LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&sies,
                                         RealIndexer &&pressures,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                          ConstRealIndexer &&temperatures,
                                          RealIndexer &&cvs,
                                          const int num,
                                          LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                             ConstRealIndexer &&sies,
                                             RealIndexer &&cvs,
                                             const int num,
                                             LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&temperatures,
                                         RealIndexer &&bmods,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&bmods,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&gm1s,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&sies,
                                               RealIndexer &&gm1s,
                                               const int num,
                                               LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename LambdaIndexer>
  inline
  void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
               RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
               const int num, const unsigned long output,
               LambdaIndexer &&lambdas) const;
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  std::string filename() const { return std::string(filename_); }
  Real lRhoOffset() const { return lRhoOffset_; }
  Real lTOffset() const { return lTOffset_; }
  Real lEOffset() const { return lEOffset_; }
  Real lRhoMin() const { return lRhoMin_; }
  Real lRhoMax() const { return lRhoMax_; }
  Real rhoMin() const { return rho_(lRhoMin_); }
  Real rhoMax() const { return rho_(lRhoMax_); }
  Real lTMin() const { return lTMin_; }
  Real lTMax() const { return lTMax_; }
  Real TMin() const { return T_(lTMin_); }
  Real TMax() const { return T_(lTMax_); }
  Real YeMin() const { return YeMin_; }
  Real YeMax() const { return YeMax_; }
  Real sieMin() const { return sieMin_; }
  Real sieMax() const { return sieMax_; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    printf("StellarCollapse parameters:\n"
           "depends on log10(rho), log10(T), Ye\n"
           "EOS file = %s\n"
           "lrho bounds = %g, %g\n"
           "lT bounds = %g, %g\n"
           "Ye bounds = %g, %g\n",
           filename_, lRhoMin_, lRhoMax_, lTMin_, lTMax_, YeMin_, YeMax_);
    return;
  }
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return _n_lambda; }
  inline RootFinding1D::Status rootStatus() const { return status_; }
  RootFinding1D::RootCounts counts;
  void Finalize();
  static std::string EosType() { return std::string("StellarCollapse"); }

 private:
  void LoadFromSP5File_(const std::string &filename);
  void LoadFromStellarCollapseFile_(const std::string &filename);
  int readSCInt_(const hid_t &file_id, const std::string &name);
  void readBounds_(const hid_t &file_id, const std::string &name, int size, Real &lo,
                   Real &hi);
  void readSCDset_(const hid_t &file_id, const std::string &name, Spiner::DataBox &db);

  void medianFilter_(Spiner::DataBox &db);
  void medianFilter_(const Spiner::DataBox &in, Spiner::DataBox &out);
  void fillMedianBuffer_(Real buffer[], int width, int iY, int iT, int irho,
                         const Spiner::DataBox &tab) const;
  Real findMedian_(Real buffer[], int size) const;
  void computeBulkModulus_();
  void computeColdAndHotCurves_();
  void setNormalValues_();

  PORTABLE_INLINE_FUNCTION __attribute__((always_inline)) void
  checkLambda_(Real *lambda) const noexcept {
    if (lambda == nullptr) {
      EOS_ERROR("StellarCollapse: lambda must contain Ye and 1 space for caching.\n");
    }
  }

  PORTABLE_INLINE_FUNCTION __attribute__((always_inline)) Real
  toLog_(const Real x, const Real offset) const noexcept {
    // return std::log10(x + offset + EPS);
    // return std::log10(std::abs(std::max(x,-offset) + offset)+EPS);
    return Math::log10(std::abs(std::max(x, -offset) + offset) + EPS);
  }
  PORTABLE_INLINE_FUNCTION Real __attribute__((always_inline))
  fromLog_(const Real lx, const Real offset) const noexcept {
    return std::pow(10., lx) - offset;
  }
  PORTABLE_INLINE_FUNCTION Real __attribute__((always_inline))
  lRho_(const Real rho) const noexcept {
    Real out = toLog_(rho, lRhoOffset_);
    return out;
    // return out < lRhoMin_ ? lRhoMin_ : out;
  }
  // TODO(JMM): Use log K2MeV so we can do addition instead of
  // multiplication?
  PORTABLE_INLINE_FUNCTION Real __attribute__((always_inline))
  lT_(const Real T) const noexcept {
    return toLog_(T * K2MeV_, lTOffset_);
  }
  PORTABLE_INLINE_FUNCTION Real __attribute__((always_inline))
  rho_(const Real lRho) const noexcept {
    Real rho = fromLog_(lRho, lRhoOffset_);
    return rho < 0 ? 0 : rho;
  }
  PORTABLE_INLINE_FUNCTION Real __attribute__((always_inline))
  T_(const Real lT) const noexcept {
    return MeV2K_ * fromLog_(lT, lTOffset_);
  }
  PORTABLE_INLINE_FUNCTION Real __attribute__((always_inline))
  le2e_(const Real le) const noexcept {
    return fromLog_(le, lEOffset_);
  }
  PORTABLE_INLINE_FUNCTION Real __attribute__((always_inline))
  e2le_(const Real e) const noexcept {
    return toLog_(e, lEOffset_);
  }
  PORTABLE_INLINE_FUNCTION Real __attribute__((always_inline))
  lP2P_(const Real lP) const noexcept {
    return fromLog_(lP, lPOffset_);
  }
  PORTABLE_INLINE_FUNCTION Real __attribute__((always_inline))
  P2lP_(const Real P) const noexcept {
    return toLog_(P, lPOffset_);
  }
  PORTABLE_INLINE_FUNCTION Real __attribute__((always_inline))
  lB2B_(const Real lB) const noexcept {
    return fromLog_(lB, lBOffset_);
  }
  PORTABLE_INLINE_FUNCTION Real __attribute__((always_inline))
  B2lB_(const Real B) const noexcept {
    return toLog_(B, lBOffset_);
  }

  PORTABLE_FUNCTION Real lTFromlRhoSie_(const Real lRho, const Real sie,
                                        Real *lambda) const noexcept;
  PORTABLE_INLINE_FUNCTION __attribute__((always_inline)) void
  getLogsFromRhoT_(const Real rho, const Real temp, Real *lambda, Real &lRho, Real &lT,
                   Real &Ye) const noexcept {
    checkLambda_(lambda);
    lRho = lRho_(rho);
    lT = lT_(temp);
    Ye = lambda[Lambda::Ye];
    lambda[Lambda::lT] = lT;
  }
  PORTABLE_INLINE_FUNCTION __attribute__((always_inline)) void
  getLogsFromRhoSie_(const Real rho, const Real sie, Real *lambda, Real &lRho, Real &lT,
                     Real &Ye) const noexcept {
    lRho = lRho_(rho);
    lT = lTFromlRhoSie_(lRho, sie, lambda);
    Ye = lambda[Lambda::Ye];
    return;
  }

  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::temperature;

  // Dependent variables
  Spiner::DataBox lP_, lE_, dPdRho_, dPdE_, dEdT_, lBMod_;
  Spiner::DataBox entropy_; // kb/baryon
  Spiner::DataBox Xa_;      // mass fraction of alpha particles
  Spiner::DataBox Xh_;      // mass fraction of heavy ions
  Spiner::DataBox Xn_;      // mass fraction of neutrons
  Spiner::DataBox Xp_;      // mass fraction of protons
  Spiner::DataBox Abar_;    // Average atomic mass
  Spiner::DataBox Zbar_;    // Average atomic number
  // Spiner::DataBox gamma_; // polytropic index. dlog(P)/dlog(rho).
  // dTdRho_, dTdE_, dEdRho_, dEdT_;

  // Bounds of dependent variables. Needed for root finding.
  Spiner::DataBox eCold_, eHot_;

  // Independent variable bounds
  int numRho_, numT_, numYe_;
  Real lRhoMin_, lRhoMax_;
  Real lTMin_, lTMax_;
  Real YeMin_, YeMax_;
  Real sieMin_, sieMax_;

  static constexpr Real MeV2GK_ = 11.604525006;
  static constexpr Real GK2MeV_ = 1. / MeV2GK_;
  static constexpr Real MeV2K_ = 1.e6 * 1.60217653e-12;
  static constexpr Real K2MeV_ = 1. / MeV2K_;
  static constexpr Real TNormal_ = 5 * GK2MeV_; // Threshold of NSE
  static constexpr Real rhoNormal_ = 2.e12;     // 1./100'th of nuclear density
  static constexpr Real YeNormal_ = 0.3517;     // Beta equilibrium value
  Real sieNormal_, PNormal_;
  Real CvNormal_, bModNormal_, dPdENormal_, dVdTNormal_;

  // offsets must be non-negative
  Real lEOffset_;
  static constexpr Real lRhoOffset_ = 0.0; // TODO(JMM): Address if this ever changes
  static constexpr Real lTOffset_ = 0.0;
  static constexpr Real lPOffset_ = 0.0;
  static constexpr Real lBOffset_ = 0.0;

  const char *filename_;
  // whereAmI_ and status_ used only for reporting. They are not thread-safe.
  mutable RootFinding1D::Status status_ = RootFinding1D::Status::SUCCESS;
  static constexpr const Real ROOT_THRESH = 1e-14; // TODO: experiment
  DataStatus memoryStatus_ = DataStatus::Deallocated;
  static constexpr const int _n_lambda = 2;
  static constexpr const char *_lambda_names[] = {"Ye", "log(T)"};

  // Stuff for median filter smoothing
  static constexpr Real EPSSMOOTH = 10.0;
  static constexpr int MF_W = 3;
  static constexpr int MF_S = (2 * MF_W + 1) * (2 * MF_W + 1) * (2 * MF_W + 1);
};
#endif // SPINER_USE_HDF

#ifdef SINGULARITY_USE_EOSPAC
// Only really works in serial
// Not really supported on device
class EOSPAC {
 public:
  EOSPAC() = default;
  EOSPAC(int matid, bool invert_at_setup = false);
  EOSPAC GetOnDevice() { return *this; }
  PORTABLE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real PressureFromDensityTemperature(const Real rho,
                                                        const Real temperature,
                                                        Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                                           Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real SpecificHeatFromDensityTemperature(const Real rho,
                                                            const Real temperature,
                                                            Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real BulkModulusFromDensityTemperature(const Real rho,
                                                           const Real temperature,
                                                           Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                 Real &cv, Real &bmod, const unsigned long output,
                                 Real *lambda = nullptr) const;
  PORTABLE_FUNCTION
  void DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                            Real *lambda, Real &rho, Real &sie) const;
  PORTABLE_FUNCTION void PTofRE(const Real rho, const Real sie, Real *lambda, Real &press,
                                Real &temp, Real &dpdr, Real &dpde, Real &dtdr,
                                Real &dtde) const;
  PORTABLE_FUNCTION void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie,
                                                Real &press, Real &cv, Real &bmod,
                                                Real &dpde, Real &dvdt,
                                                Real *lambda = nullptr) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&temperatures,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&sies,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityTemperature(ConstRealIndexer &&rhos,
                                      ConstRealIndexer &&temperatures,
                                      RealIndexer &&pressures,
                                      const int num,
                                      LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&sies,
                                         RealIndexer &&pressures,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                          ConstRealIndexer &&temperatures,
                                          RealIndexer &&cvs,
                                          const int num,
                                          LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                             ConstRealIndexer &&sies,
                                             RealIndexer &&cvs,
                                             const int num,
                                             LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                         ConstRealIndexer &&temperatures,
                                         RealIndexer &&bmods,
                                         const int num,
                                         LambdaIndexer &&lambdas) const;
  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&sies,
                                            RealIndexer &&bmods,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&gm1s,
                                            const int num,
                                            LambdaIndexer &&lambdas) const;

  template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline
  void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&sies,
                                               RealIndexer &&gm1s,
                                               const int num,
                                               LambdaIndexer &&lambdas) const;

  template<typename RealIndexer, typename LambdaIndexer>
  inline
  void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
               RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
               const int num, const unsigned long output,
               LambdaIndexer &&lambdas) const;

  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  int nlambda() const noexcept { return 0; }
  inline void Finalize() {}
  static std::string EosType() { return std::string("EOSPAC"); }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    printf("EOSPAC parameters:\nmatid = %s\n", matid_);
  }

 private:
  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::temperature;
  int matid_;
  static constexpr int NT = 5;
  EOS_INTEGER PofRT_table_;
  EOS_INTEGER TofRE_table_;
  EOS_INTEGER EofRT_table_;
  EOS_INTEGER RofPT_table_;
  EOS_INTEGER TofRP_table_;
  EOS_INTEGER tablehandle[NT];
  EOS_INTEGER EOS_Info_table_;
  static constexpr Real temp_ref_ = 293;
  Real rho_ref_ = 1;
  Real sie_ref_ = 1;
  Real press_ref_ = 1;
  Real cv_ref_ = 1;
  Real bmod_ref_ = 1;
  Real dpde_ref_ = 1;
  Real dvdt_ref_ = 1;
};
#endif // SINGULARITY_USE_EOSPAC

using EOS = Variant<
    IdealGas, Gruneisen, JWL, DavisReactants, DavisProducts, ScaledEOS<IdealGas>,
    ShiftedEOS<IdealGas>, ShiftedEOS<ScaledEOS<IdealGas>>,
    ScaledEOS<ShiftedEOS<IdealGas>>, RelativisticEOS<IdealGas>, ScaledEOS<Gruneisen>,
    ShiftedEOS<Gruneisen>, ScaledEOS<ShiftedEOS<Gruneisen>>, ScaledEOS<JWL>,
    ShiftedEOS<JWL>, ScaledEOS<ShiftedEOS<JWL>>, ScaledEOS<DavisReactants>,
    ShiftedEOS<DavisReactants>, ScaledEOS<ShiftedEOS<DavisReactants>>,
    ScaledEOS<DavisProducts>, ShiftedEOS<DavisProducts>,
    ScaledEOS<ShiftedEOS<DavisProducts>>, UnitSystem<IdealGas>
#ifdef SPINER_USE_HDF
    ,
    SpinerEOSDependsRhoT, SpinerEOSDependsRhoSie, ScaledEOS<SpinerEOSDependsRhoT>,
    ScaledEOS<SpinerEOSDependsRhoSie>, ShiftedEOS<SpinerEOSDependsRhoT>,
    ShiftedEOS<SpinerEOSDependsRhoSie>, ShiftedEOS<ScaledEOS<SpinerEOSDependsRhoT>>,
    ShiftedEOS<ScaledEOS<SpinerEOSDependsRhoSie>>,
    ScaledEOS<ShiftedEOS<SpinerEOSDependsRhoT>>,
    ScaledEOS<ShiftedEOS<SpinerEOSDependsRhoSie>>, RelativisticEOS<SpinerEOSDependsRhoT>,
    RelativisticEOS<SpinerEOSDependsRhoSie>, UnitSystem<SpinerEOSDependsRhoT>,
    UnitSystem<SpinerEOSDependsRhoSie>,
    // TODO(JMM): Might need shifted + relativistic
    // for StellarCollapse. Might not. Negative
    // energies can throw off normalization of cs2 by
    // enthalpy weirdly.
    StellarCollapse, ScaledEOS<StellarCollapse>, ShiftedEOS<StellarCollapse>,
    ShiftedEOS<ScaledEOS<StellarCollapse>>, RelativisticEOS<StellarCollapse>,
    UnitSystem<StellarCollapse>
#endif // SPINER_USE_HDF
#ifdef SINGULARITY_USE_EOSPAC
    ,
    EOSPAC, ScaledEOS<EOSPAC>, ShiftedEOS<EOSPAC>, ShiftedEOS<ScaledEOS<EOSPAC>>,
    ScaledEOS<ShiftedEOS<EOSPAC>>
#endif // SINGULARITY_USE_EOSPAC
    >;

} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_EOS_HPP_
