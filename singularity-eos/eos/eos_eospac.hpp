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

#ifndef _SINGULARITY_EOS_EOS_EOS_EOSPAC_HPP_
#define _SINGULARITY_EOS_EOS_EOS_EOSPAC_HPP_
#ifdef SINGULARITY_USE_EOSPAC

#include <algorithm>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include <eos_Interface.h>
// ports-of-call
#include <ports-of-call/portability.hpp>

#include <eospac-wrapper/eospac_wrapper.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {

using namespace eos_base;

// How do I make a destructor? How do I free the EOS memory if more than one
// points to the same table? Does EOSPAC give me multiple (reference-counted)
// handles to the same table or read in the table multiple times?

// Only really works in serial
// Not really supported on device
#ifndef PORTABILITY_STRATEGY_NONE
#if defined(__CUDACC__)
#define SG_PIF_NOWARN #pragma nv_exec_check_disable
#endif // __CUDACC__
#else
#define SG_PIF_NOWARN
#endif // PORTABILITY_STRATEGY_NONE

class EOSPAC : public EosBase<EOSPAC> {
 public:
  inline EOSPAC() = default;
  inline EOSPAC(int matid, bool invert_at_setup = false);
  inline EOSPAC GetOnDevice() { return *this; }
  SG_PIF_NOWARN
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  SG_PIF_NOWARN
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  SG_PIF_NOWARN
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  SG_PIF_NOWARN
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  SG_PIF_NOWARN
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  SG_PIF_NOWARN
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  SG_PIF_NOWARN
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  SG_PIF_NOWARN
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  SG_PIF_NOWARN
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  SG_PIF_NOWARN
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  SG_PIF_NOWARN
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  SG_PIF_NOWARN
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  SG_PIF_NOWARN
  PORTABLE_INLINE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                        Real &cv, Real &bmod, const unsigned long output,
                                        Real *lambda = nullptr) const;
  SG_PIF_NOWARN
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp, Real *lambda,
                                       Real &rho, Real &sie) const;
  SG_PIF_NOWARN
  PORTABLE_INLINE_FUNCTION void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie,
                                                       Real &press, Real &cv, Real &bmod,
                                                       Real &dpde, Real &dvdt,
                                                       Real *lambda = nullptr) const;

  // Generic (Scalar)
  using EosBase<EOSPAC>::is_raw_pointer;
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                       RealIndexer &&temperatures, const int num,
                                       LambdaIndexer &&lambdas) const {
    EosBase<EOSPAC>::TemperatureFromDensityInternalEnergy(rhos, sies, temperatures, num,
                                                          lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void
  TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                       RealIndexer &&temperatures, Real * /*scratch*/,
                                       const int num, LambdaIndexer &&lambdas) const {
    PORTABLE_WARN("EOSPAC type mismatch will cause significant performance degradation");
    EosBase<EOSPAC>::TemperatureFromDensityInternalEnergy(rhos, sies, temperatures, num,
                                                          lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&sies, const int num,
                                                   LambdaIndexer &&lambdas) const {
    EosBase<EOSPAC>::InternalEnergyFromDensityTemperature(rhos, temperatures, sies, num,
                                                          lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&sies, Real * /*scratch*/,
                                                   const int num,
                                                   LambdaIndexer &&lambdas) const {
    PORTABLE_WARN("EOSPAC type mismatch will cause significant performance degradation");
    EosBase<EOSPAC>::InternalEnergyFromDensityTemperature(rhos, temperatures, sies, num,
                                                          lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void PressureFromDensityTemperature(ConstRealIndexer &&rhos,
                                             ConstRealIndexer &&temperatures,
                                             RealIndexer &&pressures, const int num,
                                             LambdaIndexer &&lambdas) const {
    EosBase<EOSPAC>::PressureFromDensityTemperature(rhos, temperatures, pressures, num,
                                                    lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void
  PressureFromDensityTemperature(ConstRealIndexer &&rhos, ConstRealIndexer &&temperatures,
                                 RealIndexer &&pressures, Real * /*scratch*/,
                                 const int num, LambdaIndexer &&lambdas) const {
    PORTABLE_WARN("EOSPAC type mismatch will cause significant performance degradation");
    EosBase<EOSPAC>::PressureFromDensityTemperature(rhos, temperatures, pressures, num,
                                                    lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                ConstRealIndexer &&sies,
                                                RealIndexer &&pressures, const int num,
                                                LambdaIndexer &&lambdas) const {
    EosBase<EOSPAC>::PressureFromDensityInternalEnergy(rhos, sies, pressures, num,
                                                       lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void
  PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                    RealIndexer &&pressures, Real * /*scratch*/,
                                    const int num, LambdaIndexer &&lambdas) const {
    PORTABLE_WARN("EOSPAC type mismatch will cause significant performance degradation");
    EosBase<EOSPAC>::PressureFromDensityInternalEnergy(rhos, sies, pressures, num,
                                                       lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void EntropyFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&entropies, const int num,
                                            LambdaIndexer &&lambdas) const {
    EosBase<EOSPAC>::EntropyFromDensityTemperature(rhos, temperatures, entropies, num,
                                                   lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void
  EntropyFromDensityTemperature(ConstRealIndexer &&rhos, ConstRealIndexer &&temperatures,
                                RealIndexer &&entropies, Real * /*scratch*/,
                                const int num, LambdaIndexer &&lambdas) const {
    PORTABLE_WARN("EOSPAC type mismatch will cause significant performance degradation");
    EosBase<EOSPAC>::EntropyFromDensityTemperature(rhos, temperatures, entropies, num,
                                                   lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void EntropyFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&sies,
                                               RealIndexer &&entropies, const int num,
                                               LambdaIndexer &&lambdas) const {
    EosBase<EOSPAC>::EntropyFromDensityInternalEnergy(rhos, sies, entropies, num,
                                                      lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void
  EntropyFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                   RealIndexer &&entropies, Real * /*scratch*/,
                                   const int num, LambdaIndexer &&lambdas) const {
    PORTABLE_WARN("EOSPAC type mismatch will cause significant performance degradation");
    EosBase<EOSPAC>::EntropyFromDensityInternalEnergy(rhos, sies, entropies, num,
                                                      lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                                 ConstRealIndexer &&temperatures,
                                                 RealIndexer &&cvs, const int num,
                                                 LambdaIndexer &&lambdas) const {
    EosBase<EOSPAC>::SpecificHeatFromDensityTemperature(rhos, temperatures, cvs, num,
                                                        lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                                 ConstRealIndexer &&temperatures,
                                                 RealIndexer &&cvs, Real * /*scratch*/,
                                                 const int num,
                                                 LambdaIndexer &&lambdas) const {
    PORTABLE_WARN("EOSPAC type mismatch will cause significant performance degradation");
    EosBase<EOSPAC>::SpecificHeatFromDensityTemperature(rhos, temperatures, cvs, num,
                                                        lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                    ConstRealIndexer &&sies,
                                                    RealIndexer &&cvs, const int num,
                                                    LambdaIndexer &&lambdas) const {
    EosBase<EOSPAC>::SpecificHeatFromDensityInternalEnergy(rhos, sies, cvs, num, lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void
  SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                        RealIndexer &&cvs, Real * /*scratch*/,
                                        const int num, LambdaIndexer &&lambdas) const {
    PORTABLE_WARN("EOSPAC type mismatch will cause significant performance degradation");
    EosBase<EOSPAC>::SpecificHeatFromDensityInternalEnergy(rhos, sies, cvs, num, lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                                ConstRealIndexer &&temperatures,
                                                RealIndexer &&bmods, const int num,
                                                LambdaIndexer &&lambdas) const {
    EosBase<EOSPAC>::BulkModulusFromDensityTemperature(rhos, temperatures, bmods, num,
                                                       lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                                ConstRealIndexer &&temperatures,
                                                RealIndexer &&bmods, Real * /*scratch*/,
                                                const int num,
                                                LambdaIndexer &&lambdas) const {
    PORTABLE_WARN("EOSPAC type mismatch will cause significant performance degradation");
    EosBase<EOSPAC>::BulkModulusFromDensityTemperature(rhos, temperatures, bmods, num,
                                                       lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&sies,
                                                   RealIndexer &&bmods, const int num,
                                                   LambdaIndexer &&lambdas) const {
    EosBase<EOSPAC>::BulkModulusFromDensityInternalEnergy(rhos, sies, bmods, num,
                                                          lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void
  BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                       RealIndexer &&bmods, Real * /*scratch*/,
                                       const int num, LambdaIndexer &&lambdas) const {
    PORTABLE_WARN("EOSPAC type mismatch will cause significant performance degradation");
    EosBase<EOSPAC>::BulkModulusFromDensityInternalEnergy(rhos, sies, bmods, num,
                                                          lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&gm1s, const int num,
                                                   LambdaIndexer &&lambdas) const {
    EosBase<EOSPAC>::GruneisenParamFromDensityTemperature(rhos, temperatures, gm1s, num,
                                                          lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&gm1s, Real * /*scratch*/,
                                                   const int num,
                                                   LambdaIndexer &&lambdas) const {
    PORTABLE_WARN("EOSPAC type mismatch will cause significant performance degradation");
    EosBase<EOSPAC>::GruneisenParamFromDensityTemperature(rhos, temperatures, gm1s, num,
                                                          lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                      ConstRealIndexer &&sies,
                                                      RealIndexer &&gm1s, const int num,
                                                      LambdaIndexer &&lambdas) const {
    EosBase<EOSPAC>::GruneisenParamFromDensityInternalEnergy(rhos, sies, gm1s, num,
                                                             lambdas);
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                      ConstRealIndexer &&sies,
                                                      RealIndexer &&gm1s,
                                                      Real * /*scratch*/, const int num,
                                                      LambdaIndexer &&lambdas) const {
    PORTABLE_WARN("EOSPAC type mismatch will cause significant performance degradation");
    EosBase<EOSPAC>::GruneisenParamFromDensityInternalEnergy(rhos, sies, gm1s, num,
                                                             lambdas);
  }

  using EosBase<EOSPAC>::PTofRE;
  using EosBase<EOSPAC>::FillEos;
  using EosBase<EOSPAC>::EntropyIsNotEnabled;

  // EOSPAC vector implementations
  template <typename LambdaIndexer>
  inline void
  TemperatureFromDensityInternalEnergy(const Real *rhos, const Real *sies,
                                       Real *temperatures, Real *scratch, const int num,
                                       LambdaIndexer /*lambdas*/,
                                       Transform &&transform = Transform()) const {
    using namespace EospacWrapper;
    EOS_REAL *R = const_cast<EOS_REAL *>(&rhos[0]);
    EOS_REAL *E = const_cast<EOS_REAL *>(&sies[0]);
    EOS_REAL *T = &temperatures[0];
    EOS_REAL *dTdr = scratch + 0 * num;
    EOS_REAL *dTde = scratch + 1 * num;

    EOS_INTEGER table = TofRE_table_;
    EOS_INTEGER options[3];
    EOS_REAL values[3];
    EOS_INTEGER nopts = 0;

    if (transform.x.is_set()) {
      options[nopts] = EOS_X_CONVERT;
      values[nopts] = 1.0 / transform.x.get();
      ++nopts;
    }

    options[nopts] = EOS_Y_CONVERT;
    values[nopts] = sieFromSesame(1.0);
    if (transform.y.is_set()) {
      values[nopts] /= transform.y.get();
    }
    ++nopts;

    if (transform.f.is_set()) {
      options[nopts] = EOS_F_CONVERT;
      values[nopts] = transform.f.get();
      ++nopts;
    }

    eosSafeInterpolate(&table, num, R, E, T, dTdr, dTde, "TofRE", Verbosity::Quiet,
                       options, values, nopts);
  }

  template <typename LambdaIndexer>
  inline void PressureFromDensityTemperature(const Real *rhos, const Real *temperatures,
                                             Real *pressures, Real *scratch,
                                             const int num, LambdaIndexer /*lambdas*/,
                                             Transform &&transform = Transform()) const {
    using namespace EospacWrapper;
    EOS_REAL *R = const_cast<EOS_REAL *>(&rhos[0]);
    EOS_REAL *T = const_cast<EOS_REAL *>(&temperatures[0]);
    EOS_REAL *P = &pressures[0];
    EOS_REAL *dPdr = scratch + 0 * num;
    EOS_REAL *dPdT = scratch + 1 * num;

    EOS_INTEGER table = PofRT_table_;
    EOS_INTEGER options[3];
    EOS_REAL values[3];
    EOS_INTEGER nopts = 0;

    if (!transform.x.is_set() && !transform.y.is_set()) {
      options[nopts] = EOS_XY_PASSTHRU;
      values[nopts] = 1.0;
      ++nopts;
    } else {
      if (transform.x.is_set()) {
        options[nopts] = EOS_X_CONVERT;
        values[nopts] = 1.0 / transform.x.get();
        ++nopts;
      }

      if (transform.y.is_set()) {
        options[nopts] = EOS_Y_CONVERT;
        values[nopts] = 1.0 / transform.y.get();
        ++nopts;
      }
    }

    options[nopts] = EOS_F_CONVERT;
    values[nopts] = pressureFromSesame(1.0);

    if (transform.f.is_set()) {
      values[nopts] *= transform.f.get();
    }
    ++nopts;

    eosSafeInterpolate(&table, num, R, T, P, dPdr, dPdT, "PofRT", Verbosity::Quiet,
                       options, values, nopts);
  }

  template <typename LambdaIndexer>
  inline void FillEos(Real *rhos, Real *temps, Real *energies, Real *presses, Real *cvs,
                      Real *bmods, Real *scratch, const int num,
                      const unsigned long output, LambdaIndexer /*lambdas*/) const {
    static auto const name =
        singularity::mfuncname::member_func_name(typeid(EOSPAC).name(), __func__);
    static auto const cname = name.c_str();
    using namespace EospacWrapper;
    EOS_REAL *R = &rhos[0];
    EOS_REAL *T = &temps[0];
    EOS_REAL *E = &energies[0];
    EOS_REAL *P = &presses[0];
    EOS_REAL *DPDT = scratch + 0 * num;
    EOS_REAL *DPDR = scratch + 1 * num;
    EOS_REAL *DEDT = scratch + 2 * num;
    EOS_REAL *DEDR = scratch + 3 * num;
    EOS_REAL *dx = DPDT;
    EOS_REAL *dy = DPDR;

    EOS_INTEGER nxypairs = num;

    const unsigned long input = ~output;
    if (output == thermalqs::none) {
      UNDEFINED_ERROR;
    }

    if (output & thermalqs::density) {
      if (input & thermalqs::pressure && input & thermalqs::temperature) {
        EOS_INTEGER table = RofPT_table_;
        EOS_INTEGER options[]{EOS_X_CONVERT};
        EOS_REAL values[]{pressureFromSesame(1.0)};
        EOS_INTEGER nopts = 1;
        eosSafeInterpolate(&table, nxypairs, P, T, R, dx, dy, "RofPT", Verbosity::Quiet,
                           options, values, nopts);
      } else {
        UNDEFINED_ERROR;
      }
    }
    if (output & thermalqs::temperature) {
      if (input & thermalqs::density && input & thermalqs::specific_internal_energy) {
        EOS_INTEGER table = TofRE_table_;
        EOS_INTEGER options[]{EOS_Y_CONVERT};
        EOS_REAL values[]{sieFromSesame(1.0)};
        EOS_INTEGER nopts = 1;
        eosSafeInterpolate(&table, nxypairs, R, E, T, dx, dy, "TofRE", Verbosity::Quiet,
                           options, values, nopts);
      } else if (input & thermalqs::density && input & thermalqs::pressure) {
        EOS_INTEGER table = TofRP_table_;
        EOS_INTEGER options[]{EOS_Y_CONVERT};
        EOS_REAL values[]{pressureFromSesame(1.0)};
        EOS_INTEGER nopts = 1;
        eosSafeInterpolate(&table, nxypairs, R, P, T, dx, dy, "TofRP", Verbosity::Quiet,
                           options, values, nopts);
      } else {
        UNDEFINED_ERROR;
      }
    }
    if ((output & thermalqs::specific_internal_energy) ||
        (output & thermalqs::specific_heat || output & thermalqs::bulk_modulus)) {
      EOS_INTEGER table = EofRT_table_;
      EOS_INTEGER options[]{EOS_F_CONVERT};
      EOS_REAL values[]{sieFromSesame(1.0)};
      EOS_INTEGER nopts = 1;
      eosSafeInterpolate(&table, nxypairs, R, T, E, DEDR, DEDT, "EofRT", Verbosity::Quiet,
                         options, values, nopts);
    }
    if ((output & thermalqs::pressure) || (output & thermalqs::bulk_modulus)) {
      EOS_INTEGER table = PofRT_table_;
      EOS_INTEGER options[]{EOS_F_CONVERT};
      EOS_REAL values[]{pressureFromSesame(1.0)};
      EOS_INTEGER nopts = 1;
      eosSafeInterpolate(&table, nxypairs, R, T, P, DPDR, DPDT, "PofRT", Verbosity::Quiet,
                         options, values, nopts);
      // Thermodynamics: Bt = rho*dP/drho|T, Bs=Bt*Cv/Cp, Cv=dE/dT|v, and
      // Cp-Cv=-T(dP/dT|v)^2/(dP/dV|T) or Cp-Cv=-T/rho^2 (dP/dT|v)^2/(dP/drho|T)
      // Therefore: Bs=Bt Cv/(CV+above)
      // BMOD = rho*dx[0]*CV/(CV-T[0]/(rho*rho)*dy[0]*dy[0]/dx[0]);
    }
    if (output & thermalqs::specific_heat) {
      portableFor(
          cname, 0, num, PORTABLE_LAMBDA(const int i) {
            cvs[i] = std::max(DEDT[i], 0.0); // Here we do something to the data!
          });
    }
    if (output & thermalqs::bulk_modulus) {
      portableFor(
          cname, 0, num, PORTABLE_LAMBDA(const int i) {
            const Real rho = R[i];
            Real BMOD = 0.0;
            if (DEDT[i] > 0.0 && rho > 0.0) {
              const Real DPDE = DPDT[i] / DEDT[i];
              BMOD = rho * DPDR[i] + DPDE * (P[i] / rho - rho * DEDR[i]);
            } else if (rho > 0.0) { // Case: DEDT <= 0
              // We need a different DPDE call apparently????
              // In xRAGE, they call out to P(rho,e) in this case to get the
              // derivative directly from the half-inverted table.
              // But upon further review,
              // I think it will end up evaluating to BMOD_T in any case because cv will
              // be zero! See eos_eospac.f90 line 1261 BMOD =
              // BMOD_T+DPDE*DPDE*std::max(DEDT,0.0)*T[0]/rho;
              // BMOD_T = std::max(rho * DPDR[i], 0.0);
              // BMOD = BMOD_T; //+DPDE*DPDE*std::max(DEDT,0.0)*T[0]/rho;
              BMOD = std::max(rho * DPDR[i], 0.0);
            }
            bmods[i] = std::max(BMOD, 0.0);
          });
    }
  }

  template <typename LambdaIndexer>
  inline void
  InternalEnergyFromDensityTemperature(const Real *rhos, const Real *temperatures,
                                       Real *sies, Real *scratch, const int num,
                                       LambdaIndexer /*lambdas*/,
                                       Transform &&transform = Transform()) const {
    static auto const name =
        singularity::mfuncname::member_func_name(typeid(EOSPAC).name(), __func__);
    static auto const cname = name.c_str();
    using namespace EospacWrapper;
    EOS_REAL *R = const_cast<EOS_REAL *>(&rhos[0]);
    EOS_REAL *T = const_cast<EOS_REAL *>(&temperatures[0]);
    EOS_REAL *E = &sies[0];
    EOS_REAL *DEDT = scratch + 0 * num;
    EOS_REAL *DEDR = scratch + 1 * num;

    EOS_INTEGER table = EofRT_table_;
    EOS_INTEGER options[3];
    EOS_REAL values[3];
    EOS_INTEGER nopts = 0;

    if (!transform.x.is_set() && !transform.y.is_set()) {
      options[nopts] = EOS_XY_PASSTHRU;
      values[nopts] = 1.0;
      ++nopts;
    } else {
      if (transform.x.is_set()) {
        options[nopts] = EOS_X_CONVERT;
        values[nopts] = 1.0 / transform.x.get();
        ++nopts;
      }

      if (transform.y.is_set()) {
        options[nopts] = EOS_Y_CONVERT;
        values[nopts] = 1.0 / transform.y.get();
        ++nopts;
      }
    }

    options[nopts] = EOS_F_CONVERT;
    values[nopts] = sieFromSesame(1.0);

    if (transform.f.is_set()) {
      values[nopts] *= transform.f.get();
    }
    ++nopts;

    eosSafeInterpolate(&table, num, R, T, E, DEDR, DEDT, "EofRT", Verbosity::Quiet,
                       options, values, nopts);
  }

  template <typename LambdaIndexer>
  inline void PressureFromDensityInternalEnergy(
      const Real *rhos, const Real *sies, Real *pressures, Real *scratch, const int num,
      LambdaIndexer /*lambdas*/, Transform &&transform = Transform()) const {
    using namespace EospacWrapper;
    EOS_REAL *R = const_cast<EOS_REAL *>(&rhos[0]);
    EOS_REAL *E = const_cast<EOS_REAL *>(&sies[0]);
    EOS_REAL *P = &pressures[0];
    EOS_REAL *T = scratch + 0 * num;
    EOS_REAL *dTdr = scratch + 1 * num;
    EOS_REAL *dTde = scratch + 2 * num;
    EOS_REAL *dPdr = dTdr;
    EOS_REAL *dPdT = dTde;

    EOS_INTEGER table = TofRE_table_;
    {
      EOS_INTEGER options[2];
      EOS_REAL values[2];
      EOS_INTEGER nopts = 0;

      if (transform.x.is_set()) {
        options[nopts] = EOS_X_CONVERT;
        values[nopts] = 1.0 / transform.x.get();
        ++nopts;
      }

      options[nopts] = EOS_Y_CONVERT;
      values[nopts] = sieFromSesame(1.0);

      if (transform.y.is_set()) {
        values[nopts] /= transform.y.get();
      }
      ++nopts;

      eosSafeInterpolate(&table, num, R, E, T, dTdr, dTde, "TofRE", Verbosity::Quiet,
                         options, values, nopts);
    }

    table = PofRT_table_;
    {
      EOS_INTEGER options[2];
      EOS_REAL values[2];
      EOS_INTEGER nopts = 0;

      if (transform.x.is_set()) {
        options[nopts] = EOS_X_CONVERT;
        values[nopts] = 1.0 / transform.x.get();
      } else {
        options[nopts] = EOS_XY_PASSTHRU;
        values[nopts] = 1.0;
      }
      ++nopts;

      options[nopts] = EOS_F_CONVERT;
      values[nopts] = pressureFromSesame(1.0);

      if (transform.f.is_set()) {
        values[nopts] *= transform.f.get();
      }
      ++nopts;

      eosSafeInterpolate(&table, num, R, T, P, dPdr, dPdT, "PofRT", Verbosity::Quiet,
                         options, values, nopts);
    }
  }

  template <typename LambdaIndexer>
  inline void SpecificHeatFromDensityTemperature(
      const Real *rhos, const Real *temperatures, Real *cvs, Real *scratch, const int num,
      LambdaIndexer /*lambdas*/, Transform &&transform = Transform()) const {
    static auto const name =
        singularity::mfuncname::member_func_name(typeid(EOSPAC).name(), __func__);
    static auto const cname = name.c_str();
    using namespace EospacWrapper;
    EOS_REAL *R = const_cast<EOS_REAL *>(&rhos[0]);
    EOS_REAL *T = const_cast<EOS_REAL *>(&temperatures[0]);
    EOS_REAL *E = scratch + 0 * num;
    EOS_REAL *DEDR = scratch + 1 * num;
    EOS_REAL *DEDT = &cvs[0];

    EOS_INTEGER table = EofRT_table_;

    EOS_INTEGER options[3];
    EOS_REAL values[3];
    EOS_INTEGER nopts = 0;

    if (!transform.x.is_set() && !transform.y.is_set()) {
      options[nopts] = EOS_XY_PASSTHRU;
      values[nopts] = 1.0;
      ++nopts;
    } else {
      if (transform.x.is_set()) {
        options[nopts] = EOS_X_CONVERT;
        values[nopts] = 1.0 / transform.x.get();
        ++nopts;
      }
      if (transform.y.is_set()) {
        options[nopts] = EOS_Y_CONVERT;
        values[nopts] = 1.0 / transform.y.get();
        ++nopts;
      }
    }

    eosSafeInterpolate(&table, num, R, T, E, DEDR, DEDT, "EofRT", Verbosity::Quiet,
                       options, values, nopts);

    const Real y = transform.y.is_set() ? (1.0 / transform.y.get()) : 1.0;
    const Real f = transform.f.is_set() ? transform.f.get() : 1.0;

    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          cvs[i] =
              f * y *
              cvFromSesame(std::max(DEDT[i], 0.0)); // Here we do something to the data!
        });
  }

  template <typename LambdaIndexer>
  inline void SpecificHeatFromDensityInternalEnergy(
      const Real *rhos, const Real *sies, Real *cvs, Real *scratch, const int num,
      LambdaIndexer /*lambdas*/, Transform &&transform = Transform()) const {

    static auto const name =
        singularity::mfuncname::member_func_name(typeid(EOSPAC).name(), __func__);
    static auto const cname = name.c_str();
    using namespace EospacWrapper;
    EOS_REAL *R = const_cast<EOS_REAL *>(&rhos[0]);
    EOS_REAL *E = const_cast<EOS_REAL *>(&sies[0]);
    EOS_REAL *T = scratch + 0 * num;
    EOS_REAL *dTdr = scratch + 1 * num;
    EOS_REAL *dTde = scratch + 2 * num;
    EOS_REAL *NE = scratch + 3 * num;
    EOS_REAL *DEDT = dTdr;
    EOS_REAL *DEDR = dTde;

    EOS_INTEGER table = TofRE_table_;
    {
      EOS_INTEGER options[2];
      EOS_REAL values[2];
      EOS_INTEGER nopts = 0;

      if (transform.x.is_set()) {
        options[nopts] = EOS_X_CONVERT;
        values[nopts] = 1.0 / transform.x.get();
        ++nopts;
      }

      options[nopts] = EOS_Y_CONVERT;
      values[nopts] = sieFromSesame(1.0);
      if (transform.y.is_set()) {
        values[nopts] /= transform.y.get();
      }
      ++nopts;

      eosSafeInterpolate(&table, num, R, E, T, dTdr, dTde, "TofRE", Verbosity::Quiet,
                         options, values, nopts);
    }

    table = EofRT_table_;
    {
      EOS_INTEGER options[2];
      EOS_REAL values[2];
      EOS_INTEGER nopts = 0;

      if (transform.x.is_set()) {
        options[nopts] = EOS_X_CONVERT;
        values[nopts] = 1.0 / transform.x.get();
      } else {
        options[nopts] = EOS_XY_PASSTHRU;
        values[nopts] = 1.0;
      }
      ++nopts;

      eosSafeInterpolate(&table, num, R, T, NE, DEDR, DEDT, "EofRT", Verbosity::Quiet,
                         options, values, nopts);
    }

    const Real f = transform.f.is_set() ? transform.f.get() : 1.0;

    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          cvs[i] = f * cvFromSesame(
                           std::max(DEDT[i], 0.0)); // Here we do something to the data!
        });
  }

  template <typename LambdaIndexer>
  inline void
  BulkModulusFromDensityTemperature(const Real *rhos, const Real *temperatures,
                                    Real *bmods, Real *scratch, const int num,
                                    LambdaIndexer /*lambdas*/,
                                    Transform &&transform = Transform()) const {
    static auto const name =
        singularity::mfuncname::member_func_name(typeid(EOSPAC).name(), __func__);
    static auto const cname = name.c_str();
    using namespace EospacWrapper;
    EOS_REAL *R = const_cast<EOS_REAL *>(&rhos[0]);
    EOS_REAL *T = const_cast<EOS_REAL *>(&temperatures[0]);
    EOS_REAL *E = scratch + 0 * num;
    EOS_REAL *P = scratch + 1 * num;
    EOS_REAL *DPDT = scratch + 2 * num;
    EOS_REAL *DPDR = scratch + 3 * num;
    EOS_REAL *DEDT = scratch + 4 * num;
    EOS_REAL *DEDR = scratch + 5 * num;

    EOS_INTEGER options[2];
    EOS_REAL values[2];
    EOS_INTEGER nopts = 0;

    if (!transform.x.is_set() && !transform.y.is_set()) {
      options[nopts] = EOS_XY_PASSTHRU;
      values[nopts] = 1.0;
      ++nopts;
    } else {
      if (transform.x.is_set()) {
        options[nopts] = EOS_X_CONVERT;
        values[nopts] = 1.0 / transform.x.get();
        ++nopts;
      }
      if (transform.y.is_set()) {
        options[nopts] = EOS_Y_CONVERT;
        values[nopts] = 1.0 / transform.y.get();
        ++nopts;
      }
    }

    EOS_INTEGER table = EofRT_table_;
    eosSafeInterpolate(&table, num, R, T, E, DEDR, DEDT, "EofRT", Verbosity::Quiet,
                       options, values, nopts);

    table = PofRT_table_;
    eosSafeInterpolate(&table, num, R, T, P, DPDR, DPDT, "PofRT", Verbosity::Quiet,
                       options, values, nopts);

    // WARNING: DEDR and DPDR are divided by EOS_X_CONVERT
    // This is why the BMOD calculation was changed below compared to the scalar version.

    // Thermodynamics: Bt = rho*dP/drho|T, Bs=Bt*Cv/Cp, Cv=dE/dT|v, and
    // Cp-Cv=-T(dP/dT|v)^2/(dP/dV|T) or Cp-Cv=-T/rho^2 (dP/dT|v)^2/(dP/drho|T)
    // Therefore: Bs=Bt Cv/(CV+above)
    // BMOD = rho*dx[0]*CV/(CV-T[0]/(rho*rho)*dy[0]*dy[0]/dx[0]);
    const Real x = transform.x.is_set() ? transform.x.get() : 1.0;
    const Real f = transform.f.is_set() ? transform.f.get() : 1.0;

    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          const Real rho = R[i];
          Real BMOD = 0.0;
          if (DEDT[i] > 0.0 && rho > 0.0) {
            const Real DPDE = DPDT[i] / DEDT[i];
            BMOD = rho * DPDR[i] + DPDE * (P[i] / (rho * x) - rho * DEDR[i]);
          } else if (rho > 0.0) { // Case: DEDT <= 0
            // We need a different DPDE call apparently????
            // In xRAGE, they call out to P(rho,e) in this case to get the
            // derivative directly from the half-inverted table.
            // But upon further review,
            // I think it will end up evaluating to BMOD_T in any case because cv will
            // be zero! See eos_eospac.f90 line 1261 BMOD =
            // BMOD_T+DPDE*DPDE*std::max(DEDT,0.0)*T[0]/rho;
            // BMOD_T = std::max(rho * DPDR[i], 0.0);
            // BMOD = BMOD_T; //+DPDE*DPDE*std::max(DEDT,0.0)*T[0]/rho;
            BMOD = std::max(rho * DPDR[i], 0.0);
          }
          bmods[i] = f * bulkModulusFromSesame(std::max(BMOD, 0.0));
        });
  }

  template <typename LambdaIndexer>
  inline void BulkModulusFromDensityInternalEnergy(
      const Real *rhos, const Real *sies, Real *bmods, Real *scratch, const int num,
      LambdaIndexer /*lambdas*/, Transform &&transform = Transform()) const {
    static auto const name =
        singularity::mfuncname::member_func_name(typeid(EOSPAC).name(), __func__);
    static auto const cname = name.c_str();
    using namespace EospacWrapper;
    EOS_REAL *R = const_cast<EOS_REAL *>(&rhos[0]);
    EOS_REAL *E = const_cast<EOS_REAL *>(&sies[0]);
    EOS_REAL *T = scratch + 0 * num;
    EOS_REAL *dTdr = scratch + 1 * num;
    EOS_REAL *dTde = scratch + 2 * num;
    EOS_REAL *Etmp = scratch + 3 * num;
    EOS_REAL *DEDT = scratch + 4 * num;
    EOS_REAL *DEDR = scratch + 5 * num;
    EOS_REAL *P = Etmp;
    EOS_REAL *DPDT = dTdr;
    EOS_REAL *DPDR = dTde;

    EOS_INTEGER table = TofRE_table_;
    {
      EOS_INTEGER options[2];
      EOS_REAL values[2];
      EOS_INTEGER nopts = 0;
      if (transform.x.is_set()) {
        options[nopts] = EOS_X_CONVERT;
        values[nopts] = 1.0 / transform.x.get();
        ++nopts;
      }
      options[nopts] = EOS_Y_CONVERT;
      values[nopts] = sieFromSesame(1.0);
      if (transform.y.is_set()) {
        values[nopts] /= transform.y.get();
      }
      ++nopts;
      eosSafeInterpolate(&table, num, R, E, T, dTdr, dTde, "TofRE", Verbosity::Quiet,
                         options, values, nopts);
    }

    EOS_INTEGER options[1];
    EOS_REAL values[1];
    EOS_INTEGER nopts = 0;
    if (transform.x.is_set()) {
      options[nopts] = EOS_X_CONVERT;
      values[nopts] = 1.0 / transform.x.get();
    } else {
      options[nopts] = EOS_XY_PASSTHRU;
      values[nopts] = 1.0;
    }
    ++nopts;

    table = EofRT_table_;
    eosSafeInterpolate(&table, num, R, T, Etmp, DEDR, DEDT, "EofRT", Verbosity::Quiet,
                       options, values, nopts);

    table = PofRT_table_;
    eosSafeInterpolate(&table, num, R, T, P, DPDR, DPDT, "PofRT", Verbosity::Quiet,
                       options, values, nopts);

    // WARNING: DEDR and DPDR are divided by EOS_X_CONVERT
    // This is why the BMOD calculation was changed below compared to the scalar version.

    // Thermodynamics: Bt = rho*dP/drho|T, Bs=Bt*Cv/Cp, Cv=dE/dT|v, and
    // Cp-Cv=-T(dP/dT|v)^2/(dP/dV|T) or Cp-Cv=-T/rho^2 (dP/dT|v)^2/(dP/drho|T)
    // Therefore: Bs=Bt Cv/(CV+above)
    // BMOD = rho*dx[0]*CV/(CV-T[0]/(rho*rho)*dy[0]*dy[0]/dx[0]);
    const Real x = transform.x.is_set() ? transform.x.get() : 1.0;
    const Real f = transform.f.is_set() ? transform.f.get() : 1.0;

    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          const Real rho = R[i];
          Real BMOD = 0.0;
          if (DEDT[i] > 0.0 && rho > 0.0) {
            const Real DPDE = DPDT[i] / DEDT[i];
            BMOD = rho * DPDR[i] + DPDE * (P[i] / (rho * x) - rho * DEDR[i]);
          } else if (rho > 0.0) { // Case: DEDT <= 0
            // We need a different DPDE call apparently????
            // In xRAGE, they call out to P(rho,e) in this case to get the
            // derivative directly from the half-inverted table.
            // But upon further review,
            // I think it will end up evaluating to BMOD_T in any case because cv will
            // be zero! See eos_eospac.f90 line 1261 BMOD =
            // BMOD_T+DPDE*DPDE*std::max(DEDT,0.0)*T[0]/rho;
            // BMOD_T = std::max(rho * DPDR[i], 0.0);
            // BMOD = BMOD_T; //+DPDE*DPDE*std::max(DEDT,0.0)*T[0]/rho;
            BMOD = std::max(rho * DPDR[i], 0.0);
          }
          bmods[i] = f * bulkModulusFromSesame(std::max(BMOD, 0.0));
        });
  }

  template <typename LambdaIndexer>
  inline void
  GruneisenParamFromDensityTemperature(const Real *rhos, const Real *temperatures,
                                       Real *gm1s, Real *scratch, const int num,
                                       LambdaIndexer /*lambdas*/,
                                       Transform &&transform = Transform()) const {
    static auto const name =
        singularity::mfuncname::member_func_name(typeid(EOSPAC).name(), __func__);
    static auto const cname = name.c_str();
    using namespace EospacWrapper;
    EOS_REAL *R = const_cast<EOS_REAL *>(&rhos[0]);
    EOS_REAL *T = const_cast<EOS_REAL *>(&temperatures[0]);
    EOS_REAL *E = scratch + 0 * num;
    EOS_REAL *dx = scratch + 1 * num;
    EOS_REAL *DEDT = scratch + 2 * num;
    EOS_REAL *DPDT = scratch + 3 * num;
    EOS_REAL *P = E;

    EOS_INTEGER options[2];
    EOS_REAL values[2];
    EOS_INTEGER nopts = 0;

    if (!transform.x.is_set() && !transform.y.is_set()) {
      options[nopts] = EOS_XY_PASSTHRU;
      values[nopts] = 1.0;
      ++nopts;
    } else {
      if (transform.x.is_set()) {
        options[nopts] = EOS_X_CONVERT;
        values[nopts] = 1.0 / transform.x.get();
        ++nopts;
      }
      if (transform.y.is_set()) {
        options[nopts] = EOS_Y_CONVERT;
        values[nopts] = 1.0 / transform.y.get();
        ++nopts;
      }
    }
    EOS_INTEGER table = EofRT_table_;
    eosSafeInterpolate(&table, num, R, T, E, dx, DEDT, "EofRT", Verbosity::Quiet, options,
                       values, nopts);

    table = PofRT_table_;
    eosSafeInterpolate(&table, num, R, T, P, dx, DPDT, "PofRT", Verbosity::Quiet, options,
                       values, nopts);

    const Real x = transform.x.is_set() ? transform.x.get() : 1.0;
    const Real f = transform.f.is_set() ? transform.f.get() : 1.0;

    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          const Real DPDE = DPDT[i] / DEDT[i];
          gm1s[i] = f * robust::ratio(pressureFromSesame(sieToSesame(DPDE)), x * R[i]);
        });
  }

  template <typename LambdaIndexer>
  inline void GruneisenParamFromDensityInternalEnergy(
      const Real *rhos, const Real *sies, Real *gm1s, Real *scratch, const int num,
      LambdaIndexer /*lambdas*/, Transform &&transform = Transform()) const {
    static auto const name =
        singularity::mfuncname::member_func_name(typeid(EOSPAC).name(), __func__);
    static auto const cname = name.c_str();
    using namespace EospacWrapper;
    EOS_REAL *R = const_cast<EOS_REAL *>(&rhos[0]);
    EOS_REAL *E = const_cast<EOS_REAL *>(&sies[0]);
    EOS_REAL *T = scratch + 0 * num;
    EOS_REAL *P = scratch + 1 * num;
    EOS_REAL *dx = scratch + 2 * num;
    EOS_REAL *DEDT = scratch + 3 * num;
    EOS_REAL *DPDT = scratch + 4 * num;
    EOS_REAL *Etmp = P;
    EOS_REAL *dy = DEDT;

    EOS_INTEGER table = TofRE_table_;
    {
      EOS_INTEGER options[3];
      EOS_REAL values[3];
      EOS_INTEGER nopts = 0;

      if (transform.x.is_set()) {
        options[nopts] = EOS_X_CONVERT;
        values[nopts] = 1.0 / transform.x.get();
        ++nopts;
      }

      options[nopts] = EOS_Y_CONVERT;
      values[nopts] = sieFromSesame(1.0);
      if (transform.y.is_set()) {
        values[nopts] /= transform.y.get();
      }
      ++nopts;
      eosSafeInterpolate(&table, num, R, E, T, dx, dy, "TofRE", Verbosity::Quiet, options,
                         values, nopts);
    }

    EOS_INTEGER options[1];
    EOS_REAL values[1];
    EOS_INTEGER nopts = 0;
    if (transform.x.is_set()) {
      options[nopts] = EOS_X_CONVERT;
      values[nopts] = 1.0 / transform.x.get();
    } else {
      options[nopts] = EOS_XY_PASSTHRU;
      values[nopts] = 1.0;
    }
    ++nopts;
    table = EofRT_table_;
    eosSafeInterpolate(&table, num, R, T, Etmp, dx, DEDT, "EofRT", Verbosity::Quiet,
                       options, values, nopts);

    table = PofRT_table_;
    eosSafeInterpolate(&table, num, R, T, P, dx, DPDT, "PofRT", Verbosity::Quiet, options,
                       values, nopts);

    const Real x = transform.x.is_set() ? transform.x.get() : 1.0;
    const Real f = transform.f.is_set() ? transform.f.get() : 1.0;

    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          const Real DPDE = DPDT[i] / DEDT[i];
          gm1s[i] = f * robust::ratio(pressureFromSesame(sieToSesame(DPDE)), x * R[i]);
        });
  }

  template <typename LambdaIndexer>
  inline void PTofRE(Real *rhos, Real *sies, Real *presses, Real *temps, Real *dpdrs,
                     Real *dpdes, Real *dtdrs, Real *dtdes, Real *scratch, const int num,
                     LambdaIndexer lambdas) const {
    static auto const name =
        singularity::mfuncname::member_func_name(typeid(EOSPAC).name(), __func__);
    static auto const cname = name.c_str();

    PressureFromDensityInternalEnergy(rhos, sies, presses, scratch, num, lambdas);
    TemperatureFromDensityInternalEnergy(rhos, sies, temps, scratch, num, lambdas);

    Real *drho = scratch;
    Real *de = scratch + num;
    Real *Pr = scratch + 2 * num;
    Real *Pe = scratch + 3 * num;
    Real *Tr = scratch + 4 * num;
    Real *Te = scratch + 5 * num;
    Real *rho_p_drho = scratch + 6 * num;
    Real *sie_p_de = scratch + 7 * num;

    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          drho[i] = rhos[i] * 1.0e-6;
          de[i] = sies[i] * 1.0e-6;
          rho_p_drho[i] = rhos[i] + drho[i];
          sie_p_de[i] = sies[i] + de[i];
        });

    Real *R = &rhos[0];
    Real *E = &sies[0];

    Real *internal_scratch = scratch + 8 * num;

    PressureFromDensityInternalEnergy(rho_p_drho, E, Pr, internal_scratch, num, lambdas);
    PressureFromDensityInternalEnergy(R, sie_p_de, Pe, internal_scratch, num, lambdas);
    TemperatureFromDensityInternalEnergy(rho_p_drho, E, Tr, internal_scratch, num,
                                         lambdas);
    TemperatureFromDensityInternalEnergy(R, sie_p_de, Te, internal_scratch, num, lambdas);

    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          dpdrs[i] = (Pr[i] - presses[i]) / drho[i];
          dpdes[i] = (Pe[i] - presses[i]) / de[i];
          dtdrs[i] = (Tr[i] - temps[i]) / drho[i];
          dtdes[i] =
              (Te[i] - temps[i]) /
              de[i]; // Would it be better to skip the calculation of Te and return 1/cv?
        });
  }

  SG_PIF_NOWARN
  PORTABLE_INLINE_FUNCTION void PTofRE(Real &rho, Real &sie, Real *lambda, Real &press,
                                       Real &temp, Real &dpdr, Real &dpde, Real &dtdr,
                                       Real &dtde) const {

    press = PressureFromDensityInternalEnergy(rho, sie, lambda);
    temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
    const Real drho = rho * 1.0e-6;
    const Real de = sie * 1.0e-6;
    const Real Pr = PressureFromDensityInternalEnergy(rho + drho, sie, lambda);
    const Real Pe = PressureFromDensityInternalEnergy(rho, sie + de, lambda);
    const Real Tr = TemperatureFromDensityInternalEnergy(rho + drho, sie, lambda);
    const Real Te = TemperatureFromDensityInternalEnergy(rho, sie + de, lambda);
    dpdr = (Pr - press) / drho;
    dpde = (Pe - press) / de;
    dtdr = (Tr - temp) / drho;
    dtde = (Te - temp) /
           de; // Would it be better to skip the calculation of Te and return 1/cv?
    return;
  }

  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    auto nbuffers = scratch_nbuffers();
    if (nbuffers.find(method) != nbuffers.end()) {
      return sizeof(Real) * nbuffers[method] * nelements;
    }
    return 0;
  }

  static inline unsigned long max_scratch_size(unsigned int nelements) {
    auto nbuffers = scratch_nbuffers();
    auto e = std::max_element(
        nbuffers.begin(), nbuffers.end(),
        [](const auto &a, const auto &b) { return a.second < b.second; });
    return scratch_size(e->first, nelements);
  }

  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  PORTABLE_INLINE_FUNCTION int nlambda() const noexcept { return 0; }
  inline void Finalize() {
    using namespace EospacWrapper;
    eosSafeDestroy(NT, tablehandle, Verbosity::Quiet);
  }
  static std::string EosType() { return std::string("EOSPAC"); }
  static std::string EosPyType() { return EosType(); }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    printf("EOSPAC parameters:\nmatid = %i\n", matid_);
  }
  PORTABLE_FORCEINLINE_FUNCTION Real MinimumDensity() const { return rho_min_; }
  PORTABLE_FORCEINLINE_FUNCTION Real MinimumTemperature() const { return temp_min_; }

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
  // EOS_INTEGER PofRE_table_;
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
  Real rho_min_ = 0;
  Real temp_min_ = 0;

  static inline std::map<std::string, unsigned int> &scratch_nbuffers() {
    static std::map<std::string, unsigned int> nbuffers = {
        {"TemperatureFromDensityInternalEnergy", 2},
        {"PressureFromDensityTemperature", 2},
        {"FillEos", 4},
        {"InternalEnergyFromDensityTemperature", 2},
        {"PressureFromDensityInternalEnergy", 3},
        {"SpecificHeatFromDensityTemperature", 2},
        {"SpecificHeatFromDensityInternalEnergy", 4},
        {"BulkModulusFromDensityTemperature", 6},
        {"BulkModulusFromDensityInternalEnergy", 6},
        {"GruneisenParamFromDensityTemperature", 4},
        {"GruneisenParamFromDensityInternalEnergy", 5},
        {"PTofRE", 11}};
    return nbuffers;
  }
};

// ======================================================================
// Implementation details below
// ======================================================================

inline EOSPAC::EOSPAC(const int matid, bool invert_at_setup) : matid_(matid) {
  using namespace EospacWrapper;
  EOS_INTEGER tableType[NT] = {EOS_Pt_DT, EOS_T_DUt, EOS_Ut_DT, EOS_D_PtT, EOS_T_DPt};
  eosSafeLoad(NT, matid, tableType, tablehandle,
              std::vector<std::string>(
                  {"EOS_Pt_DT", "EOS_T_DUt", "EOS_Ut_DT", "EOS_D_PtT", "EOS_T_DPt"}),
              Verbosity::Quiet, invert_at_setup);
  PofRT_table_ = tablehandle[0];
  TofRE_table_ = tablehandle[1];
  EofRT_table_ = tablehandle[2];
  RofPT_table_ = tablehandle[3];
  TofRP_table_ = tablehandle[4];

  // Set reference states and table bounds
  SesameMetadata m;
  eosGetMetadata(matid, m, Verbosity::Quiet);
  rho_ref_ = m.normalDensity;
  rho_min_ = m.rhoMin;
  temp_min_ = m.TMin;

  EOS_REAL R[1] = {rho_ref_};
  EOS_REAL T[1] = {temperatureToSesame(temp_ref_)};
  EOS_REAL E[1], P[1], dx[1], dy[1];
  EOS_REAL DEDR, DEDT, DPDR, DPDT, DPDE, BMOD;
  EOS_INTEGER nxypairs = 1;
  eosSafeInterpolate(&EofRT_table_, nxypairs, R, T, E, dx, dy, "EofRT", Verbosity::Quiet);
  DEDR = dx[0];
  DEDT = dy[0];
  sie_ref_ = sieFromSesame(E[0]);
  cv_ref_ = cvFromSesame(std::max(DEDT, 0.0));

  eosSafeInterpolate(&PofRT_table_, nxypairs, R, T, P, dx, dy, "PofRT", Verbosity::Quiet);
  Real PRESS = P[0];
  press_ref_ = pressureFromSesame(PRESS);
  DPDR = dx[0];
  DPDT = dy[0];
  DPDE = robust::ratio(DPDT, DEDT);
  BMOD = rho_ref_ * DPDR + DPDE * (robust::ratio(PRESS, rho_ref_) - rho_ref_ * DEDR);
  bmod_ref_ = bulkModulusFromSesame(std::max(BMOD, 0.0));
  dpde_ref_ = pressureFromSesame(sieToSesame(DPDE));
  dvdt_ref_ =
      robust::ratio(dpde_ref_ * cv_ref_, rho_ref_ * rho_ref_ * pressureFromSesame(DPDR));
}

SG_PIF_NOWARN
PORTABLE_INLINE_FUNCTION Real EOSPAC::TemperatureFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  using namespace EospacWrapper;
  EOS_REAL R[1] = {rho}, E[1] = {sieToSesame(sie)}, T[1], dTdr[1], dTde[1];
  EOS_INTEGER nxypairs = 1;
  EOS_INTEGER table = TofRE_table_;
  eosSafeInterpolate(&table, nxypairs, R, E, T, dTdr, dTde, "TofRE", Verbosity::Quiet);
  return Real(temperatureFromSesame(T[0]));
}

SG_PIF_NOWARN
PORTABLE_INLINE_FUNCTION Real EOSPAC::PressureFromDensityTemperature(const Real rho,
                                                                     const Real temp,
                                                                     Real *lambda) const {
  using namespace EospacWrapper;
  EOS_REAL R[1] = {rho}, P[1], T[1] = {temperatureToSesame(temp)}, dPdr[1], dPdT[1];
  EOS_INTEGER nxypairs = 1;
  EOS_INTEGER table = PofRT_table_;
  eosSafeInterpolate(&table, nxypairs, R, T, P, dPdr, dPdT, "PofRT", Verbosity::Quiet);
  return Real(pressureFromSesame(P[0]));
}

SG_PIF_NOWARN
PORTABLE_INLINE_FUNCTION Real EOSPAC::EntropyFromDensityTemperature(
    const Real rho, const Real temperature, Real *lambda) const {
  EntropyIsNotEnabled("EOSPAC");
  return 1.0;
}

SG_PIF_NOWARN
PORTABLE_INLINE_FUNCTION void EOSPAC::FillEos(Real &rho, Real &temp, Real &sie,
                                              Real &press, Real &cv, Real &bmod,
                                              const unsigned long output,
                                              Real *lambda) const {
  using namespace EospacWrapper;
  EOS_REAL R[1] = {rho}, T[1] = {temperatureToSesame(temp)};
  EOS_REAL E[1] = {sie}, P[1] = {pressureToSesame(press)};
  EOS_REAL dx[1], dy[1];
  EOS_INTEGER nxypairs = 1;
  Real /*CV,*/ BMOD_T, BMOD, SIE, PRESS, DPDE, DPDT, DPDR, DEDT, DEDR;
  const unsigned long input = ~output;
  if (output == thermalqs::none) {
    UNDEFINED_ERROR;
  }
  if (output & thermalqs::density) {
    if (input & thermalqs::pressure && input & thermalqs::temperature) {
      EOS_INTEGER table = RofPT_table_;
      eosSafeInterpolate(&table, nxypairs, P, T, R, dx, dy, "RofPT", Verbosity::Quiet);
      rho = R[0];
    } else {
      UNDEFINED_ERROR;
    }
  }
  if (output & thermalqs::temperature) {
    if (input & thermalqs::density && input & thermalqs::specific_internal_energy) {
      temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
    } else if (input & thermalqs::density && input & thermalqs::pressure) {
      EOS_INTEGER table = TofRP_table_;
      eosSafeInterpolate(&table, nxypairs, R, P, T, dx, dy, "TofRP", Verbosity::Quiet);
      temp = temperatureFromSesame(T[0]);
    } else {
      UNDEFINED_ERROR;
    }
  }
  if ((output & thermalqs::specific_internal_energy) ||
      (output & thermalqs::specific_heat || output & thermalqs::bulk_modulus)) {
    EOS_INTEGER table = EofRT_table_;
    eosSafeInterpolate(&table, nxypairs, R, T, E, dx, dy, "EofRT", Verbosity::Quiet);
    SIE = E[0];
    DEDR = dx[0];
    DEDT = dy[0];
  }
  if ((output & thermalqs::pressure) || (output & thermalqs::bulk_modulus)) {
    EOS_INTEGER table = PofRT_table_;
    eosSafeInterpolate(&table, nxypairs, R, T, P, dx, dy, "PofRT", Verbosity::Quiet);
    PRESS = P[0];
    DPDR = dx[0];
    DPDT = dy[0];
    // Thermodynamics: Bt = rho*dP/drho|T, Bs=Bt*Cv/Cp, Cv=dE/dT|v, and
    // Cp-Cv=-T(dP/dT|v)^2/(dP/dV|T) or Cp-Cv=-T/rho^2 (dP/dT|v)^2/(dP/drho|T)
    // Therefore: Bs=Bt Cv/(CV+above)
    // BMOD = rho*dx[0]*CV/(CV-T[0]/(rho*rho)*dy[0]*dy[0]/dx[0]);
  }
  if (output & thermalqs::specific_internal_energy) {
    sie = sieFromSesame(SIE);
  }
  if (output & thermalqs::pressure) {
    press = pressureFromSesame(PRESS);
  }
  if (output & thermalqs::specific_heat) {
    cv = cvFromSesame(std::max(DEDT, 0.0)); // Here we do something to the data!
  }
  if (output & thermalqs::bulk_modulus) {
    BMOD_T = std::max(rho * DPDR, 0.0);
    if (DEDT > 0.0 && rho > 0.0) {
      DPDE = DPDT / DEDT;
      BMOD = rho * DPDR + DPDE * (PRESS / rho - rho * DEDR);
    } else if (rho > 0.0) { // Case: DEDT <= 0
      // We need a different DPDE call apparently????
      // In xRAGE, they call out to P(rho,e) in this case to get the
      // derivative directly from the half-inverted table.
      // But upon further review,
      // I think it will end up evaluating to BMOD_T in any case because cv will
      // be zero! See eos_eospac.f90 line 1261 BMOD =
      // BMOD_T+DPDE*DPDE*std::max(DEDT,0.0)*T[0]/rho;
      BMOD = BMOD_T; //+DPDE*DPDE*std::max(DEDT,0.0)*T[0]/rho;
    } else {
      BMOD = 0.0;
    }
    bmod = bulkModulusFromSesame(std::max(BMOD, 0.0));
  }
}

SG_PIF_NOWARN
PORTABLE_INLINE_FUNCTION Real EOSPAC::InternalEnergyFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  using namespace EospacWrapper;
  Real RHO = rho, TEMP = temp, sie, press, cv, bmod;
  const unsigned long output = thermalqs::specific_internal_energy;
  FillEos(RHO, TEMP, sie, press, cv, bmod, output, lambda);
  return sie;
}
SG_PIF_NOWARN
PORTABLE_INLINE_FUNCTION Real EOSPAC::BulkModulusFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  using namespace EospacWrapper;
  Real RHO = rho, TEMP = temp, sie, press, cv, bmod;
  const unsigned long output = thermalqs::bulk_modulus;
  FillEos(RHO, TEMP, sie, press, cv, bmod, output, lambda);
  return bmod;
}
SG_PIF_NOWARN
PORTABLE_INLINE_FUNCTION Real EOSPAC::SpecificHeatFromDensityTemperature(
    const Real rho, const Real temp, Real *lambda) const {
  using namespace EospacWrapper;
  Real RHO = rho, TEMP = temp, sie, press, cv, bmod;
  const unsigned long output = thermalqs::specific_heat;
  FillEos(RHO, TEMP, sie, press, cv, bmod, output, lambda);
  return cv;
}
SG_PIF_NOWARN
PORTABLE_INLINE_FUNCTION Real EOSPAC::PressureFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  using namespace EospacWrapper;
  Real temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  return PressureFromDensityTemperature(rho, temp, lambda);
}
SG_PIF_NOWARN
PORTABLE_INLINE_FUNCTION Real EOSPAC::EntropyFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  using namespace EospacWrapper;
  const Real temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  return EntropyFromDensityTemperature(rho, temp, lambda);
}
SG_PIF_NOWARN
PORTABLE_INLINE_FUNCTION Real EOSPAC::SpecificHeatFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  using namespace EospacWrapper;
  Real temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  return SpecificHeatFromDensityTemperature(rho, temp, lambda);
}
SG_PIF_NOWARN
PORTABLE_INLINE_FUNCTION Real EOSPAC::BulkModulusFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  using namespace EospacWrapper;
  Real temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  return BulkModulusFromDensityTemperature(rho, temp, lambda);
}

SG_PIF_NOWARN
PORTABLE_INLINE_FUNCTION Real EOSPAC::GruneisenParamFromDensityTemperature(
    const Real rho, const Real temperature, Real *lambda) const {
  using namespace EospacWrapper;
  EOS_REAL R[1] = {rho}, T[1] = {temperatureToSesame(temperature)};
  EOS_REAL E[1], P[1], dx[1], dy[1];
  EOS_INTEGER nxypairs = 1;
  Real DPDT, DEDT, DPDE;
  EOS_INTEGER table = EofRT_table_;
  eosSafeInterpolate(&table, nxypairs, R, T, E, dx, dy, "EofRT", Verbosity::Quiet);
  DEDT = dy[0];
  table = PofRT_table_;
  eosSafeInterpolate(&table, nxypairs, R, T, P, dx, dy, "PofRT", Verbosity::Quiet);
  DPDT = dy[0];
  DPDE = DPDT / DEDT;
  return robust::ratio(pressureFromSesame(sieToSesame(DPDE)), rho);
}
SG_PIF_NOWARN
PORTABLE_INLINE_FUNCTION Real EOSPAC::GruneisenParamFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  Real temperature = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  return GruneisenParamFromDensityTemperature(rho, temperature, lambda);
}

SG_PIF_NOWARN
PORTABLE_INLINE_FUNCTION void
EOSPAC::DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                             Real *lambda, Real &rho, Real &sie) const {
  using namespace EospacWrapper;
  EOS_REAL P[1] = {pressureToSesame(press)};
  EOS_REAL T[1] = {temperatureToSesame(temp)};
  EOS_REAL dx[1], dy[1], R[1], E[1];
  EOS_INTEGER nxypairs = 1;
  EOS_INTEGER table;

  table = RofPT_table_;
  eosSafeInterpolate(&table, nxypairs, P, T, R, dx, dy, "RofPT", Verbosity::Quiet);
  rho = R[0];

  table = EofRT_table_;
  eosSafeInterpolate(&table, nxypairs, R, T, E, dx, dy, "EofPRT", Verbosity::Quiet);
  sie = sieFromSesame(E[0]);
}

SG_PIF_NOWARN
PORTABLE_INLINE_FUNCTION void
EOSPAC::ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                               Real &bmod, Real &dpde, Real &dvdt, Real *lambda) const {
  using namespace EospacWrapper;
  rho = rho_ref_;
  temp = temp_ref_;
  sie = sie_ref_;
  press = press_ref_;
  cv = cv_ref_;
  bmod = bmod_ref_;
  dpde = dpde_ref_;
  dvdt = dvdt_ref_;
}

} // namespace singularity

#endif // SINGULARITY_USE_EOSPAC
#endif // _SINGULARITY_EOS_EOS_EOS_EOSPAC_HPP_
