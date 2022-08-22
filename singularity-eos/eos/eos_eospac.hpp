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

#ifndef _SINGULARITY_EOS_EOS_EOS_EOSPAC_HPP_
#define _SINGULARITY_EOS_EOS_EOS_EOSPAC_HPP_
#ifdef SINGULARITY_USE_EOSPAC

#include <cstdio>
#include <string>
#include <vector>

#include <eos_Interface.h>

#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/eos/eos_base.hpp>
#include <eospac-wrapper/eospac_wrapper.hpp>

namespace singularity {

using namespace eos_base;

// How do I make a destructor? How do I free the EOS memory if more than one
// points to the same table? Does EOSPAC give me multiple (reference-counted)
// handles to the same table or read in the table multiple times?

// Only really works in serial
// Not really supported on device
class EOSPAC : public EosBase<EOSPAC> {
 public:
  inline EOSPAC() = default;
  inline EOSPAC(int matid, bool invert_at_setup = false);
  inline EOSPAC GetOnDevice() { return *this; }
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(const Real rho,
                                                        const Real temperature,
                                                        Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                                           Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(const Real rho,
                                                            const Real temperature,
                                                            Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(const Real rho,
                                                           const Real temperature,
                                                           Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                 Real &cv, Real &bmod, const unsigned long output,
                                 Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION
  void DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                            Real *lambda, Real &rho, Real &sie) const;
  PORTABLE_INLINE_FUNCTION void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie,
                                                Real &press, Real &cv, Real &bmod,
                                                Real &dpde, Real &dvdt,
                                                Real *lambda = nullptr) const;

  // Generic functions provided by the base class. These contain
  // e.g. the vector overloads that use the scalar versions declared
  // here We explicitly list, rather than using the macro because we
  // overload some methods.
  using EosBase<EOSPAC>::TemperatureFromDensityInternalEnergy;
  using EosBase<EOSPAC>::InternalEnergyFromDensityTemperature;
  using EosBase<EOSPAC>::PressureFromDensityTemperature;
  using EosBase<EOSPAC>::PressureFromDensityInternalEnergy;
  using EosBase<EOSPAC>::SpecificHeatFromDensityTemperature;
  using EosBase<EOSPAC>::SpecificHeatFromDensityInternalEnergy;
  using EosBase<EOSPAC>::BulkModulusFromDensityTemperature;
  using EosBase<EOSPAC>::BulkModulusFromDensityInternalEnergy;
  using EosBase<EOSPAC>::GruneisenParamFromDensityTemperature;
  using EosBase<EOSPAC>::GruneisenParamFromDensityInternalEnergy;
  using EosBase<EOSPAC>::PTofRE;
  using EosBase<EOSPAC>::FillEos;

  // TODO (JHP): Change EOSPAC vector implementations to be more performant
  // template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  // inline
  // void TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
  //                                           ConstRealIndexer &&sies,
  //                                           RealIndexer &&temperatures,
  //                                           const int num,
  //                                           LambdaIndexer &&lambdas) const;
  // template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  // inline
  // void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
  //                                           ConstRealIndexer &&temperatures,
  //                                           RealIndexer &&sies,
  //                                           const int num,
  //                                           LambdaIndexer &&lambdas) const;
  // template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  // inline
  // void PressureFromDensityTemperature(ConstRealIndexer &&rhos,
  //                                     ConstRealIndexer &&temperatures,
  //                                     RealIndexer &&pressures,
  //                                     const int num,
  //                                     LambdaIndexer &&lambdas) const;
  // template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  // inline
  // void PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
  //                                        ConstRealIndexer &&sies,
  //                                        RealIndexer &&pressures,
  //                                        const int num,
  //                                        LambdaIndexer &&lambdas) const;
  // template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  // inline
  // void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
  //                                         ConstRealIndexer &&temperatures,
  //                                         RealIndexer &&cvs,
  //                                         const int num,
  //                                         LambdaIndexer &&lambdas) const;
  // template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  // inline
  // void SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos,
  //                                            ConstRealIndexer &&sies,
  //                                            RealIndexer &&cvs,
  //                                            const int num,
  //                                            LambdaIndexer &&lambdas) const;
  // template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  // inline
  // void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
  //                                        ConstRealIndexer &&temperatures,
  //                                        RealIndexer &&bmods,
  //                                        const int num,
  //                                        LambdaIndexer &&lambdas) const;
  // template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  // inline
  // void BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos,
  //                                           ConstRealIndexer &&sies,
  //                                           RealIndexer &&bmods,
  //                                           const int num,
  //                                           LambdaIndexer &&lambdas) const;

  // template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  // inline
  // void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
  //                                           ConstRealIndexer &&temperatures,
  //                                           RealIndexer &&gm1s,
  //                                           const int num,
  //                                           LambdaIndexer &&lambdas) const;

  // template<typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  // inline
  // void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
  //                                              ConstRealIndexer &&sies,
  //                                              RealIndexer &&gm1s,
  //                                              const int num,
  //                                              LambdaIndexer &&lambdas) const;

  // template<typename RealIndexer, typename LambdaIndexer>
  // inline
  // void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
  //              RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
  //              const int num, const unsigned long output,
  //              LambdaIndexer &&lambdas) const;
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  PORTABLE_INLINE_FUNCTION int nlambda() const noexcept { return 0; }
  inline void Finalize() {}
  static std::string EosType() { return std::string("EOSPAC"); }
  inline void PrintParams() const { printf("EOSPAC parameters:\nmatid = %s\n", matid_); }
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
  DPDE = DPDT / (DEDT + EPS);
  BMOD = rho_ref_ * DPDR + DPDE * (PRESS / (rho_ref_ + EPS) - rho_ref_ * DEDR);
  bmod_ref_ = bulkModulusFromSesame(std::max(BMOD, 0.0));
  dpde_ref_ = pressureFromSesame(sieToSesame(DPDE));
  dvdt_ref_ =
      dpde_ref_ * cv_ref_ / (rho_ref_ * rho_ref_ * pressureFromSesame(DPDR) + EPS);
}

PORTABLE_INLINE_FUNCTION Real EOSPAC::TemperatureFromDensityInternalEnergy(const Real rho,
                                                                    const Real sie,
                                                                    Real *lambda) const {
  using namespace EospacWrapper;
  EOS_REAL R[1] = {rho}, E[1] = {sieToSesame(sie)}, T[1], dTdr[1], dTde[1];
  EOS_INTEGER nxypairs = 1;
  EOS_INTEGER table = TofRE_table_;
  eosSafeInterpolate(&table, nxypairs, R, E, T, dTdr, dTde, "TofRE", Verbosity::Quiet);
  return Real(temperatureFromSesame(T[0]));
}

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

PORTABLE_INLINE_FUNCTION void EOSPAC::FillEos(Real &rho, Real &temp, Real &sie, Real &press,
                                       Real &cv, Real &bmod, const unsigned long output,
                                       Real *lambda) const {
using namespace EospacWrapper;
  EOS_REAL R[1] = {rho}, T[1] = {temperatureToSesame(temp)};
  EOS_REAL E[1] = {sie}, P[1] = {pressureToSesame(press)};
  EOS_REAL dx[1], dy[1];
  EOS_INTEGER nxypairs = 1;
  Real CV, BMOD_T, BMOD, SIE, PRESS, DPDE, DPDT, DPDR, DEDT, DEDR;
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

PORTABLE_INLINE_FUNCTION Real EOSPAC::InternalEnergyFromDensityTemperature(const Real rho,
                                                                    const Real temp,
                                                                    Real *lambda) const {
using namespace EospacWrapper;
  Real RHO = rho, TEMP = temp, sie, press, cv, bmod;
  const unsigned long output = thermalqs::specific_internal_energy;
  FillEos(RHO, TEMP, sie, press, cv, bmod, output, lambda);
  return sie;
}
PORTABLE_INLINE_FUNCTION Real EOSPAC::BulkModulusFromDensityTemperature(const Real rho,
                                                                 const Real temp,
                                                                 Real *lambda) const {
using namespace EospacWrapper;
  Real RHO = rho, TEMP = temp, sie, press, cv, bmod;
  const unsigned long output = thermalqs::bulk_modulus;
  FillEos(RHO, TEMP, sie, press, cv, bmod, output, lambda);
  return bmod;
}
PORTABLE_INLINE_FUNCTION Real EOSPAC::SpecificHeatFromDensityTemperature(const Real rho,
                                                                  const Real temp,
                                                                  Real *lambda) const {
using namespace EospacWrapper;
  Real RHO = rho, TEMP = temp, sie, press, cv, bmod;
  const unsigned long output = thermalqs::specific_heat;
  FillEos(RHO, TEMP, sie, press, cv, bmod, output, lambda);
  return cv;
}
PORTABLE_INLINE_FUNCTION Real EOSPAC::PressureFromDensityInternalEnergy(const Real rho,
                                                                 const Real sie,
                                                                 Real *lambda) const {
using namespace EospacWrapper;
  Real temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  return PressureFromDensityTemperature(rho, temp, lambda);
}
PORTABLE_INLINE_FUNCTION Real EOSPAC::SpecificHeatFromDensityInternalEnergy(const Real rho,
                                                                     const Real sie,
                                                                     Real *lambda) const {
using namespace EospacWrapper;
  Real temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  return SpecificHeatFromDensityTemperature(rho, temp, lambda);
}
PORTABLE_INLINE_FUNCTION Real EOSPAC::BulkModulusFromDensityInternalEnergy(const Real rho,
                                                                    const Real sie,
                                                                    Real *lambda) const {
using namespace EospacWrapper;
  Real temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  return BulkModulusFromDensityTemperature(rho, temp, lambda);
}

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
  return pressureFromSesame(sieToSesame(DPDE)) / (rho + EPS);
}
PORTABLE_INLINE_FUNCTION Real EOSPAC::GruneisenParamFromDensityInternalEnergy(
    const Real rho, const Real sie, Real *lambda) const {
  Real temperature = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  return GruneisenParamFromDensityTemperature(rho, temperature, lambda);
}

PORTABLE_INLINE_FUNCTION
void EOSPAC::DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                                  Real *lambda, Real &rho,
                                                  Real &sie) const {
using namespace EospacWrapper;
  EOS_REAL P[1] = {pressureToSesame(press)};
  EOS_REAL T[1] = {temperatureToSesame(temp)};
  EOS_REAL dx[1], dy[1], R[1], E[1];
  EOS_INTEGER nxypairs = 1;
  EOS_INTEGER errorCode;
  EOS_INTEGER table;

  table = RofPT_table_;
  eosSafeInterpolate(&table, nxypairs, P, T, R, dx, dy, "RofPT", Verbosity::Quiet);
  rho = R[0];

  table = EofRT_table_;
  eosSafeInterpolate(&table, nxypairs, R, T, E, dx, dy, "EofPRT", Verbosity::Quiet);
  sie = sieFromSesame(E[0]);
}

PORTABLE_INLINE_FUNCTION void EOSPAC::ValuesAtReferenceState(Real &rho, Real &temp, Real &sie,
                                                      Real &press, Real &cv, Real &bmod,
                                                      Real &dpde, Real &dvdt,
                                                      Real *lambda) const {
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
