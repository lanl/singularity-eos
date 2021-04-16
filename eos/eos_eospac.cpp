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

#ifdef SINGULARITY_USE_EOSPAC

#include "eos.hpp"
#include "eos_eospac.hpp"

inline Real TemperatureToSesame(const Real CodeTemp) { return CodeTemp; }
inline Real TemperatureFromSesame(const Real SesTemp) { return SesTemp; }
inline Real PressureFromSesame(const Real SesPress) { return 1e10 * SesPress; }
inline Real PressureToSesame(const Real CodePress) { return 1e-10 * CodePress; }
inline Real SieToSesame(const Real CodeSie) { return 1e-10 * CodeSie; }
inline Real SieFromSesame(const Real SesSie) { return 1e10 * SesSie; }
inline Real CvFromSesame(const Real SesCv) { return 1e10 * SesCv; }
inline Real BulkModulusFromSesame(const Real SesBmod) { return 1e10 * SesBmod; }

namespace singularity {

// How do I make a destructor? How do I free the EOS memory if more than one
// points to the same table? Does EOSPAC give me multiple (reference-counted)
// handles to the same table or read in the table multiple times?

EOSPAC::EOSPAC(const int matid) : matid_(matid) {
  const int NT = 5;
  EOS_INTEGER ntables = NT, errorCode = EOS_OK, tablehandle[NT];
  EOS_INTEGER tableType[NT] = {EOS_Pt_DT, EOS_T_DUt, EOS_Ut_DT, EOS_D_PtT, Eos_Info};
  EOS_INTEGER MATID[NT];
  for (int i = 0; i < NT; i++) {
    MATID[i] = matid;
  }
  eos_CreateTables(&ntables, tableType, MATID, tablehandle, &errorCode);
  // eos_SetOption(&tablehandle[0],&THE_OPTION,EOS_NullPtr,&errorCode);
  eos_LoadTables(&ntables, tablehandle, &errorCode);
  PofRT_table_ = tablehandle[0];
  TofRE_table_ = tablehandle[1];
  EofRT_table_ = tablehandle[2];
  RofPT_table_ = tablehandle[3];
  EOS_Info_table_ = tablehandle[4];

  // Set reference states
  EOS_INTEGER error_code = EOS_OK;
  EOS_INTEGER infoItems[1] = {EOS_Normal_Density};
  eos_GetTableInfo(EOS_Info_table_, 1, infoItems, &rho_ref_, &error_code);
  FillEos(rho_ref_, temp_ref_, sie_ref_, press_ref_, cv_ref_, bmod_ref_,
          thermalqs::specific_internal_energy|thermalqs::pressure|thermalqs::specific_heat|thermalqs::bulk_modulus);
}
Real EOSPAC::TemperatureFromDensityInternalEnergy(const Real rho,
                                                  const Real sie,
                                                  Real *lambda) {
  EOS_REAL R[1] = {rho}, E[1] = {SieToSesame(sie)}, T[1], dTdr[1], dTde[1];
  EOS_INTEGER nxypairs = 1, errorCode;
  eos_Interpolate(&TofRE_table_, &nxypairs, R, E, T, dTdr, dTde, &errorCode);
  return Real(TemperatureFromSesame(T[0]));
}
Real EOSPAC::PressureFromDensityTemperature(const Real rho, const Real temp,
                                            Real *lambda) {
  EOS_REAL R[1] = {rho}, P[1], T[1] = {TemperatureToSesame(temp)}, dPdr[1],
           dPdT[1];
  EOS_INTEGER nxypairs = 1, errorCode;
  eos_Interpolate(&PofRT_table_, &nxypairs, R, T, P, dPdr, dPdT, &errorCode);
  return Real(PressureFromSesame(P[0]));
}
void EOSPAC::FillEos(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                     Real &bmod, const unsigned long output, Real *lambda) {
  EOS_REAL R[1] = {rho}, T[1] = {TemperatureToSesame(temp)}, E[1], P[1], dx[1],
           dy[1];
  EOS_INTEGER nxypairs = 1, errorCode = EOS_OK;
  Real CV, BMOD_T, BMOD, SIE, PRESS, DPDE, DPDT, DPDR, DEDT, DEDR;
  if ((output & EOS::specific_internal_energy) ||
      (output & EOS::specific_heat || output & EOS::bulk_modulus)) {
    eos_Interpolate(&EofRT_table_, &nxypairs, R, T, E, dx, dy, &errorCode);
    SIE = E[0];
    DEDR = dx[0];
    DEDT = dy[0];
  }
  if ((output & EOS::pressure) || (output & EOS::bulk_modulus)) {
    eos_Interpolate(&PofRT_table_, &nxypairs, R, T, P, dx, dy, &errorCode);
    PRESS = P[0];
    DPDR = dx[0];
    DPDT = dy[0];
    // Thermodynamics: Bt = rho*dP/drho|T, Bs=Bt*Cv/Cp, Cv=dE/dT|v, and
    // Cp-Cv=-T(dP/dT|v)^2/(dP/dV|T) or Cp-Cv=-T/rho^2 (dP/dT|v)^2/(dP/drho|T)
    // Therefore: Bs=Bt Cv/(CV+above)
    // BMOD = rho*dx[0]*CV/(CV-T[0]/(rho*rho)*dy[0]*dy[0]/dx[0]);
  }
  if (output & EOS::specific_internal_energy) {
    sie = SieFromSesame(SIE);
  }
  if (output & EOS::pressure) {
    press = PressureFromSesame(PRESS);
  }
  if (output & EOS::specific_heat) {
    cv = CvFromSesame(std::max(DEDT, 0.0)); // Here we do something to the data!
  }
  if (output & EOS::bulk_modulus) {
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
    bmod = BulkModulusFromSesame(std::max(BMOD, 0.0));
  }
}

Real EOSPAC::InternalEnergyFromDensityTemperature(const Real rho,
                                                  const Real temp,
                                                  Real *lambda) {
  Real RHO = rho, TEMP = temp, sie, press, cv, bmod;
  const unsigned long output = EOS::specific_internal_energy;
  FillEos(RHO, TEMP, sie, press, cv, bmod, output, lambda);
  return sie;
}
Real EOSPAC::BulkModulusFromDensityTemperature(const Real rho, const Real temp,
                                               Real *lambda) {
  Real RHO = rho, TEMP = temp, sie, press, cv, bmod;
  const unsigned long output = EOS::bulk_modulus;
  FillEos(RHO, TEMP, sie, press, cv, bmod, output, lambda);
  return bmod;
}
Real EOSPAC::SpecificHeatFromDensityTemperature(const Real rho, const Real temp,
                                                Real *lambda) {
  Real RHO = rho, TEMP = temp, sie, press, cv, bmod;
  const unsigned long output = EOS::specific_heat;
  FillEos(RHO, TEMP, sie, press, cv, bmod, output, lambda);
  return cv;
}
Real EOSPAC::PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                               Real *lambda) {
  Real temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  return PressureFromDensityTemperature(rho, temp, lambda);
}
Real EOSPAC::SpecificHeatFromDensityInternalEnergy(const Real rho,
                                                   const Real sie,
                                                   Real *lambda) {
  Real temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  return SpecificHeatFromDensityTemperature(rho, temp, lambda);
}
Real EOSPAC::BulkModulusFromDensityInternalEnergy(const Real rho,
                                                  const Real sie,
                                                  Real *lambda) {
  Real temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
  return BulkModulusFromDensityTemperature(rho, temp, lambda);
}
void EOSPAC::PrintParams() { std::cout << "Matid: " << matid_ << std::endl; }

void EOSPAC::DensityEnergyFromPressureTemperature(const Real press,
                                                  const Real temp, Real *lambda,
                                                  Real &rho, Real &sie) {
  EOS_REAL P[1] = PressureToSesame(press), T[1] = TemperatureToSesame(temp),
           R[1], dx[1], dy[1];
  EOS_INTEGER nxypairs = 1, errorCode = EOS_OK;
  eos_Interpolate(&_DofPT_table, &nxypairs, P, T, R, dx, dy, &errorCode);
  rho = R[0];
  sie = InternalEnergyFromDensityTemperature(rho, temp, lambda);
}

void EOSPAC::PTofRE(const Real rho, const Real sie, Real *lambda, Real &press,
                    Real &temp, Real &dpdr, Real &dpde, Real &dtdr,
                    Real &dtde) {
  EOS_REAL R[1] = {rho}, T[1], E[1] = {SieToSesame(sie)}, P[1], dx[1], dy[1];
  EOS_INTEGER nxypairs = 1, errorCode = EOS_OK;
  Real DTDR_E, DTDE_R, DPDR_T, DPDT_R;

  eos_Interpolate(&TofRE_table_, &nxypairs, R, E, T, dx, dy, &errorCode);
  DTDR_E = dx[0];
  DTDE_R = dy[0];
  eos_Interpolate(&PofRT_table_, &nxypairs, R, T, P, dx, dy, &errorCode);
  DPDR_T = dx[0];
  DPDT_R = dy[0];

  press = PressureFromSesame(P[0]);
  temp = TemperatureFromSesame(T[0]);
  dtdr = TemperatureFromSesame(DTDR_E);
  dtde = SieToSesame(TemperatureFromSesame(DTDE_R));
  dpde = SieToSesame(PressureFromSesame(DPDT_R * DTDE_R));
  dpdr = PressureFromSesame(DRDR_T + DTDR_E * DPDT_R);
}

} // namespace singularity

#endif // SINGULARITY_USE_EOSPAC
