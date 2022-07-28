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

#include <cmath>
#include <singularity-eos/eos/eos.hpp>

namespace singularity {

PORTABLE_FUNCTION Real JWL::ReferencePressure(const Real rho) const {
  const Real x = _rho0 / rho;
  return _A * std::exp(-_R1 * x) + _B * std::exp(-_R2 * x);
}
PORTABLE_FUNCTION Real JWL::ReferenceEnergy(const Real rho) const {
  const Real x = _rho0 / rho;
  return _A / (_rho0 * _R1) * std::exp(-_R1 * x) +
         _B / (_rho0 * _R2) * std::exp(-_R2 * x);
}
PORTABLE_FUNCTION Real JWL::InternalEnergyFromDensityTemperature(const Real rho,
                                                                 const Real temp,
                                                                 Real *lambda) const {
  return ReferenceEnergy(rho) + _Cv * temp;
}
PORTABLE_FUNCTION Real JWL::PressureFromDensityInternalEnergy(const Real rho,
                                                              const Real sie,
                                                              Real *lambda) const {
  return ReferencePressure(rho) + _w * rho * (sie - ReferenceEnergy(rho));
}
PORTABLE_FUNCTION Real JWL::TemperatureFromDensityInternalEnergy(const Real rho,
                                                                 const Real sie,
                                                                 Real *lambda) const {
  return (sie - ReferenceEnergy(rho)) / _Cv;
}
PORTABLE_FUNCTION Real JWL::SpecificHeatFromDensityInternalEnergy(const Real rho,
                                                                  const Real sie,
                                                                  Real *lambda) const {
  return _Cv;
}
PORTABLE_FUNCTION Real JWL::BulkModulusFromDensityInternalEnergy(const Real rho,
                                                                 const Real sie,
                                                                 Real *lambda) const {
  const Real x = _rho0 / rho;
  // return
  // (_w+1)*(PressureFromDensityInternalEnergy(rho,sie)-ReferencePressure(rho))+x*(_A*_R1*std::exp(-_R1*x)+_B*_R2*std::exp(-_R2*x));
  return (_w + 1) * _w * rho * (sie - ReferenceEnergy(rho)) +
         x * (_A * _R1 * std::exp(-_R1 * x) + _B * _R2 * std::exp(-_R2 * x));
}
PORTABLE_FUNCTION
Real JWL::GruneisenParamFromDensityInternalEnergy(const Real rho, const Real sie,
                                                  Real *lambda) const {
  return _w;
}
// Below are "unimplemented" routines
PORTABLE_FUNCTION Real JWL::PressureFromDensityTemperature(const Real rho,
                                                           const Real temp,
                                                           Real *lambda) const {
  return PressureFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_FUNCTION Real JWL::SpecificHeatFromDensityTemperature(const Real rho,
                                                               const Real temp,
                                                               Real *lambda) const {
  return SpecificHeatFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_FUNCTION Real JWL::BulkModulusFromDensityTemperature(const Real rho,
                                                              const Real temp,
                                                              Real *lambda) const {
  return BulkModulusFromDensityInternalEnergy(
      rho, InternalEnergyFromDensityTemperature(rho, temp));
}
PORTABLE_FUNCTION
Real JWL::GruneisenParamFromDensityTemperature(const Real rho, const Real temp,
                                               Real *lambda) const {
  return _w;
}
PORTABLE_FUNCTION void JWL::DensityEnergyFromPressureTemperature(const Real press,
                                                                 const Real temp,
                                                                 Real *lambda, Real &rho,
                                                                 Real &sie) const {
  // sie = sie_r + cv*T;  Thus sie-sie_r = cv*T
  // Thus P = P_r +_w*rho*cv*T ==> Invertable?
  // Turns out not to be exactly invertible
  Real rhoguess = (rho < 1e-8) ? _rho0 : rho;
  auto PofRatT = [&](const Real r) { return _Cv * temp * r * _w + ReferencePressure(r); };
  using RootFinding1D::regula_falsi;
  using RootFinding1D::Status;
  RootFinding1D::RootCounts counts;
  auto status =
      regula_falsi(PofRatT, press, rhoguess, 1.0e-5, 1.0e3, 1.0e-8, 1.0e-8, rho, counts);
  if (status != Status::SUCCESS) {
    // Root finder failed even though the solution was bracketed... this is an error
    EOS_ERROR("JWL::DensityEnergyFromPressureTemperature: "
              "Root find failed to find a solution given P, T\n");
  }
  sie = InternalEnergyFromDensityTemperature(rho, temp);
}
PORTABLE_FUNCTION void JWL::FillEos(Real &rho, Real &temp, Real &sie, Real &press,
                                    Real &cv, Real &bmod, const unsigned long output,
                                    Real *lambda) const {
  if (output & thermalqs::pressure) press = PressureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::temperature)
    temp = TemperatureFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::bulk_modulus)
    bmod = BulkModulusFromDensityInternalEnergy(rho, sie);
  if (output & thermalqs::specific_heat)
    cv = SpecificHeatFromDensityInternalEnergy(rho, sie);
}

// TODO(JMM): pre-cache these rather than recomputing them each time
// TODO: Chad, please decide if STP is actually right here. Should it be
// based on the reference energy and pressure instead?
PORTABLE_FUNCTION
void JWL::ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                                 Real &bmod, Real &dpde, Real &dvdt, Real *lambda) const {
  rho = _rho0;
  temp = ROOM_TEMPERATURE;
  sie = InternalEnergyFromDensityTemperature(rho, temp, lambda);
  press = ATMOSPHERIC_PRESSURE;
  cv = _Cv;
  bmod = BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
  dpde = _w * _rho0;
  // TODO: chad please fix this one for me. This is wrong.
  Real gm1 = GruneisenParamFromDensityInternalEnergy(rho, sie, lambda) * rho;
  dvdt = gm1 * cv / bmod;
}

} // namespace singularity
