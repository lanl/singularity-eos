//======================================================================
// sesame2spiner tool for converting eospac to spiner
// Author: Jonah Miller (jonahm@lanl.gov)
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
//======================================================================

#ifndef _EOSPAC_WRAPPER_EOSPAC_WRAPPER_HPP_
#define _EOSPAC_WRAPPER_EOSPAC_WRAPPER_HPP_

#include <string>
#include <iostream>

#include <eos_Interface.h> // eospac API

namespace EospacWrapper {

inline Real densityToSesame(const Real CodeTemp) { return CodeTemp;}
inline Real densityFromSesame(const Real sesTemp) { return sesTemp;}
inline Real temperatureToSesame(const Real CodeTemp) { return CodeTemp;}
inline Real temperatureFromSesame(const Real SesTemp) { return SesTemp;}
inline Real pressureFromSesame(const Real SesPress) { return 1e10*SesPress;}
inline Real pressureToSesame(const Real CodePress) { return 1e-10*CodePress;}
inline Real sieToSesame(const Real CodeSie) { return 1e-10*CodeSie;}
inline Real sieFromSesame(const Real SesSie) { return 1e10*SesSie;}
inline Real cvFromSesame(const Real SesCv) { return 1e10*SesCv;}
inline Real bulkModulusFromSesame(const Real SesBmod) { return 1e10*SesBmod;}
inline Real getBulkModulus(const Real rho, const Real P,
                           const Real DPDR_T, const Real DPDE_R,
                           const Real DEDR_T) {
  return rho*DPDR_T + DPDE_R*((P/rho) - rho*DEDR_T);
}

enum class Verbosity { Quiet, Verbose, Debug };

class SesameMetadata {
public:
  int matid;
  Real exchangeCoefficient;
  Real meanAtomicMass;
  Real meanAtomicNumber;
  Real solidBulkModulus;
  Real normalDensity;
  Real rhoMin, rhoMax;
  Real TMin, TMax;
  Real sieMin, sieMax;
  Real rhoConversionFactor;
  Real TConversionFactor;
  Real sieConversionFactor;
  int numRho, numT;
  std::string comments;
  std::string name;
  friend std::ostream& operator<< (std::ostream& os, const SesameMetadata& m) {
    os << "MATID: "  << m.matid << "\n"
       << "\tname: " << m.name << "\n"
       << "\texchange coefficient  = " << m.exchangeCoefficient << "\n"
       << "\tmean atomic mass      = " << m.meanAtomicMass      << "\n"
       << "\tsolid bulk modulus    = " << m.solidBulkModulus    << "\n"
       << "\tnormal density        = " << m.normalDensity       << "\n"
       << "\t[rho min, rho max]    = "
       << "[" << m.rhoMin << ", " << m.rhoMax << "]" << "\n"
       << "\t[T min, T max]        = "
       << "[" << m.TMin << ", " << m.TMax << "]" << "\n"
       << "\t[sie min, sie max]    = "
       << "[" << m.sieMin << ", " << m.sieMax << "]" << "\n"
       << "\tnum rho               = " << m.numRho  << "\n"
       << "\tnum T                 = " << m.numT    << "\n"
       << "\tComments:\n"
       << m.comments << "\n";
    return os;
  }
};

void eosGetMetadata(int matid, SesameMetadata& metadata,
                    Verbosity eospacWarn = Verbosity::Quiet);

EOS_INTEGER eosSafeLoad(int ntables, int matid,
                        EOS_INTEGER tableType[],
                        EOS_INTEGER tableHandle[],
                        Verbosity eospacWarn,
			bool invert_at_setup=false);
EOS_INTEGER eosSafeLoad(int ntables, int matid,
                        EOS_INTEGER tableType[],
                        EOS_INTEGER tableHandle[],
                        const std::vector<std::string>& table_names,
                        Verbosity eospacWarn,
			bool invert_at_setup=false);

// output is boolean mask. 1 for no errors. 0 for errors.
bool eosSafeInterpolate(EOS_INTEGER *table,
			EOS_INTEGER nxypairs,
			EOS_REAL xVals[],
			EOS_REAL yVals[],
			EOS_REAL var[],
			EOS_REAL dx[],
			EOS_REAL dy[],
			const char tablename[],
			Verbosity eospacWarn);

void eosSafeTableInfo(EOS_INTEGER* table,
                      EOS_INTEGER numInfoItems,
                      EOS_INTEGER infoItems[],
                      EOS_REAL infoVals[],
                      Verbosity eospacWarn);

void eosSafeTableCmnts(EOS_INTEGER* table, EOS_CHAR* comments,
                       Verbosity eospacWarn);

void eosCheckError(EOS_INTEGER errorCode, const std::string& name,
                   Verbosity eospacWarn);
std::string eosErrorString(EOS_INTEGER errorCode);
void eosSafeDestroy(int ntables, EOS_INTEGER tableHandles[],
                    Verbosity eospacWarn);
std::string getName(std::string comment);

} // namespace EospacWrapper

#endif // _EOSPAC_WRAPPER_EOSPAC_WRAPPER_HPP_
