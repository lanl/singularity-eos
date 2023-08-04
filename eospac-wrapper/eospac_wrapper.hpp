//======================================================================
// sesame2spiner tool for converting eospac to spiner
// Author: Jonah Miller (jonahm@lanl.gov)
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
//======================================================================

#ifndef _EOSPAC_WRAPPER_EOSPAC_WRAPPER_HPP_
#define _EOSPAC_WRAPPER_EOSPAC_WRAPPER_HPP_

#include <iostream>
#include <string>

#include <eos_Interface.h> // eospac API

namespace EospacWrapper {
enum class eospacSplit { none = 0, splitNumProp = 1, splitIdealGas = 2, splitCowan = 3 };
enum class eospacMonotonicity {
  none = 0,
  monotonicityX = 1,
  monotonicityXY = 2,
  monotonicityY = 3
};
inline double densityToSesame(const double CodeTemp) { return CodeTemp; }
inline double densityFromSesame(const double sesTemp) { return sesTemp; }
inline double temperatureToSesame(const double CodeTemp) { return CodeTemp; }
inline double temperatureFromSesame(const double SesTemp) { return SesTemp; }
inline double pressureFromSesame(const double SesPress) { return 1e10 * SesPress; }
inline double pressureToSesame(const double CodePress) { return 1e-10 * CodePress; }
inline double sieToSesame(const double CodeSie) { return 1e-10 * CodeSie; }
inline double sieFromSesame(const double SesSie) { return 1e10 * SesSie; }
inline double cvFromSesame(const double SesCv) { return 1e10 * SesCv; }
inline double bulkModulusFromSesame(const double SesBmod) { return 1e10 * SesBmod; }
inline double getBulkModulus(const double rho, const double P, const double DPDR_T,
                             const double DPDE_R, const double DEDR_T) {
  return rho * DPDR_T + DPDE_R * ((P / rho) - rho * DEDR_T);
}

enum class Verbosity { Quiet, Verbose, Debug };

class SesameMetadata {
 public:
  int matid;
  double exchangeCoefficient;
  double meanAtomicMass;
  double meanAtomicNumber;
  double solidBulkModulus;
  double normalDensity;
  double rhoMin, rhoMax;
  double TMin, TMax;
  double sieMin, sieMax;
  double rhoConversionFactor;
  double TConversionFactor;
  double sieConversionFactor;
  int numRho, numT;
  std::string comments;
  std::string name;
  friend std::ostream &operator<<(std::ostream &os, const SesameMetadata &m) {
    os << "MATID: " << m.matid << "\n"
       << "\tname: " << m.name << "\n"
       << "\texchange coefficient  = " << m.exchangeCoefficient << "\n"
       << "\tmean atomic mass      = " << m.meanAtomicMass << "\n"
       << "\tsolid bulk modulus    = " << m.solidBulkModulus << "\n"
       << "\tnormal density        = " << m.normalDensity << "\n"
       << "\t[rho min, rho max]    = "
       << "[" << m.rhoMin << ", " << m.rhoMax << "]"
       << "\n"
       << "\t[T min, T max]        = "
       << "[" << m.TMin << ", " << m.TMax << "]"
       << "\n"
       << "\t[sie min, sie max]    = "
       << "[" << m.sieMin << ", " << m.sieMax << "]"
       << "\n"
       << "\tnum rho               = " << m.numRho << "\n"
       << "\tnum T                 = " << m.numT << "\n"
       << "\tComments:\n"
       << m.comments << "\n";
    return os;
  }
};

void eosGetMetadata(int matid, SesameMetadata &metadata,
                    Verbosity eospacWarn = Verbosity::Quiet);

EOS_INTEGER eosSafeLoad(int ntables, int matid, EOS_INTEGER tableType[],
                        EOS_INTEGER tableHandle[], Verbosity eospacWarn,
                        bool invert_at_setup = false, EOS_REAL insert_data = 0.0,
                        eospacMonotonicity monotonicity = eospacMonotonicity::none,
                        bool apply_smoothing = false,
                        eospacSplit apply_splitting = eospacSplit::none,
                        bool linear_interp = false);
EOS_INTEGER eosSafeLoad(int ntables, int matid, EOS_INTEGER tableType[],
                        EOS_INTEGER tableHandle[],
                        const std::vector<std::string> &table_names, Verbosity eospacWarn,
                        bool invert_at_setup = false, EOS_REAL insert_data = 0.0,
                        eospacMonotonicity monotonicity = eospacMonotonicity::none,
                        bool apply_smoothing = false,
                        eospacSplit apply_splitting = eospacSplit::none,
                        bool linear_interp = false);

// output is boolean mask. 1 for no errors. 0 for errors.
bool eosSafeInterpolate(EOS_INTEGER *table, EOS_INTEGER nxypairs, EOS_REAL xVals[],
                        EOS_REAL yVals[], EOS_REAL var[], EOS_REAL dx[], EOS_REAL dy[],
                        const char tablename[], Verbosity eospacWarn,
                        EOS_INTEGER options[] = nullptr,
                        EOS_REAL option_values[] = nullptr, EOS_INTEGER nopts = 0);

void eosSafeTableInfo(EOS_INTEGER *table, EOS_INTEGER numInfoItems,
                      EOS_INTEGER infoItems[], EOS_REAL infoVals[], Verbosity eospacWarn);

void eosSafeTableCmnts(EOS_INTEGER *table, EOS_CHAR *comments, Verbosity eospacWarn);

void eosCheckError(EOS_INTEGER errorCode, const std::string &name, Verbosity eospacWarn);
std::string eosErrorString(EOS_INTEGER errorCode);
void eosSafeDestroy(int ntables, EOS_INTEGER tableHandles[], Verbosity eospacWarn);
std::string getName(std::string comment);

} // namespace EospacWrapper

#endif // _EOSPAC_WRAPPER_EOSPAC_WRAPPER_HPP_
