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

#ifndef _SESAME2SPINER_IO_EOSPAC_HPP_
#define _SESAME2SPINER_IO_EOSPAC_HPP_

#include <cmath>
#include <string>
#include <iostream>
#include <limits>
#include <vector>

#include <eos_Interface.h> // eospac API

#ifndef SPINER_USE_HDF
#error "HDF5 must be enabled"
#endif

#include <spiner/ports-of-call/portability.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <fast-math/logs.hpp>

using Spiner::DataBox;
using Spiner::RegularGrid1D;

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

// For logarithmic interpolation, quantities may be negative.
// If they are, use offset to ensure negative values make sense.
class Bounds {
public:
  Bounds() {}
  
  Bounds(Real min, Real max, int N, Real offset)
    : grid(RegularGrid1D(min,max,N))
    , offset(offset)
  {}

  Bounds(Real min, Real max, int N,
         bool convertToLog = false,
         Real shrinkRange = 0,
         Real anchor_point = std::numeric_limits<Real>::signaling_NaN())
    : offset(0) {    
    if (convertToLog) {
      // should be single-precision epsilon b/c that's whats used for the logs
      constexpr Real epsilon = std::numeric_limits<float>::epsilon();
      const Real min_offset = 10*std::abs(epsilon);
      if (min <= 0) offset = 1.1*std::abs(min) + min_offset;

      min += offset;
      max += offset;

      min = std::log10(std::abs(min));
      max = std::log10(std::abs(max));
      Real delta = max - min;
      min += 0.5*shrinkRange*delta;
      max -= 0.5*shrinkRange*delta;
      
      if (!(std::isnan(anchor_point))) {
        anchor_point += offset;
        anchor_point = std::log10(std::abs(anchor_point));
      }
    }

    if (!(std::isnan(anchor_point))) {
      if (min < anchor_point && anchor_point < max) {
        Real dxguess = (max-min)/((Real)N-1);
        int Nmax = static_cast<int>((max - min)/dxguess);
        int Nanchor = static_cast<int>((anchor_point - min)/dxguess);
        Real dx = (anchor_point - min)/static_cast<Real>(Nanchor+1);
        max = dx*(Nmax) + min;
        N += 1;
      }
    }
    
    grid = RegularGrid1D(min,max,N);
  }

  inline Real log2lin(Real xl) const {
    return pow(10.,xl) - offset;
  }
  inline Real i2lin(int i) const {
    return log2lin(grid.x(i));
  }

  friend std::ostream& operator<< (std::ostream& os, const Bounds& b) {
    os << "Bounds: [" << b.grid.min() << ", " << b.grid.max() << "]"
       << " + " << b.offset << ", "
       << "[N,dx] = [" << b.grid.nPoints() << ", " << b.grid.dx() << "]"
       << "\n";
    return os;
  }

public:
  RegularGrid1D grid;
  Real offset;
};

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

void eosDataOfRhoSie(int matid,
                     const Bounds& lRhoBounds,
                     const Bounds& leBounds,
                     DataBox& P,
                     DataBox& T,
                     DataBox& bMods,
                     DataBox& dPdRho,
                     DataBox& dPdE,
                     DataBox& dTdRho,
                     DataBox& dTdE,
                     DataBox& dEdRho,
                     DataBox& mask,
                     Verbosity eospacWarn = Verbosity::Quiet
                     );

void eosDataOfRhoT(int matid,
                   const Bounds& lRhoBounds,
                   const Bounds& lTBounds,
                   DataBox& Ps,
                   DataBox& sies,
                   DataBox& bMods,
                   DataBox& dPdRho,
                   DataBox& dPdE,
                   DataBox& dTdRho,
                   DataBox& dTdE,
                   DataBox& dEdRho,
                   DataBox& dEdT,
                   DataBox& mask,
                   Verbosity eospacWarn = Verbosity::Quiet
                   );

void eosColdCurves(int matid,
                   const Bounds& lRhoBounds,
                   DataBox& Ps,
                   DataBox& sies,
                   DataBox& dPdRho,
                   DataBox& dEdRho,
                   DataBox& bMod,
                   DataBox& mask,
                   Verbosity eospacWarn = Verbosity::Quiet
                   );

void eosColdCurveMask(int matid,
                      const Bounds& lRhoBounds,
                      const int numSie,
                      const DataBox& sieColdCurve,
                      DataBox& mask,
                      Verbosity eospacWarn = Verbosity::Quiet
                      );

EOS_INTEGER eosSafeLoad(int ntables, int matid,
                        EOS_INTEGER tableType[],
                        EOS_INTEGER tableHandle[],
                        Verbosity eospacWarn);
EOS_INTEGER eosSafeLoad(int ntables, int matid,
                        EOS_INTEGER tableType[],
                        EOS_INTEGER tableHandle[],
                        const std::vector<std::string>& table_names,
                        Verbosity eospacWarn);

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

void makeInterpPoints(std::vector<EOS_REAL>& v, const Bounds& b);
std::string getName(std::string comment);

#endif // _SESAME2SPINER_IO_EOSPAC_HPP_
