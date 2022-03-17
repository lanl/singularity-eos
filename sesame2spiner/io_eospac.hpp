//======================================================================
// sesame2spiner tool for converting eospac to spiner
// Author: Jonah Miller (jonahm@lanl.gov)
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
//======================================================================

#ifndef _SESAME2SPINER_IO_EOSPAC_HPP_
#define _SESAME2SPINER_IO_EOSPAC_HPP_

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <eos_Interface.h> // eospac API

#ifndef SPINER_USE_HDF
#error "HDF5 must be enabled"
#endif

#include <fast-math/logs.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/ports-of-call/portability.hpp>

#include <eospac-wrapper/eospac_wrapper.hpp>

using EospacWrapper::Verbosity;
using Spiner::DataBox;
using Spiner::RegularGrid1D;

// For logarithmic interpolation, quantities may be negative.
// If they are, use offset to ensure negative values make sense.
class Bounds {
 public:
  Bounds() {}

  Bounds(Real min, Real max, int N, Real offset)
      : grid(RegularGrid1D(min, max, N)), offset(offset) {}

  Bounds(Real min, Real max, int N, bool convertToLog = false, Real shrinkRange = 0,
         Real anchor_point = std::numeric_limits<Real>::signaling_NaN())
      : offset(0) {
    if (convertToLog) {
      // Log scales can't handle negative numbers or exactly zero. To
      // deal with that, we offset.
      constexpr Real epsilon = std::numeric_limits<float>::epsilon();
      const Real min_offset = 10 * std::abs(epsilon);
      // 1.1 so that the y-intercept isn't identically zero
      // when min < 0.
      // min_offset to handle the case where min=0
      if (min <= 0) offset = 1.1 * std::abs(min) + min_offset;

      min += offset;
      max += offset;

      min = std::log10(std::abs(min));
      max = std::log10(std::abs(max));
      Real delta = max - min;
      min += 0.5 * shrinkRange * delta;
      max -= 0.5 * shrinkRange * delta;

      if (!(std::isnan(anchor_point))) {
        anchor_point += offset;
        anchor_point = std::log10(std::abs(anchor_point));
      }
    }

    if (!(std::isnan(anchor_point))) {
      if (min < anchor_point && anchor_point < max) {
        Real dxguess = (max - min) / ((Real)N - 1);
        int Nmax = static_cast<int>((max - min) / dxguess);
        int Nanchor = static_cast<int>((anchor_point - min) / dxguess);
        Real dx = (anchor_point - min) / static_cast<Real>(Nanchor + 1);
        int Nmax_new = static_cast<int>((max - min) / dx);
        Nmax_new = std::max(Nmax, Nmax_new);
        max = dx * (Nmax_new) + min;
        N = Nmax_new;
      }
    }

    grid = RegularGrid1D(min, max, N);
  }

  inline Real log2lin(Real xl) const { return pow(10., xl) - offset; }
  inline Real i2lin(int i) const { return log2lin(grid.x(i)); }

  friend std::ostream &operator<<(std::ostream &os, const Bounds &b) {
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

void eosDataOfRhoSie(int matid, const Bounds &lRhoBounds, const Bounds &leBounds,
                     DataBox &P, DataBox &T, DataBox &bMods, DataBox &dPdRho,
                     DataBox &dPdE, DataBox &dTdRho, DataBox &dTdE, DataBox &dEdRho,
                     DataBox &mask, Verbosity eospacWarn = Verbosity::Quiet);

void eosDataOfRhoT(int matid, const Bounds &lRhoBounds, const Bounds &lTBounds,
                   DataBox &Ps, DataBox &sies, DataBox &bMods, DataBox &dPdRho,
                   DataBox &dPdE, DataBox &dTdRho, DataBox &dTdE, DataBox &dEdRho,
                   DataBox &dEdT, DataBox &mask, Verbosity eospacWarn = Verbosity::Quiet);

void eosColdCurves(int matid, const Bounds &lRhoBounds, DataBox &Ps, DataBox &sies,
                   DataBox &dPdRho, DataBox &dEdRho, DataBox &bMod, DataBox &mask,
                   Verbosity eospacWarn = Verbosity::Quiet);

void eosColdCurveMask(int matid, const Bounds &lRhoBounds, const int numSie,
                      const DataBox &sieColdCurve, DataBox &mask,
                      Verbosity eospacWarn = Verbosity::Quiet);

void makeInterpPoints(std::vector<EOS_REAL> &v, const Bounds &b);

#endif // _SESAME2SPINER_IO_EOSPAC_HPP_
