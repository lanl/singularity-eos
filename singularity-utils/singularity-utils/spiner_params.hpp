// © 2026. Triad National Security, LLC. All rights reserved.  This
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

#ifndef SINGULARITY_UTILS_SPINER_PARAMS_HPP_
#define SINGULARITY_UTILS_SPINER_PARAMS_HPP_

#include <limits>
#include <ports-of-call/portability.hpp>

#include <singularity-utils/bounds.hpp>

namespace singularity {

namespace table_utils {

static constexpr int NGRIDS = 3;
static constexpr Real TAnchor = 298.15; // K, 25°C - for piecewise grid construction
// Grid parameters for constructing Spiner tables from generic EOS
// Defaults match sesame2spiner behavior
// Used by both SpinerEOSDependsRhoT and SpinerEOSDependsRhoSie constructors
struct SpinerTableGridParams {
  // Density bounds
  Real rhoMin, rhoMax;
  int numRho = -1; // -1 means use numRhoPerDecade
  int numRhoPerDecade = 350;
  Real shrinklRhoBounds = 0.0;

  // Temperature bounds
  Real TMin, TMax;
  int numT = -1;
  int numTPerDecade = 100;
  Real shrinklTBounds = 0.0;

  // SIE bounds (energy can be negative!)
  // Note: Only used by SpinerEOSDependsRhoSie, ignored by SpinerEOSDependsRhoT
  Real sieMin, sieMax;
  int numSie = -1;
  int numSiePerDecade = 100;
  Real shrinkleBounds = 0.0;

  // Offset control (usually automatic, but allow override)
  // Set to -1 for auto-compute (default behavior)
  Real rhoOffset = -1.0;
  Real TOffset = -1.0;
  Real sieOffset = -1.0; // Only used by SpinerEOSDependsRhoSie

  // Enforce positive minimums (like sesame2spiner does for rho/T)
  // Set to <= 0 to disable enforcement
  Real strictlyPositiveMinRho = 1e-8;
  Real strictlyPositiveMinT = 1e-2;
  Real strictlyPositiveMinSie = -1.0; // disabled for sie (can be negative)

  // Material properties
  int matid = 0;
  Real Abar = std::numeric_limits<Real>::signaling_NaN();
  Real Zbar = std::numeric_limits<Real>::signaling_NaN();
  Real rhoNormal = std::numeric_limits<Real>::signaling_NaN();

  // Piecewise grid options (advanced - follow sesame2spiner defaults)
  bool piecewiseRho = true;
  bool piecewiseT = true;
  bool piecewiseSie = true; // Only used by SpinerEOSDependsRhoSie
  Real rhoCoarseFactorLo = 3.0;
  Real rhoCoarseFactorHi = 5.0;
  Real TCoarseFactor = 1.5;
  Real sieCoarseFactor = 1.5; // Only used by SpinerEOSDependsRhoSie
  Real rhoFineDiameterDecades = 1.5;
  Real TSplitPoint = 1e4;

  // Optional: fine grid bounds override (advanced use)
  Real rhoFineMin = -1.0; // -1 means use diameter
  Real rhoFineMax = -1.0;
};

// Construct density grid from parameters
inline void constructRhoBounds(const SpinerTableGridParams &params,
                               Bounds<NGRIDS> &lRhoBounds) {
  // Apply positive minimum if requested
  Real rhoMin = params.rhoMin;
  if (params.strictlyPositiveMinRho > 0.0) {
    rhoMin = std::max(rhoMin, params.strictlyPositiveMinRho);
  }
  Real rhoMax = params.rhoMax;

  int numRho;
  if (params.numRho > 0) {
    numRho = params.numRho;
  } else {
    numRho = Bounds<NGRIDS>::getNumPointsFromPPD(rhoMin, rhoMax, params.numRhoPerDecade);
  }

  // Compute rhoAnchor for piecewise grid center point
  // Match sesame2spiner behavior: use rhoNormal if valid, else geometric mean
  Real rhoAnchor;
  if (!std::isnan(params.rhoNormal) && params.rhoNormal > 0.0 &&
      params.rhoNormal <= 1e8) {
    rhoAnchor = params.rhoNormal;
  } else {
    rhoAnchor = std::sqrt(rhoMin * rhoMax);
  }

  // Construct piecewise or uniform grid
  if (params.piecewiseRho) {
    if (params.rhoFineMin > 0 && params.rhoFineMax > 0) {
      lRhoBounds = Bounds(Bounds<NGRIDS>::ThreeGrids(), rhoMin, rhoMax, rhoAnchor,
                          params.rhoFineMin, params.rhoFineMax, params.numRhoPerDecade,
                          params.rhoCoarseFactorLo, params.rhoCoarseFactorHi, true,
                          params.shrinklRhoBounds);
    } else {
      lRhoBounds = Bounds<NGRIDS>(Bounds<NGRIDS>::ThreeGrids(), rhoMin, rhoMax, rhoAnchor,
                                  params.rhoFineDiameterDecades, params.numRhoPerDecade,
                                  params.rhoCoarseFactorLo, params.rhoCoarseFactorHi,
                                  true, params.shrinklRhoBounds);
    }
  } else {
    lRhoBounds =
        Bounds<NGRIDS>(rhoMin, rhoMax, numRho, true, params.shrinklRhoBounds, rhoAnchor);
  }

  // Override offset if user specified
  if (params.rhoOffset >= 0) {
    lRhoBounds.offset = params.rhoOffset;
  }
}

// Construct temperature grid from parameters
inline void constructTBounds(const SpinerTableGridParams &params,
                             Bounds<NGRIDS> &lTBounds) {
  // Apply positive minimum if requested
  Real TMin = params.TMin;
  if (params.strictlyPositiveMinT > 0.0) {
    TMin = std::max(TMin, params.strictlyPositiveMinT);
  }
  Real TMax = params.TMax;

  int numT;
  if (params.numT > 0) {
    numT = params.numT;
  } else {
    numT = Bounds<NGRIDS>::getNumPointsFromPPD(TMin, TMax, params.numTPerDecade);
  }

  // Construct piecewise or uniform grid
  if (params.piecewiseT) {
    lTBounds = Bounds<NGRIDS>(Bounds<NGRIDS>::TwoGrids(), TMin, TMax, TAnchor,
                              params.TSplitPoint, params.numTPerDecade,
                              params.TCoarseFactor, true, params.shrinklTBounds);
  } else {
    lTBounds = Bounds<NGRIDS>(TMin, TMax, numT, true, params.shrinklTBounds, TAnchor);
  }

  // Override offset if user specified
  if (params.TOffset >= 0) {
    lTBounds.offset = params.TOffset;
  }
}

// Construct sie (specific internal energy) grid from parameters
// Only used by SpinerEOSDependsRhoSie
// Note: This function requires the source_eos to compute sie at anchor points
template <typename EOS>
inline void constructSieBounds(const EOS &source_eos, const SpinerTableGridParams &params,
                               Real rhoAnchor, Bounds<NGRIDS> &lSieBounds) {
  Real sieMin = params.sieMin;
  Real sieMax = params.sieMax;

  int numSie;
  if (params.numSie > 0) {
    numSie = params.numSie;
  } else {
    numSie = Bounds<NGRIDS>::getNumPointsFromPPD(sieMin, sieMax, params.numSiePerDecade);
  }

  // Construct piecewise or uniform grid
  if (params.piecewiseSie) {
    // Compute sie at anchor points for split
    Real sieAnchor =
        source_eos.InternalEnergyFromDensityTemperature(rhoAnchor, table_utils::TAnchor);
    Real sieSplitPoint =
        source_eos.InternalEnergyFromDensityTemperature(rhoAnchor, params.TSplitPoint);
    lSieBounds = Bounds<NGRIDS>(Bounds<NGRIDS>::TwoGrids(), sieMin, sieMax, sieAnchor,
                                sieSplitPoint, params.numSiePerDecade,
                                params.sieCoarseFactor, true, params.shrinkleBounds);
  } else {
    lSieBounds = Bounds<NGRIDS>(sieMin, sieMax, numSie, true, params.shrinkleBounds);
  }

  // Override offset if user specified
  if (params.sieOffset >= 0) {
    lSieBounds.offset = params.sieOffset;
  }
}

} // namespace table_utils
} // namespace singularity

#endif // SINGULARITY_UTILS_SPINER_PARAMS_HPP_
