//------------------------------------------------------------------------------
// © 2021-2026. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_EOS_SPINER_CONSTRUCTION_HPP_
#define _SINGULARITY_EOS_EOS_SPINER_CONSTRUCTION_HPP_

#ifdef SINGULARITY_USE_SPINER_WITH_HDF5

#include <limits>

// ports-of-call
#include <ports-of-call/portability.hpp>

// spiner
#include <spiner/databox.hpp>
#include <spiner/spiner_types.hpp>

// singularity-eos
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/eos_concepts.hpp>
#include <singularity-eos/base/fast-math/logs.hpp>
#include <singularity-eos/base/finite_diff.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>
#include <singularity-eos/base/spiner_table_utils.hpp>

namespace singularity {
namespace spiner_table_builder {

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

// Shared constants
static constexpr int NGRIDS = 3;
static constexpr Real TAnchor = 298.15; // K, 25°C - for piecewise grid construction
using Bounds = table_utils::Bounds<NGRIDS>;
using Grid_t = Spiner::PiecewiseGrid1D<Real, NGRIDS>;
using DataBox = Spiner::DataBox<Real, Grid_t>;

// Helper functions for log transformations
PORTABLE_FORCEINLINE_FUNCTION Real to_log(const Real x, const Real offset) {
  return FastMath::log10(std::abs(std::max(x, -offset) + offset) + robust::EPS());
}

PORTABLE_FORCEINLINE_FUNCTION Real from_log(const Real lx, const Real offset) {
  return FastMath::pow10(lx) - offset;
}

// Extract material properties from source EOS with fallback to params
template <typename EOS>
inline Real extractAbar(const EOS &eos, const SpinerTableGridParams &params) {
  if constexpr (eos_concepts::has_abar_v<EOS>) {
    return std::isnan(params.Abar) ? eos.MeanAtomicMass() : params.Abar;
  } else {
    return std::isnan(params.Abar) ? 1.0 : params.Abar;
  }
}

template <typename EOS>
inline Real extractZbar(const EOS &eos, const SpinerTableGridParams &params) {
  if constexpr (eos_concepts::has_zbar_v<EOS>) {
    return std::isnan(params.Zbar) ? eos.MeanAtomicNumber() : params.Zbar;
  } else {
    return std::isnan(params.Zbar) ? 1.0 : params.Zbar;
  }
}

// Construct density grid from parameters
inline void constructRhoBounds(const SpinerTableGridParams &params, Bounds &lRhoBounds) {
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
    numRho = Bounds::getNumPointsFromPPD(rhoMin, rhoMax, params.numRhoPerDecade);
  }

  Real rhoAnchor = std::sqrt(rhoMin * rhoMax);

  // Construct piecewise or uniform grid
  if (params.piecewiseRho) {
    if (params.rhoFineMin > 0 && params.rhoFineMax > 0) {
      lRhoBounds =
          Bounds(Bounds::ThreeGrids(), rhoMin, rhoMax, rhoAnchor, params.rhoFineMin,
                 params.rhoFineMax, params.numRhoPerDecade, params.rhoCoarseFactorLo,
                 params.rhoCoarseFactorHi, true, params.shrinklRhoBounds);
    } else {
      lRhoBounds = Bounds(Bounds::ThreeGrids(), rhoMin, rhoMax, rhoAnchor,
                          params.rhoFineDiameterDecades, params.numRhoPerDecade,
                          params.rhoCoarseFactorLo, params.rhoCoarseFactorHi, true,
                          params.shrinklRhoBounds);
    }
  } else {
    lRhoBounds = Bounds(rhoMin, rhoMax, numRho, true, params.shrinklRhoBounds, rhoAnchor);
  }

  // Override offset if user specified
  if (params.rhoOffset >= 0) {
    lRhoBounds.offset = params.rhoOffset;
  }
}

// Construct temperature grid from parameters
inline void constructTBounds(const SpinerTableGridParams &params, Bounds &lTBounds) {
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
    numT = Bounds::getNumPointsFromPPD(TMin, TMax, params.numTPerDecade);
  }

  // Construct piecewise or uniform grid
  if (params.piecewiseT) {
    lTBounds =
        Bounds(Bounds::TwoGrids(), TMin, TMax, TAnchor, params.TSplitPoint,
               params.numTPerDecade, params.TCoarseFactor, true, params.shrinklTBounds);
  } else {
    lTBounds = Bounds(TMin, TMax, numT, true, params.shrinklTBounds, TAnchor);
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
                               Real rhoAnchor, Bounds &lSieBounds) {
  Real sieMin = params.sieMin;
  Real sieMax = params.sieMax;

  int numSie;
  if (params.numSie > 0) {
    numSie = params.numSie;
  } else {
    numSie = Bounds::getNumPointsFromPPD(sieMin, sieMax, params.numSiePerDecade);
  }

  // Construct piecewise or uniform grid
  if (params.piecewiseSie) {
    // Compute sie at anchor points for split
    Real sieAnchor = source_eos.InternalEnergyFromDensityTemperature(rhoAnchor, TAnchor);
    Real sieSplitPoint =
        source_eos.InternalEnergyFromDensityTemperature(rhoAnchor, params.TSplitPoint);
    lSieBounds = Bounds(Bounds::TwoGrids(), sieMin, sieMax, sieAnchor, sieSplitPoint,
                        params.numSiePerDecade, params.sieCoarseFactor, true,
                        params.shrinkleBounds);
  } else {
    lSieBounds = Bounds(sieMin, sieMax, numSie, true, params.shrinkleBounds);
  }

  // Override offset if user specified
  if (params.sieOffset >= 0) {
    lSieBounds.offset = params.sieOffset;
  }
}

// Compute cold curves (properties at minimum temperature)
// Used by both RhoT and RhoSie constructors
template <typename EOS>
inline void computeColdCurves(const EOS &source_eos, const Bounds &lRhoBounds, Real Tmin,
                              Real lRhoOffset, DataBox &PCold, DataBox &sieCold,
                              DataBox &bModCold, DataBox &dPdRhoCold) {
  const int numRho = lRhoBounds.grid.nPoints();

  // Allocate cold curve arrays
  PCold.resize(numRho);
  PCold.setRange(0, lRhoBounds.grid);
  sieCold.resize(numRho);
  sieCold.setRange(0, lRhoBounds.grid);
  bModCold.resize(numRho);
  bModCold.setRange(0, lRhoBounds.grid);
  dPdRhoCold.resize(numRho);
  dPdRhoCold.setRange(0, lRhoBounds.grid);

  // Evaluate at minimum temperature
  for (int j = 0; j < numRho; ++j) {
    Real lRho = lRhoBounds.grid.x(j);
    Real rho = from_log(lRho, lRhoOffset);

    // Internal energy and pressure at Tmin
    if constexpr (eos_concepts::has_sie_rho_T_v<EOS>) {
      sieCold(j) = source_eos.InternalEnergyFromDensityTemperature(rho, Tmin);
    } else {
      sieCold(j) = 0.0; // Fallback
    }

    if constexpr (eos_concepts::has_P_rho_T_v<EOS>) {
      PCold(j) = source_eos.PressureFromDensityTemperature(rho, Tmin);
    } else {
      PCold(j) = 0.0; // Fallback
    }

    // Compute dP/drho at constant T via finite difference
    auto PofR = [&source_eos, Tmin](Real r) {
      return source_eos.PressureFromDensityTemperature(r, Tmin);
    };
    dPdRhoCold(j) = finite_diff::centralDifference(PofR, rho);

    // Bulk modulus - compute isentropic B_S, not isothermal B_T
    // For ideal gas: B_S = (Gamma + 1) * P where Gamma is Gruneisen parameter
    // Note: GruneisenParam returns Gamma = (gamma - 1), not heat capacity ratio gamma
    if constexpr (eos_concepts::has_gamma_rho_T_v<EOS>) {
      Real Gamma = source_eos.GruneisenParamFromDensityTemperature(rho, Tmin);
      Real P = PCold(j);
      bModCold(j) = (Gamma + 1.0) * P;
    } else {
      // Fallback: isothermal (not ideal, but better than nothing)
      bModCold(j) = rho * dPdRhoCold(j);
    }
  }
}

// Populate (rho, T) indexed tables from a generic EOS
// This function is shared between SpinerEOSDependsRhoT and SpinerEOSDependsRhoSie
// constructors to ensure consistent table population
//
// Each set* callable has signature: void(int j, int i, Real value)
template <typename EOS, typename SetSie, typename SetP, typename SetBMod,
          typename SetDPdRho, typename SetDPdE, typename SetDTdRho, typename SetDTdE,
          typename SetDEdRho, typename SetDEdT>
inline void populateDependsRhoT(const EOS &source_eos, const Bounds &lRhoBounds,
                                const Bounds &lTBounds, Real lRhoOffset, Real lTOffset,
                                Real minimumRho, Real maximumRho, Real minimumT,
                                Real maximumT, Real sieMin, Real sieMax, SetSie setSie,
                                SetP setP, SetBMod setBMod, SetDPdRho setDPdRho,
                                SetDPdE setDPdE, SetDTdRho setDTdRho, SetDTdE setDTdE,
                                SetDEdRho setDEdRho, SetDEdT setDEdT) {
  const int numRho = lRhoBounds.grid.nPoints();
  const int numT = lTBounds.grid.nPoints();

  for (int j = 0; j < numRho; j++) {
    Real lRho = lRhoBounds.grid.x(j);
    Real rho = from_log(lRho, lRhoOffset);

    for (int i = 0; i < numT; i++) {
      Real lT = lTBounds.grid.x(i);
      Real T = from_log(lT, lTOffset);

      // Fill sie field
      Real sie;
      if constexpr (eos_concepts::has_sie_rho_T_v<EOS>) {
        sie = source_eos.InternalEnergyFromDensityTemperature(rho, T);
      } else if constexpr (eos_concepts::has_T_rho_sie_v<EOS>) {
        // Invert T(rho, sie) to find sie given (rho, T)
        const auto T_from_sie = [&source_eos, rho](const Real sie_trial) -> Real {
          return source_eos.TemperatureFromDensityInternalEnergy(rho, sie_trial);
        };
        Real sie_lower = sieMin;
        Real sie_upper = sieMax;
        Real sie_guess = 0.5 * (sie_lower + sie_upper);
        Real sie_solution;
        RootFinding1D::Status status =
            SP_ROOT_FINDER(T_from_sie, T, sie_guess, sie_lower, sie_upper, robust::EPS(),
                           robust::EPS(), sie_solution, nullptr);
        PORTABLE_ALWAYS_REQUIRE(status == RootFinding1D::Status::SUCCESS,
                                "Failed to invert TemperatureFromDensityInternalEnergy "
                                "during table construction");
        sie = sie_solution;
      } else {
        static_assert(eos_concepts::dependent_false_v<EOS>,
                      "Source eos must provide either: \n"
                      "InternalEnergyFromDensityTemperature or "
                      "TemperatureFromDensityInternalEnergy.\n");
      }
      setSie(j, i, sie);

      // Fill pressure field
      Real P;
      if constexpr (eos_concepts::has_P_rho_T_v<EOS>) {
        P = source_eos.PressureFromDensityTemperature(rho, T);
      } else if (eos_concepts::has_P_rho_sie_v<EOS>) {
        P = source_eos.PressureFromDensityInternalEnergy(rho, sie);
      } else {
        static_assert(eos_concepts::dependent_false_v<EOS>,
                      "Source eos must provide either: \n"
                      "PressureFromDensityTemperature or "
                      "PressureFromDensityInternalEnergy.\n");
      }
      setP(j, i, P);

      // Fill dPdE (at constant rho)
      Real dPdE;
      if constexpr (eos_concepts::has_gamma_rho_T_v<EOS>) {
        Real gamma = source_eos.GruneisenParamFromDensityTemperature(rho, T);
        dPdE = rho * gamma;
      } else if constexpr (eos_concepts::has_gamma_rho_sie_v<EOS>) {
        Real gamma = source_eos.GruneisenParamFromDensityInternalEnergy(rho, sie);
        dPdE = rho * gamma;
      } else if constexpr (eos_concepts::has_P_rho_sie_v<EOS>) {
        auto PofE = [&source_eos, rho](double sie) {
          return source_eos.PressureFromDensityInternalEnergy(rho, sie);
        };
        dPdE = finite_diff::centralDifference(PofE, sie);
      } else if constexpr (eos_concepts::has_T_rho_sie_v<EOS>) {
        auto PofE = [&source_eos, rho](Real sie) {
          auto T = source_eos.TemperatureFromDensityInternalEnergy(rho, sie);
          return source_eos.PressureFromDensityTemperature(rho, T);
        };
        dPdE = finite_diff::centralDifference(PofE, sie);
      } else {
        auto PofT = [&source_eos, rho](Real T) {
          return source_eos.PressureFromDensityTemperature(rho, T);
        };
        Real dPdT = finite_diff::finiteDifference(PofT, T, minimumT, maximumT);
        auto EofT = [&source_eos, rho](Real T) {
          return source_eos.InternalEnergyFromDensityTemperature(rho, T);
        };
        Real dEdT = finite_diff::finiteDifference(EofT, T, minimumT, maximumT);
        dPdE = robust::ratio(dPdT, dEdT);
      }
      setDPdE(j, i, dPdE);

      // Fill dTdE (at constant rho)
      Real dTdE;
      if constexpr (eos_concepts::has_cv_rho_T_v<EOS>) {
        Real cv = source_eos.SpecificHeatFromDensityTemperature(rho, T);
        dTdE = robust::ratio(1.0, cv);
      } else if constexpr (eos_concepts::has_cv_rho_sie_v<EOS>) {
        Real cv = source_eos.SpecificHeatFromDensityInternalEnergy(rho, sie);
        dTdE = robust::ratio(1.0, cv);
      } else if constexpr (eos_concepts::has_T_rho_sie_v<EOS>) {
        auto TofE = [&source_eos, rho](Real sie) {
          return source_eos.TemperatureFromDensityInternalEnergy(rho, sie);
        };
        dTdE = finite_diff::centralDifference(TofE, sie);
      } else {
        auto EofT = [&source_eos, rho](Real T) {
          return source_eos.InternalEnergyFromDensityTemperature(rho, T);
        };
        Real dEdT_inv = finite_diff::finiteDifference(EofT, T, minimumT, maximumT);
        dTdE = robust::ratio(1.0, dEdT_inv);
      }
      setDTdE(j, i, dTdE);

      // Fill dEdT (at constant rho)
      Real dEdT;
      {
        auto EofT = [&source_eos, rho](Real T) {
          return source_eos.InternalEnergyFromDensityTemperature(rho, T);
        };
        dEdT = finite_diff::finiteDifference(EofT, T, minimumT, maximumT);
      }
      setDEdT(j, i, dEdT);

      // Fill dEdRho (at constant T)
      Real dEdRho;
      {
        auto EofR = [&source_eos, T](Real rho) {
          return source_eos.InternalEnergyFromDensityTemperature(rho, T);
        };
        dEdRho = finite_diff::finiteDifference(EofR, rho, minimumRho, maximumRho);
      }
      setDEdRho(j, i, dEdRho);

      // Fill dPdRho (at constant E, not constant T!)
      // Compute dPdRho_T first, then convert using chain rule
      Real dPdRho_E;
      {
        auto PofR = [&source_eos, T](Real rho) {
          return source_eos.PressureFromDensityTemperature(rho, T);
        };
        Real dPdRho_T = finite_diff::finiteDifference(PofR, rho, minimumRho, maximumRho);
        // Convert to constant E: (dP/drho)_E = (dP/drho)_T - (dP/dE)_rho * (dE/drho)_T
        dPdRho_E = dPdRho_T - dPdE * dEdRho;
      }
      setDPdRho(j, i, dPdRho_E);

      // Fill dTdRho (at constant E)
      // Use chain rule: dTdRho_E = -dEdRho_T * dTdE
      Real dTdRho_E = -dEdRho * dTdE;
      setDTdRho(j, i, dTdRho_E);

      // Bulk modulus - initialize to small value, will be computed by fixBulkModulus_()
      // or calcBMod_() in the calling constructor
      setBMod(j, i, robust::EPS());
    }
  }
}

} // namespace spiner_table_builder
} // namespace singularity

#endif // SINGULARITY_USE_SPINER_WITH_HDF5
#endif // _SINGULARITY_EOS_EOS_SPINER_CONSTRUCTION_HPP_
