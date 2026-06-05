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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/eos/eos.hpp>

// This file was generated in part with the assistance of generative AI

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <test/eos_unit_test_helpers.hpp>

#ifdef SINGULARITY_USE_SPINER_WITH_HDF5

using singularity::IdealGas;
using singularity::SpinerEOSDependsRhoSie;
using singularity::SpinerEOSDependsRhoT;
using singularity::SpinerTableGridParams;

// Minimal EOS classes for testing fallback paths
// These deliberately provide only the bare minimum required interface
// to force the constructor to use finite differences, root finding, etc.

// MinimalEOS1: Provides basic methods but NO gamma, cv, or bulk modulus
// Forces: all derivatives via finite differences, no optimizations
class MinimalEOS1 {
 private:
  IdealGas ideal_; // Use IdealGas internally for physics

 public:
  MinimalEOS1(Real gm1, Real Cv) : ideal_(gm1, Cv) {}

  // Provide core thermodynamic methods (required by constructor)
  PORTABLE_INLINE_FUNCTION
  Real TemperatureFromDensityInternalEnergy(Real rho, Real sie) const {
    return ideal_.TemperatureFromDensityInternalEnergy(rho, sie);
  }

  PORTABLE_INLINE_FUNCTION
  Real InternalEnergyFromDensityTemperature(Real rho, Real T) const {
    return ideal_.InternalEnergyFromDensityTemperature(rho, T);
  }

  PORTABLE_INLINE_FUNCTION
  Real PressureFromDensityTemperature(Real rho, Real T) const {
    return ideal_.PressureFromDensityTemperature(rho, T);
  }

  // Deliberately DO NOT provide:
  // - GruneisenParamFromDensityTemperature (forces finite diff for dPdE)
  // - SpecificHeatFromDensityTemperature (forces finite diff for dTdE)
  // - BulkModulusFromDensityTemperature (forces calcBMod_ computation)

  // Provide mean atomic properties
  PORTABLE_INLINE_FUNCTION
  Real MeanAtomicMass() const { return ideal_.MeanAtomicMass(); }

  PORTABLE_INLINE_FUNCTION
  Real MeanAtomicNumber() const { return ideal_.MeanAtomicNumber(); }
};

// MinimalEOS2: Provides E(rho,T) and P(rho,T) but NOT P(rho,sie)
// Forces: chain rule conversions for dependsRhoSie_ derivatives
class MinimalEOS2 {
 private:
  IdealGas ideal_;

 public:
  MinimalEOS2(Real gm1, Real Cv) : ideal_(gm1, Cv) {}

  PORTABLE_INLINE_FUNCTION
  Real InternalEnergyFromDensityTemperature(Real rho, Real T) const {
    return ideal_.InternalEnergyFromDensityTemperature(rho, T);
  }

  PORTABLE_INLINE_FUNCTION
  Real PressureFromDensityTemperature(Real rho, Real T) const {
    return ideal_.PressureFromDensityTemperature(rho, T);
  }

  PORTABLE_INLINE_FUNCTION
  Real TemperatureFromDensityInternalEnergy(Real rho, Real sie) const {
    return ideal_.TemperatureFromDensityInternalEnergy(rho, sie);
  }

  PORTABLE_INLINE_FUNCTION
  Real MeanAtomicMass() const { return ideal_.MeanAtomicMass(); }

  PORTABLE_INLINE_FUNCTION
  Real MeanAtomicNumber() const { return ideal_.MeanAtomicNumber(); }
};

SCENARIO("SpinerEOS construction from IdealGas", "[SpinerEOS][Constructor][IdealGas]") {
  GIVEN("An IdealGas EOS and grid parameters") {
    // Create IdealGas
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.4; // gamma - 1, so gamma = 1.4
    IdealGas ideal_eos(gm1, Cv);

    // Set up grid parameters
    SpinerTableGridParams params;
    params.rhoMin = 1e-3;
    params.rhoMax = 1e3;
    params.TMin = 1e2;
    params.TMax = 1e5;
    params.sieMin = Cv * params.TMin;
    params.sieMax = Cv * params.TMax;

    // Use coarser grid for testing (faster)
    params.numRhoPerDecade = 50;
    params.numTPerDecade = 50;
    params.numSiePerDecade = 50;

    // Disable piecewise grids for simpler testing
    params.piecewiseRho = false;
    params.piecewiseT = false;
    params.piecewiseSie = false;

    params.matid = 1001;

    WHEN("We construct a SpinerEOS from the IdealGas") {
      SpinerEOSDependsRhoSie spiner_eos(ideal_eos, params);

      THEN("The SpinerEOS should interpolate correctly at grid points") {
        // Test at several grid points
        Real rho = 1.0;
        Real T = 1000.0;

        Real P_ideal = ideal_eos.PressureFromDensityTemperature(rho, T);
        Real P_spiner = spiner_eos.PressureFromDensityTemperature(rho, T);

        Real sie_ideal = ideal_eos.InternalEnergyFromDensityTemperature(rho, T);
        Real sie_spiner = spiner_eos.InternalEnergyFromDensityTemperature(rho, T);

        INFO("IdealGas P: " << P_ideal << "  Spiner P: " << P_spiner);
        INFO("IdealGas sie: " << sie_ideal << "  Spiner sie: " << sie_spiner);

        // Should be very close (interpolation error)
        CHECK(isClose(P_spiner, P_ideal, 1e-2)); // 1% tolerance
        CHECK(isClose(sie_spiner, sie_ideal, 1e-2));
      }

      AND_THEN("The SpinerEOS should give correct temperature from sie") {
        Real rho = 10.0;
        Real sie = Cv * 5000.0; // corresponds to T = 5000

        Real T_ideal = ideal_eos.TemperatureFromDensityInternalEnergy(rho, sie);
        Real T_spiner = spiner_eos.TemperatureFromDensityInternalEnergy(rho, sie);

        INFO("IdealGas T: " << T_ideal << "  Spiner T: " << T_spiner);
        CHECK(isClose(T_spiner, T_ideal, 1e-2));
      }

      AND_THEN("Bulk modulus should be reasonable") {
        Real rho = 1.0;
        Real T = 1000.0;

        Real bmod_ideal = ideal_eos.BulkModulusFromDensityTemperature(rho, T);
        Real bmod_spiner = spiner_eos.BulkModulusFromDensityTemperature(rho, T);

        INFO("IdealGas bmod: " << bmod_ideal << "  Spiner bmod: " << bmod_spiner);
        CHECK(bmod_spiner > 0);                       // Should be positive
        CHECK(isClose(bmod_spiner, bmod_ideal, 0.1)); // 10% tolerance (derivatives)
      }

      AND_THEN("Metadata should be set correctly") {
        CHECK(spiner_eos.matid() == params.matid);
        CHECK(spiner_eos.MinimumDensity() > 0);
        CHECK(spiner_eos.MaximumDensity() > spiner_eos.MinimumDensity());
      }
    }
  }
}

SCENARIO("SpinerEOS construction with piecewise grids",
         "[SpinerEOS][Constructor][PiecewiseGrid]") {
  GIVEN("An IdealGas and piecewise grid parameters") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.4;
    IdealGas ideal_eos(gm1, Cv);

    SpinerTableGridParams params;
    params.rhoMin = 1e-3;
    params.rhoMax = 1e3;
    params.TMin = 1e2;
    params.TMax = 1e5;
    params.sieMin = Cv * params.TMin;
    params.sieMax = Cv * params.TMax;

    // Use piecewise grids (default)
    params.piecewiseRho = true;
    params.piecewiseT = true;
    params.piecewiseSie = true;

    // Coarser for faster testing
    params.numRhoPerDecade = 30;
    params.numTPerDecade = 30;
    params.numSiePerDecade = 30;

    params.matid = 1002;

    WHEN("We construct with piecewise grids") {
      SpinerEOSDependsRhoSie spiner_eos(ideal_eos, params);

      THEN("Interpolation should still work correctly") {
        Real rho = 1.0;
        Real T = 1000.0;

        Real P_ideal = ideal_eos.PressureFromDensityTemperature(rho, T);
        Real P_spiner = spiner_eos.PressureFromDensityTemperature(rho, T);

        INFO("Piecewise grid - IdealGas P: " << P_ideal << "  Spiner P: " << P_spiner);
        CHECK(isClose(P_spiner, P_ideal, 1e-2));
      }
    }
  }
}

SCENARIO("SpinerEOS accuracy test", "[SpinerEOS][Constructor][Accuracy]") {
  GIVEN("An IdealGas and fine grid") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.4;
    IdealGas ideal_eos(gm1, Cv);

    SpinerTableGridParams params;
    params.rhoMin = 1e-2;
    params.rhoMax = 1e2;
    params.TMin = 1e2;
    params.TMax = 1e4;
    params.sieMin = Cv * params.TMin;
    params.sieMax = Cv * params.TMax;

    // Fine grid for better accuracy
    params.numRhoPerDecade = 100;
    params.numTPerDecade = 100;
    params.numSiePerDecade = 100;
    params.piecewiseRho = false;
    params.piecewiseT = false;
    params.piecewiseSie = false;

    params.matid = 1003;

    WHEN("We construct with a fine grid") {
      SpinerEOSDependsRhoSie spiner_eos(ideal_eos, params);

      THEN("Interpolation should be very accurate") {
        // Test at random points within bounds
        constexpr int ntest = 10;
        for (int i = 0; i < ntest; i++) {
          // Random points (geometric spacing in log)
          Real log_rho = std::log10(params.rhoMin) +
                         (std::log10(params.rhoMax) - std::log10(params.rhoMin)) *
                             static_cast<Real>(i) / ntest;
          Real rho = std::pow(10.0, log_rho);

          Real log_T = std::log10(params.TMin) +
                       (std::log10(params.TMax) - std::log10(params.TMin)) *
                           static_cast<Real>(i) / ntest;
          Real T = std::pow(10.0, log_T);

          Real P_ideal = ideal_eos.PressureFromDensityTemperature(rho, T);
          Real P_spiner = spiner_eos.PressureFromDensityTemperature(rho, T);

          Real rel_error = std::abs(P_spiner - P_ideal) / std::abs(P_ideal);

          INFO("Test point " << i << ": rho=" << rho << " T=" << T);
          INFO("P_ideal=" << P_ideal << " P_spiner=" << P_spiner
                          << " rel_error=" << rel_error);

          CHECK(rel_error < 0.01); // < 1% error
        }
      }
    }
  }
}

#ifdef SINGULARITY_USE_EOSPAC

using singularity::EOSPAC;

SCENARIO("SpinerEOS construction from EOSPAC", "[SpinerEOS][Constructor][EOSPAC]") {
  GIVEN("An EOSPAC EOS for aluminum") {
    // Aluminum matid = 3720 (commonly available in EOSPAC)
    constexpr int matid = 3720;
    EOSPAC eospac_eos(matid);

    // Set up grid parameters - re-grid EOSPAC onto custom grid
    SpinerTableGridParams params;
    params.rhoMin = 0.1;   // g/cc
    params.rhoMax = 100.0; // g/cc
    params.TMin = 1e2;     // K
    params.TMax = 1e6;     // K

    // Estimate sie bounds from EOSPAC
    Real sie_at_min =
        eospac_eos.InternalEnergyFromDensityTemperature(params.rhoMin, params.TMin);
    Real sie_at_max =
        eospac_eos.InternalEnergyFromDensityTemperature(params.rhoMax, params.TMax);

    params.sieMin = std::min(sie_at_min, sie_at_max) * 0.9;
    params.sieMax = std::max(sie_at_min, sie_at_max) * 1.1;

    // Coarser grid for faster testing
    params.numRhoPerDecade = 50;
    params.numTPerDecade = 50;
    params.numSiePerDecade = 50;

    // Use piecewise grids (more realistic)
    params.piecewiseRho = true;
    params.piecewiseT = true;
    params.piecewiseSie = true;

    params.matid = 3720;
    params.rhoNormal = 2.7; // Normal density of aluminum

    WHEN("We construct a SpinerEOS from EOSPAC") {
      SpinerEOSDependsRhoSie spiner_eos(eospac_eos, params);

      THEN("The SpinerEOS should interpolate consistently with EOSPAC") {
        // Test at several points
        Real rho = 2.7;  // Normal density
        Real T = 3000.0; // 3000 K

        Real P_eospac = eospac_eos.PressureFromDensityTemperature(rho, T);
        Real P_spiner = spiner_eos.PressureFromDensityTemperature(rho, T);

        Real sie_eospac = eospac_eos.InternalEnergyFromDensityTemperature(rho, T);
        Real sie_spiner = spiner_eos.InternalEnergyFromDensityTemperature(rho, T);

        INFO("EOSPAC P: " << P_eospac << "  Spiner P: " << P_spiner);
        INFO("EOSPAC sie: " << sie_eospac << "  Spiner sie: " << sie_spiner);

        // Tolerance higher than IdealGas due to table interpolation differences
        CHECK(isClose(P_spiner, P_eospac, 0.05)); // 5% tolerance
        CHECK(isClose(sie_spiner, sie_eospac, 0.05));
      }

      AND_THEN("Temperature inversion should work") {
        Real rho = 5.0;
        Real T_orig = 5000.0;

        // Get sie from EOSPAC
        Real sie = eospac_eos.InternalEnergyFromDensityTemperature(rho, T_orig);

        // Invert using SpinerEOS
        Real T_spiner = spiner_eos.TemperatureFromDensityInternalEnergy(rho, sie);
        Real T_eospac = eospac_eos.TemperatureFromDensityInternalEnergy(rho, sie);

        INFO("Original T: " << T_orig);
        INFO("EOSPAC T: " << T_eospac << "  Spiner T: " << T_spiner);

        CHECK(isClose(T_spiner, T_eospac, 0.05));
      }

      AND_THEN("Bulk modulus should be reasonable") {
        Real rho = 2.7;
        Real T = 3000.0;

        Real bmod_eospac = eospac_eos.BulkModulusFromDensityTemperature(rho, T);
        Real bmod_spiner = spiner_eos.BulkModulusFromDensityTemperature(rho, T);

        INFO("EOSPAC bmod: " << bmod_eospac << "  Spiner bmod: " << bmod_spiner);
        CHECK(bmod_spiner > 0);
        // Bulk modulus can differ more due to derivative approximations
        CHECK(isClose(bmod_spiner, bmod_eospac, 0.2)); // 20% tolerance
      }

      AND_THEN("Material properties should be preserved") {
        CHECK(spiner_eos.matid() == params.matid);

        // EOSPAC should provide mean atomic mass and number
        Real abar_eospac = eospac_eos.MeanAtomicMass();
        Real abar_spiner = spiner_eos.MeanAtomicMass();

        Real zbar_eospac = eospac_eos.MeanAtomicNumber();
        Real zbar_spiner = spiner_eos.MeanAtomicNumber();

        INFO("EOSPAC Abar: " << abar_eospac << " Spiner Abar: " << abar_spiner);
        INFO("EOSPAC Zbar: " << zbar_eospac << " Spiner Zbar: " << zbar_spiner);

        CHECK(isClose(abar_spiner, abar_eospac, 1e-6));
        CHECK(isClose(zbar_spiner, zbar_eospac, 1e-6));
      }
    }
  }
}

#endif // SINGULARITY_USE_EOSPAC

SCENARIO("SpinerEOS construction from minimal EOS with root finding",
         "[SpinerEOS][Constructor][MinimalEOS][RootFinding]") {
  GIVEN("A minimal EOS that only provides T(rho,sie) and P(rho,T)") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.4;
    MinimalEOS1 minimal_eos(gm1, Cv);

    // Also create full IdealGas for comparison
    IdealGas ideal_eos(gm1, Cv);

    SpinerTableGridParams params;
    params.rhoMin = 1e-2;
    params.rhoMax = 1e2;
    params.TMin = 1e2;
    params.TMax = 1e4;
    params.sieMin = Cv * params.TMin;
    params.sieMax = Cv * params.TMax;

    params.numRhoPerDecade = 50;
    params.numTPerDecade = 50;
    params.numSiePerDecade = 50;
    params.piecewiseRho = false;
    params.piecewiseT = false;
    params.piecewiseSie = false;

    params.matid = 2001;

    WHEN("We construct a SpinerEOS from the minimal EOS") {
      SpinerEOSDependsRhoSie spiner_eos(minimal_eos, params);

      THEN("Root finding should successfully invert T to get sie") {
        Real rho = 1.0;
        Real T = 1000.0;

        Real P_spiner = spiner_eos.PressureFromDensityTemperature(rho, T);
        Real P_ideal = ideal_eos.PressureFromDensityTemperature(rho, T);

        Real sie_spiner = spiner_eos.InternalEnergyFromDensityTemperature(rho, T);
        Real sie_ideal = ideal_eos.InternalEnergyFromDensityTemperature(rho, T);

        INFO("Root-finding path - Ideal P: " << P_ideal << " Spiner P: " << P_spiner);
        INFO("Root-finding path - Ideal sie: " << sie_ideal
                                               << " Spiner sie: " << sie_spiner);

        CHECK(isClose(P_spiner, P_ideal, 0.02));
        CHECK(isClose(sie_spiner, sie_ideal, 0.02));
      }

      AND_THEN("Finite difference derivatives should give reasonable results") {
        Real rho = 1.0;
        Real sie = Cv * 1000.0;

        Real T_spiner = spiner_eos.TemperatureFromDensityInternalEnergy(rho, sie);
        Real T_ideal = ideal_eos.TemperatureFromDensityInternalEnergy(rho, sie);

        INFO("FD derivatives - Ideal T: " << T_ideal << " Spiner T: " << T_spiner);
        CHECK(isClose(T_spiner, T_ideal, 0.02));
      }

      AND_THEN("Material properties should be set") {
        CHECK(spiner_eos.matid() == params.matid);
        CHECK(spiner_eos.MeanAtomicMass() == ideal_eos.MeanAtomicMass());
        CHECK(spiner_eos.MeanAtomicNumber() == ideal_eos.MeanAtomicNumber());
      }
    }
  }
}

SCENARIO("SpinerEOS construction from minimal EOS with chain rule",
         "[SpinerEOS][Constructor][MinimalEOS][ChainRule]") {
  GIVEN("A minimal EOS without P(rho,sie) method") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.4;
    MinimalEOS2 minimal_eos(gm1, Cv);

    IdealGas ideal_eos(gm1, Cv);

    SpinerTableGridParams params;
    params.rhoMin = 1e-2;
    params.rhoMax = 1e2;
    params.TMin = 1e2;
    params.TMax = 1e4;
    params.sieMin = Cv * params.TMin;
    params.sieMax = Cv * params.TMax;

    params.numRhoPerDecade = 50;
    params.numTPerDecade = 50;
    params.numSiePerDecade = 50;
    params.piecewiseRho = false;
    params.piecewiseT = false;
    params.piecewiseSie = false;

    params.matid = 2002;

    WHEN("We construct a SpinerEOS forcing chain rule fallback") {
      SpinerEOSDependsRhoSie spiner_eos(minimal_eos, params);

      THEN("Chain rule conversions should work for dependsRhoSie derivatives") {
        Real rho = 1.0;
        Real sie = Cv * 1000.0;

        Real P_spiner = spiner_eos.PressureFromDensityInternalEnergy(rho, sie);
        Real P_ideal = ideal_eos.PressureFromDensityInternalEnergy(rho, sie);

        Real T_spiner = spiner_eos.TemperatureFromDensityInternalEnergy(rho, sie);
        Real T_ideal = ideal_eos.TemperatureFromDensityInternalEnergy(rho, sie);

        INFO("Chain rule path - Ideal P: " << P_ideal << " Spiner P: " << P_spiner);
        INFO("Chain rule path - Ideal T: " << T_ideal << " Spiner T: " << T_spiner);

        // Tolerance higher due to chain rule approximations
        CHECK(isClose(P_spiner, P_ideal, 0.05));
        CHECK(isClose(T_spiner, T_ideal, 0.02));
      }

      AND_THEN("Bulk modulus should be computed correctly") {
        Real rho = 1.0;
        Real T = 1000.0;

        Real bmod_spiner = spiner_eos.BulkModulusFromDensityTemperature(rho, T);
        Real bmod_ideal = ideal_eos.BulkModulusFromDensityTemperature(rho, T);

        INFO("Ideal bmod: " << bmod_ideal << " Spiner bmod: " << bmod_spiner);
        CHECK(bmod_spiner > 0);
        // Higher tolerance for bulk modulus with chain rule derivatives
        CHECK(isClose(bmod_spiner, bmod_ideal, 0.15));
      }
    }
  }
}

// ============================================================================
// SpinerEOSDependsRhoT Constructor Tests
// ============================================================================

SCENARIO("SpinerEOSDependsRhoT construction from IdealGas",
         "[SpinerEOS][ConstructorRhoT]") {
  GIVEN("An IdealGas EOS") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.4;
    IdealGas ideal_eos(gm1, Cv);

    WHEN("We construct a SpinerEOSDependsRhoT from it") {
      SpinerTableGridParams params;
      params.rhoMin = 1e-3;
      params.rhoMax = 1e3;
      params.TMin = 1e2;
      params.TMax = 1e5;
      params.numRhoPerDecade = 50; // Coarser for faster test
      params.numTPerDecade = 50;
      params.matid = 1001;

      SpinerEOSDependsRhoT spiner_eos(ideal_eos, params);

      THEN("The table is created successfully") {
        // Check metadata
        REQUIRE(spiner_eos.matid() == 1001);
        CHECK(isClose(spiner_eos.rhoMin(), params.rhoMin, 1e-6));
        CHECK(isClose(spiner_eos.rhoMax(), params.rhoMax, 1e-6));
        CHECK(isClose(spiner_eos.TMin(), params.TMin, 1e-6));
        CHECK(isClose(spiner_eos.TMax(), params.TMax, 1e-6));

        // Check material properties
        CHECK(isClose(spiner_eos.MeanAtomicMass(), ideal_eos.MeanAtomicMass(), 1e-6));
        CHECK(isClose(spiner_eos.MeanAtomicNumber(), ideal_eos.MeanAtomicNumber(), 1e-6));
      }

      AND_THEN("We can interpolate pressure and sie from the table") {
        // Test several points within the table
        const Real rho = 1.0;  // g/cm^3
        const Real T = 1000.0; // K

        Real P_ideal = ideal_eos.PressureFromDensityTemperature(rho, T);
        Real P_spiner = spiner_eos.PressureFromDensityTemperature(rho, T);

        Real sie_ideal = ideal_eos.InternalEnergyFromDensityTemperature(rho, T);
        Real sie_spiner = spiner_eos.InternalEnergyFromDensityTemperature(rho, T);

        // Should match closely (within interpolation error)
        CHECK(isClose(P_spiner, P_ideal, 0.01)); // 1% error
        CHECK(isClose(sie_spiner, sie_ideal, 0.01));
      }

      AND_THEN("Temperature inversion works via root finding") {
        const Real rho = 1.0;
        const Real T_target = 1000.0;
        Real sie = ideal_eos.InternalEnergyFromDensityTemperature(rho, T_target);

        Real T_inverted = spiner_eos.TemperatureFromDensityInternalEnergy(rho, sie);

        // Should recover the temperature
        CHECK(isClose(T_inverted, T_target, 0.01));
      }

      AND_THEN("Bulk modulus is computed correctly") {
        const Real rho = 1.0;
        const Real T = 1000.0;

        Real bmod_ideal = ideal_eos.BulkModulusFromDensityTemperature(rho, T);
        Real bmod_spiner = spiner_eos.BulkModulusFromDensityTemperature(rho, T);

        // Bulk modulus computation involves derivatives, so allow more error
        CHECK(isClose(bmod_spiner, bmod_ideal, 0.05)); // 5% error
      }
    }
  }
}

SCENARIO("SpinerEOSDependsRhoT with piecewise grids", "[SpinerEOS][ConstructorRhoT]") {
  GIVEN("An IdealGas EOS and piecewise grid parameters") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.4;
    IdealGas ideal_eos(gm1, Cv);

    SpinerTableGridParams params;
    params.rhoMin = 1e-3;
    params.rhoMax = 1e3;
    params.TMin = 1e2;
    params.TMax = 1e5;
    params.matid = 2002;

    // Enable piecewise grids
    params.piecewiseRho = true;
    params.piecewiseT = true;
    params.numRhoPerDecade = 30;
    params.numTPerDecade = 30;

    WHEN("We construct with piecewise grids") {
      SpinerEOSDependsRhoT spiner_eos(ideal_eos, params);

      THEN("The table works correctly") {
        const Real rho = 1.0;
        const Real T = 1000.0;

        Real P_ideal = ideal_eos.PressureFromDensityTemperature(rho, T);
        Real P_spiner = spiner_eos.PressureFromDensityTemperature(rho, T);

        CHECK(isClose(P_spiner, P_ideal, 0.02));
      }
    }
  }
}

SCENARIO("SpinerEOSDependsRhoT accuracy test", "[SpinerEOS][ConstructorRhoT]") {
  GIVEN("An IdealGas EOS with fine grid") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.4;
    IdealGas ideal_eos(gm1, Cv);

    SpinerTableGridParams params;
    params.rhoMin = 0.1;
    params.rhoMax = 10.0;
    params.TMin = 100.0;
    params.TMax = 10000.0;
    params.numRhoPerDecade = 100; // Fine grid for accuracy
    params.numTPerDecade = 100;
    params.matid = 3003;

    SpinerEOSDependsRhoT spiner_eos(ideal_eos, params);

    WHEN("We test multiple points across the table") {
      int num_tests = 10;
      Real max_P_error = 0.0;
      Real max_sie_error = 0.0;
      Real max_bmod_error = 0.0;

      for (int i = 0; i < num_tests; ++i) {
        Real rho = params.rhoMin *
                   std::pow(params.rhoMax / params.rhoMin, Real(i) / (num_tests - 1));
        Real T =
            params.TMin * std::pow(params.TMax / params.TMin, Real(i) / (num_tests - 1));

        Real P_ideal = ideal_eos.PressureFromDensityTemperature(rho, T);
        Real P_spiner = spiner_eos.PressureFromDensityTemperature(rho, T);
        Real P_error = std::abs(P_spiner - P_ideal) / std::abs(P_ideal);
        max_P_error = std::max(max_P_error, P_error);

        Real sie_ideal = ideal_eos.InternalEnergyFromDensityTemperature(rho, T);
        Real sie_spiner = spiner_eos.InternalEnergyFromDensityTemperature(rho, T);
        Real sie_error = std::abs(sie_spiner - sie_ideal) / std::abs(sie_ideal);
        max_sie_error = std::max(max_sie_error, sie_error);

        Real bmod_ideal = ideal_eos.BulkModulusFromDensityTemperature(rho, T);
        Real bmod_spiner = spiner_eos.BulkModulusFromDensityTemperature(rho, T);
        Real bmod_error = std::abs(bmod_spiner - bmod_ideal) / std::abs(bmod_ideal);
        max_bmod_error = std::max(max_bmod_error, bmod_error);
      }

      THEN("Interpolation errors are small") {
        // With fine grid, errors should be very small
        REQUIRE(max_P_error < 0.005);   // < 0.5%
        REQUIRE(max_sie_error < 0.005); // < 0.5%
        REQUIRE(max_bmod_error < 0.05); // < 5% (uses source EOS method when available)
      }
    }
  }
}

#ifdef SINGULARITY_USE_EOSPAC
SCENARIO("SpinerEOSDependsRhoT construction from EOSPAC",
         "[SpinerEOS][ConstructorRhoT][EOSPAC]") {
  GIVEN("An EOSPAC EOS for aluminum") {
    constexpr int matid = 3720; // Aluminum
    singularity::EOSPAC eospac_eos(matid);

    WHEN("We create a SpinerEOSDependsRhoT from it") {
      SpinerTableGridParams params;
      params.rhoMin = 0.5; // g/cm^3
      params.rhoMax = 20.0;
      params.TMin = 300.0; // K
      params.TMax = 50000.0;
      params.numRhoPerDecade = 40;
      params.numTPerDecade = 40;
      params.matid = matid;

      SpinerEOSDependsRhoT spiner_eos(eospac_eos, params);

      THEN("The table is created and interpolates correctly") {
        const Real rho = 2.7; // Near aluminum normal density
        const Real T = 1000.0;

        Real P_eospac = eospac_eos.PressureFromDensityTemperature(rho, T);
        Real P_spiner = spiner_eos.PressureFromDensityTemperature(rho, T);

        // EOSPAC tables can be complex, so allow reasonable error
        CHECK(isClose(P_spiner, P_eospac, 0.05));
      }

      AND_THEN("Material properties are preserved") {
        CHECK(isClose(spiner_eos.MeanAtomicMass(), eospac_eos.MeanAtomicMass(), 1e-6));
        CHECK(
            isClose(spiner_eos.MeanAtomicNumber(), eospac_eos.MeanAtomicNumber(), 1e-6));
      }
    }
  }
}
#endif // SINGULARITY_USE_EOSPAC

SCENARIO("SpinerEOSDependsRhoT vs SpinerEOSDependsRhoSie comparison",
         "[SpinerEOS][ConstructorRhoT]") {
  GIVEN("An IdealGas EOS") {
    constexpr Real Cv = 2.0;
    constexpr Real gm1 = 0.4;
    IdealGas ideal_eos(gm1, Cv);

    SpinerTableGridParams params;
    params.rhoMin = 1e-2;
    params.rhoMax = 1e2;
    params.TMin = 1e2;
    params.TMax = 1e4;
    params.sieMin = Cv * params.TMin;
    params.sieMax = Cv * params.TMax;
    params.numRhoPerDecade = 50;
    params.numTPerDecade = 50;
    params.numSiePerDecade = 50;
    params.matid = 4004;

    WHEN("We construct both RhoT and RhoSie variants") {
      SpinerEOSDependsRhoT rhot_eos(ideal_eos, params);
      SpinerEOSDependsRhoSie rhosie_eos(ideal_eos, params);

      THEN("Their (rho,T) lookups should match") {
        const Real rho = 1.0;
        const Real T = 1000.0;

        Real P_rhot = rhot_eos.PressureFromDensityTemperature(rho, T);
        Real P_rhosie = rhosie_eos.PressureFromDensityTemperature(rho, T);

        Real sie_rhot = rhot_eos.InternalEnergyFromDensityTemperature(rho, T);
        Real sie_rhosie = rhosie_eos.InternalEnergyFromDensityTemperature(rho, T);

        Real bmod_rhot = rhot_eos.BulkModulusFromDensityTemperature(rho, T);
        Real bmod_rhosie = rhosie_eos.BulkModulusFromDensityTemperature(rho, T);

        // Both use the same (rho,T) table construction, should be very close
        CHECK(isClose(P_rhot, P_rhosie, 1e-6));
        CHECK(isClose(sie_rhot, sie_rhosie, 1e-6));
        CHECK(isClose(bmod_rhot, bmod_rhosie, 1e-6));
      }

      AND_THEN("RhoSie has additional (rho,sie) capabilities") {
        // RhoSie can do this efficiently with table lookup
        const Real rho = 1.0;
        const Real sie = Cv * 1000.0;

        Real P_rhosie = rhosie_eos.PressureFromDensityInternalEnergy(rho, sie);

        // RhoT can also do this but via root finding
        // (TemperatureFromDensityInternalEnergy)
        Real T_rhot = rhot_eos.TemperatureFromDensityInternalEnergy(rho, sie);
        Real P_rhot = rhot_eos.PressureFromDensityTemperature(rho, T_rhot);

        // Should give same answer, though RhoSie is more direct
        CHECK(isClose(P_rhot, P_rhosie, 0.01));
      }
    }
  }
}

SCENARIO("SpinerEOSDependsRhoT from minimal EOS", "[SpinerEOS][ConstructorRhoT]") {
  GIVEN("A minimal EOS that only provides basic methods") {
    // This class intentionally lacks gamma, cv, and other optimizations
    // to test that finite differences work correctly
    class MinimalEOS {
     private:
      IdealGas ideal_;

     public:
      MinimalEOS(Real gm1, Real Cv) : ideal_(gm1, Cv) {}

      PORTABLE_INLINE_FUNCTION
      Real InternalEnergyFromDensityTemperature(Real rho, Real T) const {
        return ideal_.InternalEnergyFromDensityTemperature(rho, T);
      }

      PORTABLE_INLINE_FUNCTION
      Real TemperatureFromDensityInternalEnergy(Real rho, Real sie) const {
        return ideal_.TemperatureFromDensityInternalEnergy(rho, sie);
      }

      PORTABLE_INLINE_FUNCTION
      Real PressureFromDensityTemperature(Real rho, Real T) const {
        return ideal_.PressureFromDensityTemperature(rho, T);
      }

      // No gamma, cv, or bulk modulus methods - forces finite differences
    };

    MinimalEOS minimal_eos(0.4, 2.0);

    WHEN("We construct a SpinerEOSDependsRhoT from it") {
      SpinerTableGridParams params;
      params.rhoMin = 0.1;
      params.rhoMax = 10.0;
      params.TMin = 100.0;
      params.TMax = 5000.0;
      params.numRhoPerDecade = 40;
      params.numTPerDecade = 40;
      params.matid = 5005;

      SpinerEOSDependsRhoT spiner_eos(minimal_eos, params);

      THEN("The table is created successfully using finite differences") {
        const Real rho = 1.0;
        const Real T = 1000.0;

        // Create reference IdealGas for comparison
        IdealGas ideal_ref(0.4, 2.0);
        Real P_ref = ideal_ref.PressureFromDensityTemperature(rho, T);
        Real P_spiner = spiner_eos.PressureFromDensityTemperature(rho, T);

        // Should still be accurate despite using finite differences
        CHECK(isClose(P_spiner, P_ref, 0.02));
      }
    }
  }
}

#endif // SINGULARITY_USE_SPINER_WITH_HDF5
