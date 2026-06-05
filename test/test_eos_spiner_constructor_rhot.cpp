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

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <test/eos_unit_test_helpers.hpp>

#ifdef SINGULARITY_USE_SPINER_WITH_HDF5

using singularity::IdealGas;
using singularity::SpinerEOSDependsRhoT;
using singularity::SpinerEOSDependsRhoSie;
using singularity::SpinerTableGridParams;

SCENARIO("SpinerEOSDependsRhoT construction from IdealGas", "[SpinerEOS][ConstructorRhoT]") {
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
        CHECK(isClose(P_spiner, P_ideal, 0.01));   // 1% error
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
        Real rho = params.rhoMin * std::pow(params.rhoMax / params.rhoMin,
                                             Real(i) / (num_tests - 1));
        Real T = params.TMin * std::pow(params.TMax / params.TMin,
                                        Real(i) / (num_tests - 1));

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
        REQUIRE(max_P_error < 0.005);    // < 0.5%
        REQUIRE(max_sie_error < 0.005);  // < 0.5%
        REQUIRE(max_bmod_error < 0.05);  // < 5% (uses source EOS method when available)
      }
    }
  }
}

#ifdef SINGULARITY_USE_EOSPAC
SCENARIO("SpinerEOSDependsRhoT construction from EOSPAC", "[SpinerEOS][ConstructorRhoT][EOSPAC]") {
  GIVEN("An EOSPAC EOS for aluminum") {
    constexpr int matid = 3720; // Aluminum
    singularity::EOSPAC eospac_eos(matid);

    WHEN("We create a SpinerEOSDependsRhoT from it") {
      SpinerTableGridParams params;
      params.rhoMin = 0.5;   // g/cm^3
      params.rhoMax = 20.0;
      params.TMin = 300.0;   // K
      params.TMax = 50000.0;
      params.numRhoPerDecade = 40;
      params.numTPerDecade = 40;
      params.matid = matid;

      SpinerEOSDependsRhoT spiner_eos(eospac_eos, params);

      THEN("The table is created and interpolates correctly") {
        const Real rho = 2.7;     // Near aluminum normal density
        const Real T = 1000.0;

        Real P_eospac = eospac_eos.PressureFromDensityTemperature(rho, T);
        Real P_spiner = spiner_eos.PressureFromDensityTemperature(rho, T);

        // EOSPAC tables can be complex, so allow reasonable error
        CHECK(isClose(P_spiner, P_eospac, 0.05));
      }

      AND_THEN("Material properties are preserved") {
        CHECK(isClose(spiner_eos.MeanAtomicMass(), eospac_eos.MeanAtomicMass(), 1e-6));
        CHECK(isClose(spiner_eos.MeanAtomicNumber(), eospac_eos.MeanAtomicNumber(), 1e-6));
      }
    }
  }
}
#endif // SINGULARITY_USE_EOSPAC

SCENARIO("SpinerEOSDependsRhoT vs SpinerEOSDependsRhoSie comparison", "[SpinerEOS][ConstructorRhoT]") {
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

        // RhoT can also do this but via root finding (TemperatureFromDensityInternalEnergy)
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
