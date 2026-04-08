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
using singularity::SpinerEOSDependsRhoSie;
using singularity::SpinerTableGridParams;

SCENARIO("SpinerEOS construction from IdealGas",
         "[SpinerEOS][Constructor][IdealGas]") {
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
        CHECK(bmod_spiner > 0); // Should be positive
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

#endif // SINGULARITY_USE_SPINER_WITH_HDF5
