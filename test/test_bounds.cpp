//------------------------------------------------------------------------------
// © 2021-2024. Triad National Security, LLC. All rights reserved.  This
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

#ifdef SINGULARITY_USE_SPINER

#include <cstdio>
#include <limits>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/fast-math/logs.hpp>
#include <singularity-eos/base/spiner_table_utils.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

constexpr int NGRIDS = 3;
using Bounds = singularity::table_utils::Bounds<NGRIDS>;
using RegularGrid1D = Spiner::RegularGrid1D<Real>;
using Grid_t = Bounds::Grid_t;

constexpr Real rho_normal = 3; // g/cm^3
constexpr Real T_normal = 300; // Kelvin
constexpr Real rho_min = 1e-6;
constexpr Real rho_max = 1e3;
constexpr Real T_min = 1;
constexpr Real T_max = 1e9;
constexpr Real T_split = 1e4;
constexpr int N_per_decade_fine = 200;
constexpr Real N_factor = 5;

constexpr Real REAL_TOL = std::numeric_limits<Real>::epsilon() * 1e3;

SCENARIO("Bounds can compute number of points from points per decade", "[Bounds]") {
  WHEN("We compute the number of points from points per decade") {
    int np = Bounds::getNumPointsFromPPD(T_min, T_max, N_per_decade_fine);
    constexpr int NDECADES = 9;
    THEN("We get the right number") { REQUIRE(np == NDECADES * N_per_decade_fine); }
  }
}

SCENARIO("Linear bounds in the bounds object", "[Bounds]") {
  WHEN("We compute linear bounds ") {
    constexpr Real min = -2;
    constexpr Real max = 5;
    constexpr int N = 21;
    Bounds lin_bounds(min, max, N);
    THEN("The min and max and point number correct") {
      REQUIRE(lin_bounds.grid.min() == min);
      REQUIRE(lin_bounds.grid.max() == max);
      REQUIRE(lin_bounds.grid.nPoints() == N);
    }
  }
}

SCENARIO("Logarithmic, single-grid bounds in the bounds object", "[Bounds]") {
  WHEN("We compute logarithmic single-grid bounds") {
    int np = Bounds::getNumPointsFromPPD(rho_min, rho_max, N_per_decade_fine);
    Bounds lRhoBounds(rho_min, rho_max, np, true, 0.0, rho_normal);
    THEN("The lower and upper bounds are right") {
      REQUIRE(std::abs(lRhoBounds.grid.min() - singularity::FastMath::log10(rho_min)) <=
              REAL_TOL);
      REQUIRE(lRhoBounds.grid.max() <=
              singularity::FastMath::log10(rho_max)); // shifted due to anchor
      AND_THEN("The anchor is on the mesh") {
        Real lanchor = singularity::FastMath::log10(rho_normal);
        int ianchor;
        Spiner::weights_t<Real> w;
        lRhoBounds.grid.weights(lanchor, ianchor, w);
        printf("%.14e, %.14e, %.14e\n",
               std::abs(w[0] - 1), std::abs(w[1]), REAL_TOL);
        REQUIRE(std::abs(w[0] - 1) <= REAL_TOL);
        REQUIRE(std::abs(w[1]) <= REAL_TOL);
      }
    }
  }
}

SCENARIO("Logarithmic, piecewise bounds in bounds object", "[Bounds]") {
  WHEN("We compute a piecewise bounds object with three grids") {
    Bounds bnds(Bounds::ThreeGrids(), rho_min, rho_max, rho_normal, 0.5,
                N_per_decade_fine, N_factor, N_factor, true);
    THEN("The bounds are right") {
      Real lrmin = singularity::FastMath::log10(rho_min);
      Real lrmax = singularity::FastMath::log10(rho_max);
      REQUIRE(std::abs(bnds.grid.min() - lrmin) <= REAL_TOL);
      REQUIRE(std::abs(bnds.grid.max() - lrmax) <= REAL_TOL);
      REQUIRE(bnds.grid.nGrids() == 3);
      AND_THEN(
          "The total number of points is less than a uniform fine spacing would imply") {
        REQUIRE(bnds.grid.nPoints() < N_per_decade_fine * (lrmax - lrmin));
        AND_THEN("The anchor is on the mesh") {
          Real lanchor = singularity::FastMath::log10(rho_normal);
          int ianchor;
          Spiner::weights_t<Real> w;
          bnds.grid.weights(lanchor, ianchor, w);
          REQUIRE(std::abs(w[0] - 1) <= REAL_TOL);
          REQUIRE(std::abs(w[1]) <= REAL_TOL);
        }
      }
    }
  }
  WHEN("We compute a piecewise bounds object with two grids") {
    Bounds bnds(Bounds::TwoGrids(), T_min, T_max, T_normal, T_split, N_per_decade_fine,
                N_factor, true);
    THEN("The bounds are right") {
      Real ltmin = singularity::FastMath::log10(T_min);
      Real ltmax = singularity::FastMath::log10(T_max);
      REQUIRE(std::abs(bnds.grid.min() - ltmin) <= REAL_TOL);
      REQUIRE(std::abs(bnds.grid.max() - ltmax) <= REAL_TOL);
      REQUIRE(bnds.grid.nGrids() == 2);
      AND_THEN(
          "The total number of points is less than a uniform fine spacing would imply") {
        REQUIRE(bnds.grid.nPoints() < N_per_decade_fine * (ltmax - ltmin));
        AND_THEN("The anchor is on the mesh") {
          Real lanchor = singularity::FastMath::log10(T_normal);
          int ianchor;
          Spiner::weights_t<Real> w;
          bnds.grid.weights(lanchor, ianchor, w);
          REQUIRE(std::abs(w[0] - 1) <= REAL_TOL);
          REQUIRE(std::abs(w[1]) <= REAL_TOL);
        }
      }
    }
  }
}
#endif // SINGULARITY_USE_SPINER
