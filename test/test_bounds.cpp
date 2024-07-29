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

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/table_bounds.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

constexpr int NGRIDS = 3;
using Bounds = singularity::Bounds<NGRIDS>;
using RegularGrid1D = Spiner::RegularGrid1D<Real>;
using Grid_t = Bounds::Grid_t;

constexpr Real rho_normal = 3; // g/cm^3
constexpr Real T_normal = 300; // Kelvin
constexpr Real rho_min = 1e-6;
constexpr Real rho_max = 1e3;
constexpr Real T_min = 1;
constexpr Real T_max = 1e9;
constexpr int N_per_decade_fine = 10;
constexpr Real N_factor = 5;

SCENARIO("linear bounds") {
  WHEN("We compute linear bounds ") {
    constexpr Real min = -2;
    constexpr Real max = 5;
    constexpr int N = 21;
    Bounds lin_bounds(min, max, N);
    THEN("The min and max and point number correct") {
      REQUIRE( lin_bounds.grid.min() == min );
      REQUIRE( lin_bounds.grid.max() == max );
      REQUIRE( lin_bounds.grid.nPoints() == N );
    }
  }
}
