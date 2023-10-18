//------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------

#ifndef _SINGULARITY_EOS_TEST_TEST_HELPERS_
#define _SINGULARITY_EOS_TEST_TEST_HELPERS_

#ifndef CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#endif
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <limits>
#include <sstream>
#include <iomanip>

// typename demangler
#ifdef __GNUG__
#include <cstdlib>
#include <cxxabi.h>
#include <memory>

#include <ports-of-call/portability.hpp>

inline std::string demangle(const char *name) {

  int status = -4; // some arbitrary value to eliminate the compiler warning

  // enable c++11 by passing the flag -std=c++11 to g++
  std::unique_ptr<char, void (*)(void *)> res{
      abi::__cxa_demangle(name, NULL, NULL, &status), std::free};

  return (status == 0) ? res.get() : name;
}

#else

// does nothing if not g++
inline std::string demangle(const char *name) { return name; }

#endif

#include <ports-of-call/portability.hpp>

using singularity::EOS;

PORTABLE_INLINE_FUNCTION bool isClose(Real a, Real b, Real eps = 5e-2) {
  return fabs(b - a) / (fabs(a + b) + 1e-20) <= eps;
}

PORTABLE_INLINE_FUNCTION Real myAtan(Real x, Real shift, Real scale, Real offset) {
  return scale * atan(x - shift) + offset;
}

// Function for comparing arrays and outputting the important information
template <typename X, typename Y, typename Z, typename ZT, typename XN, typename YN>
inline void array_compare(int num, X &&x, Y &&y, Z &&z, ZT &&ztrue, XN xname, YN yname,
                          Real tol = 1e-12) {
  assert(num > 0);
  using underlying_t = typename std::remove_cv<typename std::remove_reference<decltype(x[0])>::type>::type;
  for (int i = 0; i < num; i++) {
    auto s = std::ostringstream{};
    s << std::setprecision(std::numeric_limits<underlying_t>::max_digits10) 
      << std::scientific
      << "i: " << i << ", " << xname << ": " << x[i] << ", " 
      << yname << ": " << y[i] << ", Value: " << z[i] 
      << ", True Value: " << ztrue[i];
    INFO(s.str());
    CHECK(isClose(z[i], ztrue[i], 1e-12));
  }
}

inline void compare_two_eoss(const EOS &test_e, const EOS &ref_e) {
  // compare all individual member functions with 1 as inputs,
  // this function is meant to catch mis-implementations of
  // modifiers that can be initialized in such a way as to
  // be equivalent of an unmodified eos. Best used with analytic
  // eoss.
  REQUIRE(isClose(test_e.TemperatureFromDensityInternalEnergy(1, 1),
                  ref_e.TemperatureFromDensityInternalEnergy(1, 1), 1.e-15));
  REQUIRE(isClose(test_e.InternalEnergyFromDensityTemperature(1, 1),
                  ref_e.InternalEnergyFromDensityTemperature(1, 1), 1.e-15));
  REQUIRE(isClose(test_e.PressureFromDensityInternalEnergy(1, 1),
                  ref_e.PressureFromDensityInternalEnergy(1, 1), 1.e-15));
  REQUIRE(isClose(test_e.SpecificHeatFromDensityInternalEnergy(1, 1),
                  ref_e.SpecificHeatFromDensityInternalEnergy(1, 1), 1.e-15));
  REQUIRE(isClose(test_e.BulkModulusFromDensityInternalEnergy(1, 1),
                  ref_e.BulkModulusFromDensityInternalEnergy(1, 1), 1.e-15));
  REQUIRE(isClose(test_e.GruneisenParamFromDensityInternalEnergy(1, 1),
                  ref_e.GruneisenParamFromDensityInternalEnergy(1, 1), 1.e-15));
  REQUIRE(isClose(test_e.PressureFromDensityTemperature(1, 1),
                  ref_e.PressureFromDensityTemperature(1, 1), 1.e-15));
  REQUIRE(isClose(test_e.SpecificHeatFromDensityTemperature(1, 1),
                  ref_e.SpecificHeatFromDensityTemperature(1, 1), 1.e-15));
  REQUIRE(isClose(test_e.BulkModulusFromDensityTemperature(1, 1),
                  ref_e.BulkModulusFromDensityTemperature(1, 1), 1.e-15));
  REQUIRE(isClose(test_e.GruneisenParamFromDensityTemperature(1, 1),
                  ref_e.GruneisenParamFromDensityTemperature(1, 1), 1.e-15));
  REQUIRE(isClose(test_e.MinimumDensity(), ref_e.MinimumDensity(), 1.e-15));
  REQUIRE(isClose(test_e.MinimumTemperature(), ref_e.MinimumTemperature(), 1.e-15));
  return;
}

// Macro that checks for an exception or is a no-op depending on
// whether or not a non-serial backend is supplied
#ifdef PORTABILITY_STRATEGY_NONE
#define REQUIRE_MAYBE_THROWS(...) REQUIRE_THROWS(__VA_ARGS__)
#else
#define REQUIRE_MAYBE_THROWS(...) ((void)0)
#endif // PORTABILITY_STRATEGY_NONE

#endif // _SINGULARITY_EOS_TEST_TEST_HELPERS_
