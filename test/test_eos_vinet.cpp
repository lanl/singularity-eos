//------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#ifndef CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#endif

#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <test/eos_unit_test_helpers.hpp>

using singularity::EOS;
using singularity::Vinet;

PORTABLE_INLINE_FUNCTION Real QuadFormulaMinus(Real a, Real b, Real c) {
  return (-b - std::sqrt(b * b - 4 * a * c)) / (2 * a);
}

SCENARIO("Vinet EOS", "[VectorEOS][VinetEOS]") {
  GIVEN("Parameters for a Vinet EOS") {
    // Unit conversions
    constexpr Real cm = 1.;
    constexpr Real us = 1e-06;
    constexpr Real Mbcc_per_g = 1e12;
    // Vinet parameters for copper
    constexpr Real rho0 = 8.93;
    constexpr Real T0 = 298.0;
    constexpr Real B0 = 1.3448466 * Mbcc_per_g;
    constexpr Real BP0 = 4.956;
    constexpr Real A0 = 5.19245e-05;
    constexpr Real Cv0 = 0.383e-05 * Mbcc_per_g;
    constexpr Real E0 = 0.0;
    constexpr Real S0 = 5.05e-04 * Mbcc_per_g;
    constexpr int dnum = 38;
    constexpr Real d2to40[38]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	    0.,0.,0.,0.,0.,0.,0.,0.};
    // Create the EOS
    EOS host_eos = Vinet(rho0, T0, B0, BP0, A0, Cv0, E0, S0, d2to40);
    EOS eos = host_eos.GetOnDevice();

    eos.PrintParams();
    
    GIVEN("Densities and energies") {
      constexpr int num = 4;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create Kokkos views on device for the input arrays
      Kokkos::View<Real[num]> v_density("density");
      Kokkos::View<Real[num]> v_energy("density");
#else
      // Otherwise just create arrays to contain values and create pointers to
      // be passed to the functions in place of the Kokkos views
      std::array<Real, num> density;
      std::array<Real, num> energy;
      auto v_density = density.data();
      auto v_energy = energy.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

      // Populate the input arrays
      portableFor(
          "Initialize density and energy", 0, 1, PORTABLE_LAMBDA(int i) {
            v_density[0] = 8.0;
            v_density[1] = 8.5;
            v_density[2] = 9.0;
            v_density[3] = 8.93;
            v_energy[0] = 2.e8;
            v_energy[1] = 4.6e9;
            v_energy[2] = 4.7e9;
            v_energy[3] = 0.0;
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create host-side mirrors of the inputs and copy the inputs. These are
      // just used for the comparisons
      auto density = Kokkos::create_mirror_view(v_density);
      auto energy = Kokkos::create_mirror_view(v_energy);
      Kokkos::deep_copy(density, v_density);
      Kokkos::deep_copy(energy, v_energy);
#endif // PORTABILITY_STRATEGY_KOKKOS

      // Gold standard values for a subset of lookups
      constexpr std::array<Real, num> pressure_true{
          -1.282094800000000e+11, 1.998504088912181e+10, 9.595823319513451e+10,
           B0 - B0*B0*A0*A0/Cv0/rho0*T0};
      constexpr std::array<Real, num> bulkmodulus_true{
          9.990648504000005e+11, 1.460692677162573e+12, 1.851227213843747e+12,
          (B0 - B0*B0*A0*A0/Cv0/rho0*T0)};
      constexpr std::array<Real, num> temperature_true{
          5.590966057441253e+02, 4.285483028720627e+02, 3.241096605744125e+02, T0};
      constexpr std::array<Real, num> gruneisen_true{B0*A0/Cv0/rho0, 2.007944444444444e+00,
                                                     1.927000000000000e+00, B0*A0/Cv0/rho0};

#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create device views for outputs and mirror those views on the host
      Kokkos::View<Real[num]> v_temperature("Temperature");
      Kokkos::View<Real[num]> v_pressure("Pressure");
      Kokkos::View<Real[num]> v_bulkmodulus("bmod");
      Kokkos::View<Real[num]> v_gruneisen("Gamma");
      auto h_temperature = Kokkos::create_mirror_view(v_temperature);
      auto h_pressure = Kokkos::create_mirror_view(v_pressure);
      auto h_bulkmodulus = Kokkos::create_mirror_view(v_bulkmodulus);
      auto h_gruneisen = Kokkos::create_mirror_view(v_gruneisen);
#else
      // Create arrays for the outputs and then pointers to those arrays that
      // will be passed to the functions in place of the Kokkos views
      std::array<Real, num> h_temperature;
      std::array<Real, num> h_pressure;
      std::array<Real, num> h_bulkmodulus;
      std::array<Real, num> h_gruneisen;
      // Just alias the existing pointers
      auto v_temperature = h_temperature.data();
      auto v_pressure = h_pressure.data();
      auto v_bulkmodulus = h_bulkmodulus.data();
      auto v_gruneisen = h_gruneisen.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

      WHEN("A P(rho, e) lookup is performed") {
        eos.PressureFromDensityInternalEnergy(v_density, v_energy, v_pressure, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_pressure, v_pressure);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned P(rho, e) should be equal to the true value") {
          array_compare(num, density, energy, h_pressure, pressure_true, "Density",
                        "Energy");
        }
      }

      WHEN("A B_S(rho, e) lookup is performed") {
        eos.BulkModulusFromDensityInternalEnergy(v_density, v_energy, v_bulkmodulus, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_bulkmodulus, v_bulkmodulus);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned B_S(rho, e) should be equal to the true value") {
          array_compare(num, density, energy, h_bulkmodulus, bulkmodulus_true, "Density",
                        "Energy");
        }
      }

      WHEN("A T(rho, e) lookup is performed") {
        eos.TemperatureFromDensityInternalEnergy(v_density, v_energy, v_temperature, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_temperature, v_temperature);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned T(rho, e) should be equal to the true value") {
          array_compare(num, density, energy, h_temperature, temperature_true, "Density",
                        "Energy");
        }
      }

      WHEN("A Gamma(rho, e) lookup is performed") {
        eos.GruneisenParamFromDensityInternalEnergy(v_density, v_energy, v_gruneisen,
                                                    num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_gruneisen, v_gruneisen);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned Gamma(rho, e) should be equal to the true value") {
          array_compare(num, density, energy, h_gruneisen, gruneisen_true, "Density",
                        "Energy");
        }
      }
    }
  }
}
