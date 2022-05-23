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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#ifndef CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#endif

#include <singularity-eos/eos/eos.hpp>
#include <test/eos_unit_test_helpers.hpp>

using singularity::EOS;
using singularity::Gruneisen;

PORTABLE_INLINE_FUNCTION Real QuadFormulaMinus(Real a, Real b, Real c) {
  return (-b - std::sqrt(b * b - 4 * a * c)) / (2 * a);
}

SCENARIO("Gruneisen EOS", "[VectorEOS][GruneisenEOS]") {
  GIVEN("Parameters for a Gruneisen EOS") {
    // Unit conversions
    constexpr Real cm = 1.;
    constexpr Real us = 1e-06;
    constexpr Real Mbcc_per_g = 1e12;
    // Gruneisen parameters for copper
    constexpr Real C0 = 0.394 * cm / us;
    constexpr Real S1 = 1.489;
    constexpr Real S2 = 0.;
    constexpr Real S3 = 0.;
    constexpr Real Gamma0 = 2.02;
    constexpr Real b = 0.47;
    constexpr Real rho0 = 8.93;
    constexpr Real T0 = 298.;
    constexpr Real P0 = 0.;
    constexpr Real Cv = 0.383e-05 * Mbcc_per_g;
    // Create the EOS
    EOS host_eos = Gruneisen(C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv);
    EOS eos = host_eos.GetOnDevice();
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
            v_density[1] = 9.0;
            v_density[2] = 9.5;
            v_density[3] = 0.;
            v_energy[0] = 1.e9;
            v_energy[1] = 5.e8;
            v_energy[2] = 1.e8;
            v_energy[3] = 0.;
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
          P0 - C0 * C0 * rho0};
      constexpr std::array<Real, num> bulkmodulus_true{
          9.990648504000005e+11, 1.460692677162573e+12, 1.851227213843747e+12,
          Gamma0 * (P0 - C0 * C0 * rho0)};
      constexpr std::array<Real, num> temperature_true{
          5.590966057441253e+02, 4.285483028720627e+02, 3.241096605744125e+02, T0};
      constexpr std::array<Real, num> gruneisen_true{Gamma0, 2.007944444444444e+00,
                                                     1.927000000000000e+00, Gamma0};

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
        THEN("The returned B_S(rho, e) should be equal to the true value") {
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

SCENARIO("Aluminum Gruneisen EOS Sound Speed and Pressure Comparison", "[GruneisenEOS]") {
  GIVEN("Parameters for a Gruneisen EOS") {
    // Unit conversions
    constexpr Real mm = 10.;
    constexpr Real cm = 1.;
    constexpr Real us = 1.e-06;
    constexpr Real Mbar = 1.e12;
    constexpr Real Mbcc_per_g = 1e12;
    // Gruneisen parameters for copper
    constexpr Real C0 = 0.535 * cm / us;
    constexpr Real S1 = 1.34;
    constexpr Real S2 = 0.;
    constexpr Real S3 = 0.;
    constexpr Real Gamma0 = 1.97;
    constexpr Real b = 0.;
    constexpr Real rho0 = 2.714000;
    constexpr Real T0 = 298.;
    constexpr Real P0 = 1e-06 * Mbar;
    constexpr Real Cv = 0.383e-05 * Mbcc_per_g;
    // Create the EOS
    EOS host_eos = Gruneisen(C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv);
    EOS eos = host_eos.GetOnDevice();
    GIVEN("Density and energy") {
      constexpr Real density = 5.92418956756592;            // g/cm^3
      constexpr Real energy = 792486007.804619;             // erg/g
      constexpr Real true_pres = 2.620656373250729;         // Mbar
      constexpr Real true_sound_speed = 1.5247992468363685; // cm/us
      WHEN("A P(rho, e) lookup is performed") {
        Real pres = eos.PressureFromDensityInternalEnergy(density, energy);
        THEN("The correct pressure should be returned") {
          pres = pres / Mbar;
          INFO("Density: " << density << "  Energy: " << energy << "  Pressure: " << pres
                           << " Mbar  True pressure: " << true_pres << " Mbar");
          REQUIRE(isClose(pres, true_pres, 1e-12));
        }
      }
      WHEN("A B_S(rho, e) lookup is performed") {
        const Real bulk_modulus =
            eos.BulkModulusFromDensityInternalEnergy(density, energy);
        THEN("The correct sound speed should be computed") {
          const Real sound_speed = std::sqrt(bulk_modulus / density) / (cm / us);
          INFO("Density: " << density << "  Energy: " << energy
                           << "  Sound speed: " << sound_speed
                           << " cm/us  True sound speed: " << true_sound_speed
                           << " cm/us");
          REQUIRE(isClose(sound_speed, true_sound_speed, 1e-12));
        }
      }
      WHEN("A finite difference approximation is used for the bulk modulus") {
        // Bulk modulus approximation:
        //  B_S = rho * dPdr_e + P / rho * dPde_r
        constexpr Real drho = 1e-06 * density;
        constexpr Real de = 1e-06 * energy;
        const Real P1 = eos.PressureFromDensityInternalEnergy(density, energy);
        Real P2 = eos.PressureFromDensityInternalEnergy(density + drho, energy);
        const Real dPdr_e = (P2 - P1) / drho;
        P2 = eos.PressureFromDensityInternalEnergy(density, energy + de);
        const Real dPde_r = (P2 - P1) / de;
        const Real bmod_approx = density * dPdr_e + P1 / density * dPde_r;
        THEN("The finite difference solution should approximate the exact solution") {
          const Real bulk_modulus =
              eos.BulkModulusFromDensityInternalEnergy(density, energy);
          const Real ss_approx = std::sqrt(bmod_approx / density);
          const Real sound_speed = std::sqrt(bulk_modulus / density);
          INFO("Density: " << density << "  Energy: " << energy
                           << "  Sound speed: " << sound_speed
                           << " cm/us  Approximate sound speed: " << ss_approx
                           << " cm/us");
          REQUIRE(isClose(sound_speed, ss_approx, 1e-5));
        }
      }
    }
    GIVEN("A particle velocity and the same Hugoniot fit used in the EOS") {
      // Use Rankine-Hugoniot jump conditions to calculate a consistent point on the EOS
      constexpr Real up = 2 * cm / us; // 20 km/s is a pretty strong shock
      constexpr Real Us = C0 + S1 * up;
      constexpr Real e0 = 0;
      THEN("We have the density, energy, and pressure at this point") {
        constexpr Real density = Us * rho0 / (Us - up);
        constexpr Real true_pres = P0 + rho0 * Us * up;
        constexpr Real energy =
            e0 + 1. / 2. * Us * up * (1 - rho0 / density) + P0 * (1 / rho0 - 1 / density);
        WHEN("A P(rho, e) lookup is performed for the Hugoniot energy and density") {
          const Real pres = eos.PressureFromDensityInternalEnergy(density, energy);
          THEN("The pressure should agree with that given by the jump conditions") {
            INFO("Density: " << density << "  Energy: " << energy
                             << "  Pressure: " << pres / Mbar
                             << " Mbar  True pressure: " << true_pres / Mbar << " Mbar");
            REQUIRE(isClose(pres, true_pres, 1e-12));
          }
        }
      }
    }
  }
}

SCENARIO("Gruneisen EOS density limit") {
  /* These tests all test the functionality that finds the roots of the polynomial that
     makes up the denominator of the reference pressure curve, i.e.
      P(x) = 1 - s1 * x - s2 * x**2 - s3 * x**3,
     in order to find the maximum compression allowed by the EOS.

     It's very probable that many of these cases may not be physical because they produce
     unstable shocks (i.e. Us - up < c) or cannot propagate shocks (Us - up < 0) but we 
     admit these cases for the moment since a more comprehensive set of bounds on the EOS
     parameters is not currently available. For a linear Us-up relationship, the point at
     which up > Us occurs also corresponds to the singularity in the reference pressure.
  */
  GIVEN("Parameters for a Gruneisen EOS") {
    // Unit conversions
    constexpr Real mm = 10.;
    constexpr Real cm = 1.;
    constexpr Real us = 1.e-06;
    constexpr Real Mbar = 1.e12;
    constexpr Real Mbcc_per_g = 1e12;
    // Gruneisen parameters for copper
    constexpr Real C0 = 0.535 * cm / us;
    constexpr Real Gamma0 = 1.97;
    constexpr Real b = 0.;
    constexpr Real rho0 = 2.714000;
    constexpr Real T0 = 298.;
    constexpr Real P0 = 1e-06 * Mbar;
    constexpr Real Cv = 0.383e-05 * Mbcc_per_g;
    // No density limit
    constexpr Real no_rho_max = std::numeric_limits<Real>::infinity();
    WHEN("A small rho_max parameter is provided") {
      constexpr Real set_rho_max = 1.1 * rho0;
      // Linear Hugoniot fit
      constexpr Real S1 = 1.34;
      constexpr Real S2 = 0.;
      constexpr Real S3 = 0.;
      // Create the EOS
      Gruneisen host_eos = Gruneisen{C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv};
      auto eos = host_eos.GetOnDevice();
      THEN("The provided rho_max parameter should be less than the calculated rho_max") {
        const Real rho_max = eos.ComputeRhoMax(S1, S2, S3, rho0);
        INFO("Provided rho_max: " << set_rho_max << ", Calculated rho_max:" << rho_max);
        REQUIRE(rho_max > set_rho_max);
      }
    }
    WHEN("The rho_max parameter is not specified") {
      WHEN("A linear Hugoniot fit is used") {
        // Linear Hugoniot fit
        constexpr Real S1 = 1.34;
        constexpr Real S2 = 0.;
        constexpr Real S3 = 0.;
        // Create the EOS
        Gruneisen host_eos = Gruneisen{C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv};
        auto eos = host_eos.GetOnDevice();
        THEN("The generated rho_max parameter should be properly set") {
          constexpr Real eta_max = 1 / S1;
          const Real rho_max_true = rho0 / (1 - eta_max);
          const Real rho_max = eos.ComputeRhoMax(S1, S2, S3, rho0);
          INFO("True rho_max: " << rho_max_true << ", Calculated rho_max:" << rho_max);
          REQUIRE(isClose(rho_max, rho_max_true, 1e-12));
        }
      }
      WHEN("A quadratic Hugoniot fit is used") {
        WHEN("The fit is simple") {
          // Quadratic Hugoniot fit
          constexpr Real S1 = 2.0;
          constexpr Real S2 = 0.1;
          constexpr Real S3 = 0.;
          // Create the EOS
          Gruneisen host_eos = Gruneisen{C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv};
          auto eos = host_eos.GetOnDevice();
          THEN("The generated rho_max parameter should be properly set") {
            const Real eta_max = QuadFormulaMinus(-S2, -S1, 1.);
            const Real rho_max_true = rho0 / (1 - eta_max);
            const Real rho_max = eos.ComputeRhoMax(S1, S2, S3, rho0);
            INFO("True rho_max: " << rho_max_true << ", Calculated rho_max:" << rho_max);
            REQUIRE(isClose(rho_max, rho_max_true, 1e-12));
          }
        }
        WHEN("The quadratic oot multiplicity is 2") {
          // Quadratic Hugoniot fit
          constexpr Real S1 = 3;
          constexpr Real S2 = -2.25;
          constexpr Real S3 = 0.;
          // Create the EOS
          Gruneisen host_eos = Gruneisen{C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv};
          auto eos = host_eos.GetOnDevice();
          THEN("The generated rho_max parameter should be properly set") {
            const Real eta_max = QuadFormulaMinus(-S2, -S1, 1.);
            const Real rho_max_true = rho0 / (1 - eta_max);
            const Real rho_max = eos.ComputeRhoMax(S1, S2, S3, rho0);
            INFO("True rho_max: " << rho_max_true << ", Calculated rho_max:" << rho_max);
            REQUIRE(isClose(rho_max, rho_max_true, 1e-12));
          }
        }
        WHEN("No root exists") {
          // Quadratic Hugoniot fit without a pressure singularity
          constexpr Real S1 = 2;
          constexpr Real S2 = -1.1;
          constexpr Real S3 = 0.;
          // Create the EOS
          Gruneisen host_eos = Gruneisen{C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv};
          auto eos = host_eos.GetOnDevice();
          THEN("The generated rho_max parameter should be properly set") {
            constexpr Real rho_max_true = no_rho_max; // No maximum (see source)
            const Real rho_max = eos.ComputeRhoMax(S1, S2, S3, rho0);
            INFO("True rho_max: " << rho_max_true << ", Calculated rho_max:" << rho_max);
            REQUIRE(rho_max == rho_max_true);
          }
        }
        WHEN("The root is out of bounds") {
          // Quadratic Hugoniot fit without a pressure singularity
          // THIS ISN'T REAL... it might be better to specify what legal values
          // for the parameters are...
          constexpr Real S1 = 0.1;
          constexpr Real S2 = 0.1;
          constexpr Real S3 = 0.;
          // Create the EOS
          Gruneisen host_eos = Gruneisen{C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv};
          auto eos = host_eos.GetOnDevice();
          THEN("The generated rho_max parameter should be properly set") {
            constexpr Real rho_max_true = no_rho_max; // No maximum (see source)
            const Real rho_max = eos.ComputeRhoMax(S1, S2, S3, rho0);
            INFO("True rho_max: " << rho_max_true << ", Calculated rho_max:" << rho_max);
            REQUIRE(rho_max == rho_max_true);
          }
        }
      }
      WHEN("A cubic Hugoniot fit is used") {
        WHEN("Only one root exists") {
          // Cubic Hugoniot fit
          constexpr Real S1 = 2;
          constexpr Real S2 = -1.1;
          constexpr Real S3 = 0.2;
          // Create the EOS
          Gruneisen host_eos = Gruneisen{C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv};
          auto eos = host_eos.GetOnDevice();
          THEN("The generated rho_max parameter should be properly set") {
            constexpr Real eta_max = 0.8025706630669670284; // Wolfram alpha
            constexpr Real rho_max_true = rho0 / (1 - eta_max);
            const Real rho_max = eos.ComputeRhoMax(S1, S2, S3, rho0);
            INFO("True rho_max: " << rho_max_true << ", Calculated rho_max:" << rho_max);
            REQUIRE(isClose(rho_max, rho_max_true, 1e-8)); // root-find tolerance
          }
        }
        WHEN("A single root of multiplicity 3 exists") {
          // Cubic Hugoniot fit with a single root at 0.5 with multiplicity 3
          constexpr Real S1 = 6;
          constexpr Real S2 = -12;
          constexpr Real S3 = 8;
          // Create the EOS
          Gruneisen host_eos = Gruneisen{C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv};
          auto eos = host_eos.GetOnDevice();
          THEN("The generated rho_max parameter should be properly set") {
            constexpr Real eta_max = 0.5;
            constexpr Real rho_max_true = rho0 / (1 - eta_max);
            const Real rho_max = eos.ComputeRhoMax(S1, S2, S3, rho0);
            INFO("True rho_max: " << rho_max_true << ", Calculated rho_max:" << rho_max);
            REQUIRE(isClose(rho_max, rho_max_true, 1e-8));
          }
        }
        WHEN("Two roots exist") {
          // Cubic Hugoniot fit with roots at 0.5 (multiplicy 2) and -1
          constexpr Real S1 = 3;
          constexpr Real S2 = 0;
          constexpr Real S3 = -4;
          // Create the EOS
          Gruneisen host_eos = Gruneisen{C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv};
          auto eos = host_eos.GetOnDevice();
          THEN("The generated rho_max parameter should be properly set") {
            constexpr Real eta_max = 0.5;
            constexpr Real rho_max_true = rho0 / (1 - eta_max);
            const Real rho_max = eos.ComputeRhoMax(S1, S2, S3, rho0);
            INFO("True rho_max: " << rho_max_true << ", Calculated rho_max:" << rho_max);
            REQUIRE(isClose(rho_max, rho_max_true, 1e-8));
          }
        }
        WHEN("The cubic is decreasing and two roots exist") {
          // Cubic Hugoniot fit with roots at 0.8 and -0.4 (multiplicity 2)
          constexpr Real S1 = -3.75;
          constexpr Real S2 = 0;
          constexpr Real S3 = 7.8125;
          // Create the EOS
          Gruneisen host_eos = Gruneisen{C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv};
          auto eos = host_eos.GetOnDevice();
          THEN("The generated rho_max parameter should be properly set") {
            constexpr Real eta_max = 0.8;
            constexpr Real rho_max_true = rho0 / (1 - eta_max);
            const Real rho_max = eos.ComputeRhoMax(S1, S2, S3, rho0);
            INFO("True rho_max: " << rho_max_true << ", Calculated rho_max:" << rho_max);
            REQUIRE(isClose(rho_max, rho_max_true, 1e-8));
          }
        }
        WHEN("Three negative roots exist") {
          // Cubic Hugoniot fit with roots at -0.5, -1, and -2
          constexpr Real S1 = -3.5;
          constexpr Real S2 = -3.5;
          constexpr Real S3 = -1;
          // Create the EOS
          Gruneisen host_eos = Gruneisen{C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv};
          auto eos = host_eos.GetOnDevice();
          THEN("The generated rho_max parameter should be properly set (Cubic three "
               "negative roots)") {
            constexpr Real rho_max_true = no_rho_max; // No maximum (see source)
            const Real rho_max = eos.ComputeRhoMax(S1, S2, S3, rho0);
            INFO("True rho_max: " << rho_max_true << ", Calculated rho_max:" << rho_max);
            REQUIRE(rho_max == rho_max_true);
          }
        }
        WHEN("One root of three is bounded") {
          // Cubic Hugoniot fit with roots at 0.5, -1, and 2
          constexpr Real S1 = 2.5;
          constexpr Real S2 = -0.5;
          constexpr Real S3 = -1;
          // Create the EOS
          Gruneisen host_eos = Gruneisen{C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv};
          auto eos = host_eos.GetOnDevice();
          THEN("The generated rho_max parameter should be properly set (Cubic one "
               "bounded root)") {
            constexpr Real eta_max = 0.5;
            constexpr Real rho_max_true = rho0 / (1 - eta_max);
            const Real rho_max = eos.ComputeRhoMax(S1, S2, S3, rho0);
            INFO("True rho_max: " << rho_max_true << ", Calculated rho_max:" << rho_max);
            REQUIRE(isClose(rho_max, rho_max_true, 1e-8));
          }
        }
        WHEN("Two roots of three are bounded") {
          // Cubic Hugoniot fit with roots at 0.5, 0.8, and -2.5
          constexpr Real S1 = 2.85;
          constexpr Real S2 = -1.2;
          constexpr Real S3 = -1;
          // Create the EOS
          Gruneisen host_eos = Gruneisen{C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv};
          auto eos = host_eos.GetOnDevice();
          THEN("The generated rho_max parameter should be properly set (Cubic two "
               "bounded roots)") {
            constexpr Real eta_max = 0.5;
            constexpr Real rho_max_true = rho0 / (1 - eta_max);
            const Real rho_max = eos.ComputeRhoMax(S1, S2, S3, rho0);
            INFO("True rho_max: " << rho_max_true << ", Calculated rho_max:" << rho_max);
            REQUIRE(isClose(rho_max, rho_max_true, 1e-8));
          }
        }
        WHEN("Two roots of three are bounded and the cubic is decreasing") {
          // Cubic Hugoniot fit with roots at 0.2, 0.8, and 6.25
          constexpr Real S1 = 3.65;
          constexpr Real S2 = -3.8;
          constexpr Real S3 = 1;
          // Create the EOS
          Gruneisen host_eos = Gruneisen{C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv};
          auto eos = host_eos.GetOnDevice();
          THEN("The generated rho_max parameter should be properly set") {
            constexpr Real eta_max = 0.5;
            constexpr Real rho_max_true = rho0 / (1 - eta_max);
            const Real rho_max = eos.ComputeRhoMax(S1, S2, S3, rho0);
            INFO("True rho_max: " << rho_max_true << ", Calculated rho_max:" << rho_max);
            REQUIRE(isClose(rho_max, rho_max_true, 1e-8));
          }
        }
        WHEN("The three roots are all beyond 1") {
          // Cubic Hugoniot fit with roots at 1.25, 1.5, and 1.6
          constexpr Real S1 = 2.0916666667;
          constexpr Real S2 = -1.45;
          constexpr Real S3 = 1. / 3.;
          // Create the EOS
          Gruneisen host_eos = Gruneisen{C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv};
          auto eos = host_eos.GetOnDevice();
          THEN("The generated rho_max parameter should be properly set") {
            constexpr Real rho_max_true = no_rho_max; // No maximum (see source)
            const Real rho_max = eos.ComputeRhoMax(S1, S2, S3, rho0);
            INFO("True rho_max: " << rho_max_true << ", Calculated rho_max:" << rho_max);
            REQUIRE(rho_max == rho_max_true);
          }
        }
      }
    }
  }
}
