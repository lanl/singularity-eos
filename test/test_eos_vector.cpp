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

#include <array>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/eos/eos.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include "catch2/catch.hpp"
#endif

#include <test/eos_unit_test_helpers.hpp>

using singularity::EOS;
using singularity::IdealGas;

SCENARIO("Vector EOS", "[VectorEOS][IdealGas]") {

  GIVEN("Parameters for an ideal gas") {
    // Create ideal gas EOS ojbect
    constexpr Real Cv = 5.0;
    constexpr Real gm1 = 0.4;
    EOS host_eos = IdealGas(gm1, Cv);
    EOS eos = host_eos.GetOnDevice();

    GIVEN("Energies and densities") {
      constexpr int num = 3;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create Kokkos views on device for the input arrays
      Kokkos::View<Real[num]> v_density("density");
      Kokkos::View<Real[num]> v_energy("density");
#else
      // Otherwise just create arrays to contain values and create pointers to
      // be passed to the functions in place of the Kokkos views
      std::array<Real, num> density;
      std::array<Real, num> energy;
      Real *v_density = density.data();
      Real *v_energy = energy.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

      // Populate the input arrays
      portableFor(
          "Initialize density and energy", 0, 1, PORTABLE_LAMBDA(int i) {
            v_density[0] = 1.0;
            v_density[1] = 2.0;
            v_density[2] = 5.0;
            v_energy[0] = 5.0;
            v_energy[1] = 10.0;
            v_energy[2] = 15.0;
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create host-side mirrors of the inputs and copy the inputs. These are
      // just used for the comparisons
      auto density = Kokkos::create_mirror_view(v_density);
      auto energy = Kokkos::create_mirror_view(v_energy);
      Kokkos::deep_copy(density, v_density);
      Kokkos::deep_copy(energy, v_energy);
#endif // PORTABILITY_STRATEGY_KOKKOS

      // Gold standard values
      constexpr std::array<Real, num> pressure_true{2.0, 8.0, 30.0};
      constexpr std::array<Real, num> temperature_true{1., 2., 3.};
      constexpr std::array<Real, num> bulkmodulus_true{2.8, 11.2, 42.};
      constexpr std::array<Real, num> heatcapacity_true{Cv, Cv, Cv};
      constexpr std::array<Real, num> gruneisen_true{gm1, gm1, gm1};

      // Gold standard entropy doesn't produce round numbers so we need to
      // calculate it from the device views so this requires a bit more work
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<Real[num]> v_entropy_true("True Entropy");
#else
      std::array<Real, num> entropy_true;
      Real *v_entropy_true = entropy_true.data();
#endif
      constexpr Real P0 = 1e6;                    // microbar
      constexpr Real T0 = 293;                    // K
      constexpr Real rho0 = P0 / (gm1 * Cv * T0); // g/cm^3
      portableFor(
          "Calculate true entropy", 0, num, PORTABLE_LAMBDA(const int i) {
            v_entropy_true[i] =
                Cv * log(v_energy[i] / Cv / T0) + gm1 * Cv * log(rho0 / v_density[i]);
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      auto entropy_true = Kokkos::create_mirror_view(v_entropy_true);
      Kokkos::deep_copy(entropy_true, v_entropy_true);
#endif

#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create device views for outputs and mirror those views on the host
      Kokkos::View<Real[num]> v_temperature("Temperature");
      Kokkos::View<Real[num]> v_pressure("Pressure");
      Kokkos::View<Real[num]> v_entropy("Entropy");
      Kokkos::View<Real[num]> v_heatcapacity("Cv");
      Kokkos::View<Real[num]> v_bulkmodulus("bmod");
      Kokkos::View<Real[num]> v_gruneisen("Gamma");
      auto h_temperature = Kokkos::create_mirror_view(v_temperature);
      auto h_pressure = Kokkos::create_mirror_view(v_pressure);
      auto h_entropy = Kokkos::create_mirror_view(v_entropy);
      auto h_heatcapacity = Kokkos::create_mirror_view(v_heatcapacity);
      auto h_bulkmodulus = Kokkos::create_mirror_view(v_bulkmodulus);
      auto h_gruneisen = Kokkos::create_mirror_view(v_gruneisen);
#else
      // Create arrays for the outputs and then pointers to those arrays that
      // will be passed to the functions in place of the Kokkos views
      std::array<Real, num> h_temperature;
      std::array<Real, num> h_pressure;
      std::array<Real, num> h_entropy;
      std::array<Real, num> h_heatcapacity;
      std::array<Real, num> h_bulkmodulus;
      std::array<Real, num> h_gruneisen;
      // Just alias the existing pointers
      auto v_temperature = h_temperature.data();
      auto v_pressure = h_pressure.data();
      auto v_entropy = h_entropy.data();
      auto v_heatcapacity = h_heatcapacity.data();
      auto v_bulkmodulus = h_bulkmodulus.data();
      auto v_gruneisen = h_gruneisen.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

      WHEN("A T(rho, e) lookup is performed") {
        eos.TemperatureFromDensityInternalEnergy(v_density, v_energy, v_temperature, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_temperature, v_temperature);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned T(rho, e) should be equal to the true "
             "temperature") {
          array_compare(num, density, energy, h_temperature, temperature_true, "Density",
                        "Energy");
        }
      }

      WHEN("A P(rho, e) lookup is performed") {
        eos.PressureFromDensityInternalEnergy(v_density, v_energy, v_pressure, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_pressure, v_pressure);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned P(rho, e) should be equal to the true pressure") {
          array_compare(num, density, energy, h_pressure, pressure_true, "Density",
                        "Energy");
        }
      }

      WHEN("An S(rho, e) lookup is performed") {
        eos.EntropyFromDensityInternalEnergy(v_density, v_energy, v_entropy, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_entropy, v_entropy);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned S(rho, e) should be equal to the true entropy") {
          array_compare(num, density, energy, h_entropy, entropy_true, "Density",
                        "Energy");
        }
      }

      WHEN("A C_v(rho, e) lookup is performed") {
        eos.SpecificHeatFromDensityInternalEnergy(v_density, v_energy, v_heatcapacity,
                                                  num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_heatcapacity, v_heatcapacity);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned C_v(rho, e) should be constant") {
          array_compare(num, density, energy, h_heatcapacity, heatcapacity_true,
                        "Density", "Energy");
        }
      }

      WHEN("A B_S(rho, e) lookup is performed") {
        eos.BulkModulusFromDensityInternalEnergy(v_density, v_energy, v_bulkmodulus, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_bulkmodulus, v_bulkmodulus);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned B_S(rho, e) should be equal to the true bulk "
             "modulus") {
          array_compare(num, density, energy, h_bulkmodulus, bulkmodulus_true, "Density",
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
        THEN("The returned Gamma(rho, e) should be constant") {
          array_compare(num, density, energy, h_gruneisen, gruneisen_true, "Density",
                        "Energy");
        }
      }
    }
    GIVEN("Densities and temperatures") {
      constexpr int num = 3;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create Kokkos views on device for the input arrays
      Kokkos::View<Real[num]> v_density("density");
      Kokkos::View<Real[num]> v_temperature("density");
#else
      // Otherwise just create arrays to contain values and create pointers to
      // be passed to the functions in place of the Kokkos views
      std::array<Real, num> density;
      std::array<Real, num> temperature;
      Real *v_density = density.data();
      Real *v_temperature = temperature.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

      // Populate the input arrays
      portableFor(
          "Initialize density and energy", 0, 1, PORTABLE_LAMBDA(int i) {
            v_density[0] = 1.0;
            v_density[1] = 2.0;
            v_density[2] = 5.0;
            v_temperature[0] = 50.0;
            v_temperature[1] = 100.0;
            v_temperature[2] = 150.0;
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create host-side mirrors of the inputs and copy the inputs. These are
      // just used for the comparisons
      auto density = Kokkos::create_mirror_view(v_density);
      auto temperature = Kokkos::create_mirror_view(v_temperature);
      Kokkos::deep_copy(density, v_density);
      Kokkos::deep_copy(temperature, v_temperature);
#endif // PORTABILITY_STRATEGY_KOKKOS

      // Gold standard values
      constexpr std::array<Real, num> energy_true{250., 500., 750.};
      constexpr std::array<Real, num> pressure_true{100., 400., 1500.};
      constexpr std::array<Real, num> bulkmodulus_true{140., 560., 2100.};
      constexpr std::array<Real, num> heatcapacity_true{Cv, Cv, Cv};
      constexpr std::array<Real, num> gruneisen_true{gm1, gm1, gm1};

      // Gold standard entropy doesn't produce round numbers so we need to
      // calculate it from the device views so this requires a bit more work
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<Real[num]> v_entropy_true("True Entropy");
#else
      std::array<Real, num> entropy_true;
      Real *v_entropy_true = entropy_true.data();
#endif
      constexpr Real P0 = 1e6;                    // microbar
      constexpr Real T0 = 293;                    // K
      constexpr Real rho0 = P0 / (gm1 * Cv * T0); // g/cm^3
      portableFor(
          "Calculate true entropy", 0, num, PORTABLE_LAMBDA(const int i) {
            v_entropy_true[i] =
                Cv * log(v_temperature[i] / T0) + gm1 * Cv * log(rho0 / v_density[i]);
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      auto entropy_true = Kokkos::create_mirror_view(v_entropy_true);
      Kokkos::deep_copy(entropy_true, v_entropy_true);
#endif

#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create device views for outputs and mirror those views on the host
      Kokkos::View<Real[num]> v_energy("Energy");
      Kokkos::View<Real[num]> v_pressure("Pressure");
      Kokkos::View<Real[num]> v_entropy("Entropy");
      Kokkos::View<Real[num]> v_heatcapacity("Cv");
      Kokkos::View<Real[num]> v_bulkmodulus("bmod");
      Kokkos::View<Real[num]> v_gruneisen("Gamma");
      auto h_energy = Kokkos::create_mirror_view(v_energy);
      auto h_pressure = Kokkos::create_mirror_view(v_pressure);
      auto h_entropy = Kokkos::create_mirror_view(v_entropy);
      auto h_heatcapacity = Kokkos::create_mirror_view(v_heatcapacity);
      auto h_bulkmodulus = Kokkos::create_mirror_view(v_bulkmodulus);
      auto h_gruneisen = Kokkos::create_mirror_view(v_gruneisen);
#else
      // Create arrays for the outputs and then pointers to those arrays that
      // will be passed to the functions in place of the Kokkos views
      std::array<Real, num> h_energy;
      std::array<Real, num> h_pressure;
      std::array<Real, num> h_entropy;
      std::array<Real, num> h_heatcapacity;
      std::array<Real, num> h_bulkmodulus;
      std::array<Real, num> h_gruneisen;
      // Just alias the existing pointers
      auto v_energy = h_energy.data();
      auto v_pressure = h_pressure.data();
      auto v_entropy = h_entropy.data();
      auto v_heatcapacity = h_heatcapacity.data();
      auto v_bulkmodulus = h_bulkmodulus.data();
      auto v_gruneisen = h_gruneisen.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

      WHEN("A e(rho, T) lookup is performed") {
        eos.InternalEnergyFromDensityTemperature(v_density, v_temperature, v_energy, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_energy, v_energy);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned e(rho, T) should be equal to the true energy") {
          array_compare(num, density, temperature, h_energy, energy_true, "Density",
                        "Temperature");
        }
      }

      WHEN("A P(rho, T) lookup is performed") {
        eos.PressureFromDensityTemperature(v_density, v_temperature, v_pressure, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_pressure, v_pressure);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned P(rho, T) should be equal to the true pressure") {
          array_compare(num, density, temperature, h_pressure, pressure_true, "Density",
                        "Temperature");
        }
      }

      WHEN("An S(rho, T) lookup is performed") {
        eos.EntropyFromDensityTemperature(v_density, v_temperature, v_entropy, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_entropy, v_entropy);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned S(rho, T) should be equal to the true entropy") {
          array_compare(num, density, temperature, h_entropy, entropy_true, "Density",
                        "Temperature");
        }
      }

      WHEN("A C_v(rho, T) lookup is performed") {
        eos.SpecificHeatFromDensityTemperature(v_density, v_temperature, v_heatcapacity,
                                               num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_heatcapacity, v_heatcapacity);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned C_v(rho, T) should be constant") {
          array_compare(num, density, temperature, h_heatcapacity, heatcapacity_true,
                        "Density", "Temperature");
        }
      }

      WHEN("A B_S(rho, T) lookup is performed") {
        eos.BulkModulusFromDensityTemperature(v_density, v_temperature, v_bulkmodulus,
                                              num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_bulkmodulus, v_bulkmodulus);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned B_S(rho, T) should be equal to the true bulk "
             "modulus") {
          array_compare(num, density, temperature, h_bulkmodulus, bulkmodulus_true,
                        "Density", "Temperature");
        }
      }

      WHEN("A Gamma(rho, T) lookup is performed") {
        eos.GruneisenParamFromDensityTemperature(v_density, v_temperature, v_gruneisen,
                                                 num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_gruneisen, v_gruneisen);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned Gamma(rho, T) should be constant") {
          array_compare(num, density, temperature, h_gruneisen, gruneisen_true, "Density",
                        "Temperature");
        }
      }
    }
  }
}
