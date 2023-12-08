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
#include <cmath>
#include <cstdio>
#include <cstdlib>
#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <test/eos_unit_test_helpers.hpp>

using singularity::Vinet;
using EOS = singularity::Variant<Vinet>;

SCENARIO("Vinet EOS rho sie", "[VectorEOS][VinetEOS]") {
  GIVEN("Parameters for a Vinet EOS") {
    // Unit conversions
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
    constexpr Real d2to40[39] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    // Create the EOS
    EOS host_eos = Vinet(rho0, T0, B0, BP0, A0, Cv0, E0, S0, d2to40);
    EOS eos = host_eos.GetOnDevice();
    Vinet host_eos2 = Vinet(rho0, T0, B0, BP0, A0, Cv0, E0, S0, d2to40);
    Vinet eos2 = host_eos2.GetOnDevice();

    eos.PrintParams();

    GIVEN("Densities and energies") {
      constexpr int num = 4;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create Kokkos views on device for the input arrays
      Kokkos::View<Real[num]> v_density("density");
      // AEM: I think it should not be "density" again below
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
      constexpr std::array<Real, num> temperature_true{
          6.638294350641033e+1, 1.42266691572769e+03, 1.52867838544081e+03, T0};
      constexpr std::array<Real, num> pressure_true{
          -1.283459624233147e+11, 1.98538351697004e+10, 9.66446220551959e+10, 0.};
      constexpr std::array<Real, num> entropy_true{
          5.0015771521255749e+08, 5.1138262492866594e+08, 5.1120147992457777e+08, S0};
      constexpr std::array<Real, num> cv_true{Cv0, Cv0, Cv0, Cv0};
      constexpr std::array<Real, num> bulkmodulus_true{
          7.44411028807209e+11, 1.254880120063067e+12, 1.613822819346433e+12,
          (B0 + B0 * B0 * A0 * A0 / Cv0 / rho0 * T0)};
      constexpr std::array<Real, num> gruneisen_true{
          B0 * A0 / Cv0 / 8.0, B0 * A0 / Cv0 / 8.5, B0 * A0 / Cv0 / 9.0,
          B0 * A0 / Cv0 / rho0};

#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create device views for outputs and mirror those views on the host
      Kokkos::View<Real[num]> v_temperature("Temperature");
      Kokkos::View<Real[num]> v_pressure("Pressure");
      Kokkos::View<Real[num]> v_entropy("Entropy");
      Kokkos::View<Real[num]> v_cv("Cv");
      Kokkos::View<Real[num]> v_bulkmodulus("bmod");
      Kokkos::View<Real[num]> v_gruneisen("Gamma");
      auto h_temperature = Kokkos::create_mirror_view(v_temperature);
      auto h_pressure = Kokkos::create_mirror_view(v_pressure);
      auto h_entropy = Kokkos::create_mirror_view(v_entropy);
      auto h_cv = Kokkos::create_mirror_view(v_cv);
      auto h_bulkmodulus = Kokkos::create_mirror_view(v_bulkmodulus);
      auto h_gruneisen = Kokkos::create_mirror_view(v_gruneisen);
#else
      // Create arrays for the outputs and then pointers to those arrays that
      // will be passed to the functions in place of the Kokkos views
      std::array<Real, num> h_temperature;
      std::array<Real, num> h_pressure;
      std::array<Real, num> h_entropy;
      std::array<Real, num> h_cv;
      std::array<Real, num> h_bulkmodulus;
      std::array<Real, num> h_gruneisen;
      // Just alias the existing pointers
      auto v_temperature = h_temperature.data();
      auto v_pressure = h_pressure.data();
      auto v_entropy = h_entropy.data();
      auto v_cv = h_cv.data();
      auto v_bulkmodulus = h_bulkmodulus.data();
      auto v_gruneisen = h_gruneisen.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

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

      WHEN("A S(rho, e) lookup is performed") {
        portableFor(
            "Test entropy", 0, num, PORTABLE_LAMBDA(const int i) {
              v_entropy[i] =
                  eos2.EntropyFromDensityInternalEnergy(v_density[i], v_energy[i]);
            });
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_entropy, v_entropy);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned S(rho, e) should be equal to the true value") {
          array_compare(num, density, energy, h_entropy, entropy_true, "Density",
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

      WHEN("A Cv(rho, e) lookup is performed") {
        eos.SpecificHeatFromDensityInternalEnergy(v_density, v_energy, v_cv, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_cv, v_cv);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned T(rho, e) should be equal to the true value") {
          array_compare(num, density, energy, h_cv, cv_true, "Density", "Energy");
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
SCENARIO("Vinet EOS rho T", "[VectorEOS][VinetEOS]") {
  GIVEN("Parameters for a Vinet EOS") {
    // Unit conversions
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
    constexpr Real d2to40[39] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    // Create the EOS
    EOS host_eos = Vinet(rho0, T0, B0, BP0, A0, Cv0, E0, S0, d2to40);
    EOS eos = host_eos.GetOnDevice();
    Vinet host_eos2 = Vinet(rho0, T0, B0, BP0, A0, Cv0, E0, S0, d2to40);
    Vinet eos2 = host_eos2.GetOnDevice();

    eos.PrintParams();

    GIVEN("Densities and temperatures") {
      constexpr int num = 4;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create Kokkos views on device for the input arrays
      Kokkos::View<Real[num]> v_density("density");
      Kokkos::View<Real[num]> v_temperature("temperature");
#else
      // Otherwise just create arrays to contain values and create pointers to
      // be passed to the functions in place of the Kokkos views
      std::array<Real, num> density;
      std::array<Real, num> temperature;
      auto v_density = density.data();
      auto v_temperature = temperature.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

      // Populate the input arrays
      portableFor(
          "Initialize density and temperature", 0, 1, PORTABLE_LAMBDA(int i) {
            v_density[0] = 8.0;
            v_density[1] = 8.5;
            v_density[2] = 9.0;
            v_density[3] = 8.93;
            v_temperature[0] = 66.383;
            v_temperature[1] = 1422.7;
            v_temperature[2] = 1528.7;
            v_temperature[3] = 298.0;
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create host-side mirrors of the inputs and copy the inputs. These are
      // just used for the comparisons
      auto density = Kokkos::create_mirror_view(v_density);
      auto temperature = Kokkos::create_mirror_view(v_temperature);
      Kokkos::deep_copy(density, v_density);
      Kokkos::deep_copy(temperature, v_temperature);
#endif // PORTABILITY_STRATEGY_KOKKOS

      // Gold standard values for a subset of lookups
      constexpr std::array<Real, num> sie_true{
          2.0000021637044835e+8, 4.6001267127629471e+09, 4.7000827837617149e+09, 0.0};
      constexpr std::array<Real, num> pressure_true{
          -1.283459584783398e+11, 1.98561454605571e+10, 9.66461314103968e+10, 0.};
      constexpr std::array<Real, num> entropy_true{
          5.0015771847198445e+08, 5.1138271399469268e+08, 5.1120153407800680e+08, S0};
      constexpr std::array<Real, num> cv_true{Cv0, Cv0, Cv0, Cv0};
      constexpr std::array<Real, num> tbulkmodulus_true{
          7.3384631127398950e+11, 1.0417839336794381e+12, 1.3975684023296028e+12, B0};
      constexpr std::array<Real, num> alpha_true{
          9.5156828083623124e-05, 6.7029721830195901e-05, 4.9965702691402983e-05, A0};
      //    B0 * A0 / tbulkmodulus_true[0], B0 * A0 / tbulkmodulus_true[1], B0 * A0 /
      //    tbulkmodulus_true[2]

#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create device views for outputs and mirror those views on the host
      Kokkos::View<Real[num]> v_energy("Energy");
      Kokkos::View<Real[num]> v_pressure("Pressure");
      Kokkos::View<Real[num]> v_entropy("Entropy");
      Kokkos::View<Real[num]> v_cv("Cv");
      Kokkos::View<Real[num]> v_tbulkmodulus("tbmod");
      Kokkos::View<Real[num]> v_alpha("alpha");
      auto h_energy = Kokkos::create_mirror_view(v_energy);
      auto h_pressure = Kokkos::create_mirror_view(v_pressure);
      auto h_entropy = Kokkos::create_mirror_view(v_entropy);
      auto h_cv = Kokkos::create_mirror_view(v_cv);
      auto h_tbulkmodulus = Kokkos::create_mirror_view(v_tbulkmodulus);
      auto h_alpha = Kokkos::create_mirror_view(v_alpha);
#else
      // Create arrays for the outputs and then pointers to those arrays that
      // will be passed to the functions in place of the Kokkos views
      std::array<Real, num> h_energy;
      std::array<Real, num> h_pressure;
      std::array<Real, num> h_entropy;
      std::array<Real, num> h_cv;
      std::array<Real, num> h_tbulkmodulus;
      std::array<Real, num> h_alpha;
      // Just alias the existing pointers
      auto v_energy = h_energy.data();
      auto v_pressure = h_pressure.data();
      auto v_entropy = h_entropy.data();
      auto v_cv = h_cv.data();
      auto v_tbulkmodulus = h_tbulkmodulus.data();
      auto v_alpha = h_alpha.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

      WHEN("A E(rho, T) lookup is performed") {
        eos.InternalEnergyFromDensityTemperature(v_density, v_temperature, v_energy, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_energy, v_energy);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned E(rho, T) should be equal to the true value") {
          array_compare(num, density, temperature, h_energy, sie_true, "Density",
                        "Temperature");
        }
      }

      WHEN("A P(rho, T) lookup is performed") {
        eos.PressureFromDensityTemperature(v_density, v_temperature, v_pressure, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_pressure, v_pressure);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned P(rho, T) should be equal to the true value") {
          array_compare(num, density, temperature, h_pressure, pressure_true, "Density",
                        "Temperature");
        }
      }

      portableFor(
          "entropy", 0, num, PORTABLE_LAMBDA(const int i) {
            v_entropy[i] =
                eos2.EntropyFromDensityInternalEnergy(v_density[i], v_energy[i]);
          });
      WHEN("A S(rho, T) lookup is performed") {
        portableFor(
            "entropy", 0, num, PORTABLE_LAMBDA(const int i) {
              v_entropy[i] =
                  eos2.EntropyFromDensityTemperature(v_density[i], v_temperature[i]);
            });
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_entropy, v_entropy);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned S(rho, T) should be equal to the true value") {
          array_compare(num, density, temperature, h_entropy, entropy_true, "Density",
                        "Temperature");
        }
      }

      WHEN("A B_T(rho, T) lookup is performed") {
        portableFor(
            "bulk modulus", 0, num, PORTABLE_LAMBDA(const int i) {
              v_tbulkmodulus[i] =
                  eos2.TBulkModulusFromDensityTemperature(v_density[i], v_temperature[i]);
            });
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_tbulkmodulus, v_tbulkmodulus);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned B_T(rho, T) should be equal to the true value") {
          array_compare(num, density, temperature, h_tbulkmodulus, tbulkmodulus_true,
                        "Density", "Temperature");
        }
      }

      WHEN("A Cv(rho, T) lookup is performed") {
        eos.SpecificHeatFromDensityTemperature(v_density, v_temperature, v_cv, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_cv, v_cv);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned T(rho, e) should be equal to the true value") {
          array_compare(num, density, temperature, h_cv, cv_true, "Density",
                        "Temperature");
        }
      }

      WHEN("A alpha(rho, T) lookup is performed") {
        portableFor(
            "alpha", 0, num, PORTABLE_LAMBDA(const int i) {
              v_alpha[i] = eos2.TExpansionCoeffFromDensityTemperature(v_density[i],
                                                                      v_temperature[i]);
            });
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_alpha, v_alpha);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned alpha(rho, T) should be equal to the true value") {
          array_compare(num, density, temperature, h_alpha, alpha_true, "Density",
                        "Temperature");
        }
      }
    }
  }
}
SCENARIO("Vinet EOS SetUp", "[VectorEOS][VinetEOS]") {
  GIVEN("Parameters for a Vinet EOS") {
    // Unit conversions
    constexpr Real Mbcc_per_g = 1e12;
    // Vinet parameters for copper
    Real rho0 = 8.93;
    Real T0 = 298.0;
    Real B0 = 1.3448466 * Mbcc_per_g;
    Real BP0 = 4.956;
    Real A0 = 5.19245e-05;
    Real Cv0 = 0.383e-05 * Mbcc_per_g;
    constexpr Real E0 = 0.0;
    constexpr Real S0 = 5.05e-04 * Mbcc_per_g;
    constexpr Real d2to40[39] = {0., 0.,      0.0000000000001,
                                 0., 0.,      0.,
                                 0., 0.,      0.,
                                 0., 0.,      0.,
                                 0., 0.,      0.,
                                 0., 0.,      0.,
                                 0., 0.,      0.,
                                 0., 0.,      0.,
                                 0., 0.,      0.,
                                 0., 0.,      0.,
                                 0., 0.,      0.,
                                 0., 1.0e-26, 0.,
                                 0., 0.,      0.};
    // Create the EOS
    EOS host_eos = Vinet(rho0, T0, B0, BP0, A0, Cv0, E0, S0, d2to40);
    EOS eos = host_eos.GetOnDevice();
    eos.Finalize();

    WHEN("Faulty/not set parameter rho0 is given") {
      rho0 = -1.0;
      REQUIRE_MAYBE_THROWS(Vinet(rho0, T0, B0, BP0, A0, Cv0, E0, S0, d2to40));
      THEN("An error message should be written out") {}
    }
    WHEN("Faulty/not set parameter T0 is given") {
      T0 = -1.0;
      REQUIRE_MAYBE_THROWS(Vinet(rho0, T0, B0, BP0, A0, Cv0, E0, S0, d2to40));
      THEN("An error message should be written out") {}
    }
    WHEN("Faulty/not set parameter B0 is given") {
      B0 = -1.0;
      REQUIRE_MAYBE_THROWS(Vinet(rho0, T0, B0, BP0, A0, Cv0, E0, S0, d2to40));
      THEN("An error message should be written out") {}
    }
    WHEN("Faulty/not set parameter BP0 is given") {
      BP0 = 0.0;
      REQUIRE_MAYBE_THROWS(Vinet(rho0, T0, B0, BP0, A0, Cv0, E0, S0, d2to40));
      THEN("An error message should be written out") {}
    }
    WHEN("Faulty/not set parameter Cv0 is given") {
      Cv0 = -1.0;
      REQUIRE_MAYBE_THROWS(Vinet(rho0, T0, B0, BP0, A0, Cv0, E0, S0, d2to40));
      THEN("An error message should be written out") {}
    }
    WHEN("Faulty/not set parameter A0 is given") {
      A0 = -10000000.0;
      REQUIRE_MAYBE_THROWS(Vinet(rho0, T0, B0, BP0, A0, Cv0, E0, S0, d2to40));
      THEN("An error message should be written out") {}
    }
  }
}
