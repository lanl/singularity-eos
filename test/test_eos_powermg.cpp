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

using singularity::PowerMG;
using EOS = singularity::Variant<PowerMG>;

SCENARIO("PowerMG EOS rho sie", "[VectorEOS][PowerMGEOS]") {
  GIVEN("Parameters for a PowerMG EOS") {
    // Unit conversions
    constexpr Real Mbcc_per_g = 1e12;
    // PowerMG parameters for tin beta phase
    constexpr Real rho0 = 7.285;
    constexpr Real T0 = 298.0;
    constexpr Real Cs = 2766.0e2;
    constexpr Real s = 1.5344;
    constexpr Real G0 = 2.4659;
    constexpr Real Cv0 = 0.2149e-05 * Mbcc_per_g;
    constexpr Real E0 = 0.658e-03 * Mbcc_per_g;
    constexpr Real S0 = 0.4419e-05 * Mbcc_per_g;
    constexpr Real Pmin = -0.001;
    Real K0toK40[41] = {Cs * Cs * rho0,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.};

    for (int ind = 1; ind <= 20; ind++) {
      K0toK40[ind] = (ind + 1) * pow(s, ind);
    }
    // Create the EOS
    EOS host_eos = PowerMG(rho0, T0, G0, Cv0, E0, S0, Pmin, K0toK40);
    EOS eos = host_eos.GetOnDevice();
    PowerMG host_eos2 = PowerMG(rho0, T0, G0, Cv0, E0, S0, Pmin, K0toK40);
    PowerMG eos2 = host_eos2.GetOnDevice();

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
            v_density[0] = 7.0;
            v_density[1] = 8.0;
            v_density[2] = 9.0;
            v_density[3] = 7.285;
            v_energy[0] = 4.6e8;
            v_energy[1] = 1.146e9;
            v_energy[2] = 3.28e9;
            v_energy[3] = 0.658e9;
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
          1.77398550710815783e+02, 4.26645103516999e+02, 7.08167099610162836e+02, T0};
      constexpr std::array<Real, num> hugtemperature_true{
          2.69534427862980863e+02, 3.90542065829047e+02, 7.78960970661031411e+02, T0};
      constexpr std::array<Real, num> pressure_true{
          -3.55688813700104904e+09, 6.82999362849669e+10, 2.09379236075995941e+11, 0.};
      constexpr std::array<Real, num> hugpressure_true{-1.00000e-3, 6.690618533331688e+10,
                                                       2.12112220098032013e+11, 0.};
      constexpr std::array<Real, num> hugenergy_true{E0, 1.068414572008593e+09,
                                                     3.43213602888831663e+09, E0};
      constexpr std::array<Real, num> entropy_true{
          3.52007925459420029e+06, 4.7165703738323096e+06, 5.26934995463277772e+06, S0};
      constexpr std::array<Real, num> cv_true{Cv0, Cv0, Cv0, Cv0};
      constexpr std::array<Real, num> bulkmodulus_true{
          -9.12803262563614082e+09, 8.776332438691490e+11, 1.46522209717692139e+12,
          Cs * Cs * rho0};
      constexpr std::array<Real, num> gruneisen_true{G0 * rho0 / 7.0, G0 * rho0 / 8.0,
                                                     G0 * rho0 / 9.0, G0};

#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create device views for outputs and mirror those views on the host
      Kokkos::View<Real[num]> v_temperature("Temperature");
      Kokkos::View<Real[num]> v_hugtemperature("HugTemperature");
      Kokkos::View<Real[num]> v_pressure("Pressure");
      Kokkos::View<Real[num]> v_hugpressure("HugPressure");
      Kokkos::View<Real[num]> v_hugenergy("HugEnergy");
      Kokkos::View<Real[num]> v_entropy("Entropy");
      Kokkos::View<Real[num]> v_cv("Cv");
      Kokkos::View<Real[num]> v_bulkmodulus("bmod");
      Kokkos::View<Real[num]> v_gruneisen("Gamma");
      auto h_temperature = Kokkos::create_mirror_view(v_temperature);
      auto h_hugtemperature = Kokkos::create_mirror_view(v_hugtemperature);
      auto h_pressure = Kokkos::create_mirror_view(v_pressure);
      auto h_hugpressure = Kokkos::create_mirror_view(v_hugpressure);
      auto h_hugenergy = Kokkos::create_mirror_view(v_hugenergy);
      auto h_entropy = Kokkos::create_mirror_view(v_entropy);
      auto h_cv = Kokkos::create_mirror_view(v_cv);
      auto h_bulkmodulus = Kokkos::create_mirror_view(v_bulkmodulus);
      auto h_gruneisen = Kokkos::create_mirror_view(v_gruneisen);
#else
      // Create arrays for the outputs and then pointers to those arrays that
      // will be passed to the functions in place of the Kokkos views
      std::array<Real, num> h_temperature;
      std::array<Real, num> h_hugtemperature;
      std::array<Real, num> h_pressure;
      std::array<Real, num> h_hugpressure;
      std::array<Real, num> h_hugenergy;
      std::array<Real, num> h_entropy;
      std::array<Real, num> h_cv;
      std::array<Real, num> h_bulkmodulus;
      std::array<Real, num> h_gruneisen;
      // Just alias the existing pointers
      auto v_temperature = h_temperature.data();
      auto v_hugtemperature = h_hugtemperature.data();
      auto v_pressure = h_pressure.data();
      auto v_hugpressure = h_hugpressure.data();
      auto v_hugenergy = h_hugenergy.data();
      auto v_entropy = h_entropy.data();
      auto v_cv = h_cv.data();
      auto v_bulkmodulus = h_bulkmodulus.data();
      auto v_gruneisen = h_gruneisen.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

      // temporary

      WHEN("A PH(rho) lookup is performed") {
        portableFor(
            "Test Hugoniot pressure", 0, num, PORTABLE_LAMBDA(const int i) {
              v_hugpressure[i] = eos2.AllHugPressureFromDensity(v_density[i]);
            });
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_hugpressure, v_hugpressure);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned PH(rho) should be equal to the true value") {
          array_compare(num, density, energy, h_hugpressure, hugpressure_true, "Density",
                        "Energy");
        }
      }

      WHEN("A EH(rho) lookup is performed") {
        portableFor(
            "Test Hugoniot internal energy", 0, num, PORTABLE_LAMBDA(const int i) {
              v_hugenergy[i] = eos2.AllHugInternalEnergyFromDensity(v_density[i]);
            });
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_hugenergy, v_hugenergy);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned EH(rho) should be equal to the true value") {
          array_compare(num, density, energy, h_hugenergy, hugenergy_true, "Density",
                        "Energy");
        }
      }

      WHEN("A TH(rho) lookup is performed") {
        portableFor(
            "Test Hugoniot temperature", 0, num, PORTABLE_LAMBDA(const int i) {
              v_hugtemperature[i] = eos2.AllHugTemperatureFromDensity(v_density[i]);
            });
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_hugtemperature, v_hugtemperature);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned TH(rho) should be equal to the true value") {
          array_compare(num, density, energy, h_hugtemperature, hugtemperature_true,
                        "Density", "Energy");
        }
      }
      // end temporary

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

      WHEN("a [rho,sie](P, T) lookup is performed") {
        int nwrong = 0;
        portableReduce(
            "Check density energy from pressure temperature", 0, 1,
            PORTABLE_LAMBDA(const int i, int &nw) {
              nw += !CheckRhoSieFromPT(eos, 7.285, 7.78960970661031411e+02);
            },
            nwrong);
        THEN("There are no errors") { REQUIRE(nwrong == 0); }
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
SCENARIO("PowerMG EOS rho T", "[VectorEOS][PowerMGEOS]") {
  GIVEN("Parameters for a PowerMG EOS") {
    // Unit conversions
    constexpr Real Mbcc_per_g = 1e12;
    // PowerMG parameters for copper
    constexpr Real rho0 = 7.285;
    constexpr Real T0 = 298.0;
    constexpr Real Cs = 2766.0e2;
    constexpr Real s = 1.5344;
    constexpr Real G0 = 2.4659;
    constexpr Real Cv0 = 0.2149e-05 * Mbcc_per_g;
    constexpr Real E0 = 0.658e-03 * Mbcc_per_g;
    constexpr Real S0 = 0.4419e-05 * Mbcc_per_g;
    constexpr Real Pmin = -0.001;
    Real K0toK40[41] = {Cs * Cs * rho0,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.};

    for (int ind = 1; ind <= 20; ind++) {
      K0toK40[ind] = (ind + 1) * pow(s, ind);
    }
    // Create the EOS
    EOS host_eos = PowerMG(rho0, T0, G0, Cv0, E0, S0, Pmin, K0toK40);
    EOS eos = host_eos.GetOnDevice();
    PowerMG host_eos2 = PowerMG(rho0, T0, G0, Cv0, E0, S0, Pmin, K0toK40);
    PowerMG eos2 = host_eos2.GetOnDevice();

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
            v_density[0] = 7.0;
            v_density[1] = 8.0;
            v_density[2] = 9.0;
            v_density[3] = 7.285;
            v_temperature[0] = 150.215;
            v_temperature[1] = 426.645;
            v_temperature[2] = 708.165;
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
      constexpr std::array<Real, num> sie_true{4.01582549522456884e+08,
                                               1.1459997775419691e+09,
                                               3.27999548793775988e+09, 6.58e+08};
      constexpr std::array<Real, num> pressure_true{
          -4.60630397840184784e+09, 6.82999322887127e+10, 2.09379155020942139e+11, 0.};
      constexpr std::array<Real, num> entropy_true{
          3.1626332929526418e+06, 4.7165698524198849e+06, 5.2693435831578374e+06, S0};
      constexpr std::array<Real, num> cv_true{Cv0, Cv0, Cv0, Cv0};
      constexpr std::array<Real, num> tbulkmodulus_true{
          -2.67031598038854675e+10, 8.4064844786200464e+11, 1.41065388996176587e+12,
          (Cs * Cs - G0 * G0 * T0 * Cv0) * rho0};
      constexpr std::array<Real, num> alpha_true{
          -1.44570198534642227e-03, 4.5922657969193220e-05, 2.73666073713845830e-05,
          G0 * Cv0 / (Cs * Cs - G0 * G0 * T0 * Cv0)};

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
SCENARIO("PowerMG EOS SetUp", "[VectorEOS][PowerMGEOS]") {
  GIVEN("Parameters for a PowerMG EOS") {
    // Unit conversions
    constexpr Real Mbcc_per_g = 1e12;
    // PowerMG parameters for copper
    Real rho0 = 7.285;
    Real T0 = 298.0;
    constexpr Real Cs = 2766.0e2;
    // constexpr Real s = 1.5344;
    constexpr Real G0 = 2.4659;
    Real Cv0 = 0.2149e-05 * Mbcc_per_g;
    constexpr Real E0 = 0.658e-03 * Mbcc_per_g;
    constexpr Real S0 = 0.4419e-05 * Mbcc_per_g;
    Real Pmin = -0.001;
    Real K0toK40[41] = {Cs * Cs * rho0,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.,
                        0.};

    EOS host_eos = PowerMG(rho0, T0, G0, Cv0, E0, S0, Pmin, K0toK40);
    EOS eos = host_eos.GetOnDevice();
    eos.PrintParams();
    eos.Finalize();

    WHEN("Faulty/not set parameter rho0 is given") {
      rho0 = -1.0;
      REQUIRE_MAYBE_THROWS(PowerMG(rho0, T0, G0, Cv0, E0, S0, Pmin, K0toK40));
      THEN("An error message should be written out") {}
    }
    WHEN("Faulty/not set parameter T0 is given") {
      T0 = -1.0;
      REQUIRE_MAYBE_THROWS(PowerMG(rho0, T0, G0, Cv0, E0, S0, Pmin, K0toK40));
      THEN("An error message should be written out") {}
    }
    WHEN("Faulty/not set parameter Cv0 is given") {
      Cv0 = -1.0;
      REQUIRE_MAYBE_THROWS(PowerMG(rho0, T0, G0, Cv0, E0, S0, Pmin, K0toK40));
      THEN("An error message should be written out") {}
    }
    WHEN("Faulty/not set parameter Pmin is given") {
      Real dens = 6.0;
      Real temp = 100.0;
      Pmin = 1.0;
      EOS host_eos1 = PowerMG(rho0, T0, G0, Cv0, E0, S0, Pmin, K0toK40);
      EOS eos1 = host_eos1.GetOnDevice();
      Real p1 = eos1.PressureFromDensityTemperature(dens, temp);
      //      Pmin = -0.01 * K0toK40[0]; should fail the test
      Pmin = -1000.0 * K0toK40[0];
      EOS host_eos2 = PowerMG(rho0, T0, G0, Cv0, E0, S0, Pmin, K0toK40);
      EOS eos2 = host_eos2.GetOnDevice();
      Real p2 = eos2.PressureFromDensityTemperature(dens, temp);
      THEN("A warning message should be written out and Pmin set to default") {
        INFO("Pmin set to default failed p1: " << p1 << " p2: " << p2 << "\n");
        CHECK(isClose(p1, p2, 1.0e-15));
      }
      eos1.Finalize();
      eos2.Finalize();
    }
    WHEN("Faulty/not set parameter K0 is given") {
      K0toK40[0] = -1.0;
      REQUIRE_MAYBE_THROWS(PowerMG(rho0, T0, G0, Cv0, E0, S0, Pmin, K0toK40));
      THEN("An error message should be written out") {}
    }
  }
}
