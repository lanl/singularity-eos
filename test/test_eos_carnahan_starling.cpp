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
#include <limits>
#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <test/eos_unit_test_helpers.hpp>

using singularity::CarnahanStarling;
using singularity::IdealGas;
using EOS = singularity::Variant<IdealGas, CarnahanStarling>;

SCENARIO("CarnahanStarling1", "[CarnahanStarling][CarnahanStarling1]") {
  GIVEN("Parameters for a CarnahanStarling EOS") {
    constexpr Real gm1 = 0.4;
    constexpr Real Cv = 7180000.0;
    constexpr Real bb = 1.e-3;
    constexpr Real qq = 0.0;
    // Create the EOS
    EOS host_eos = CarnahanStarling(gm1, Cv, bb, qq);
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
            v_density[0] = 1.1833012071291069e-03;
            v_density[1] = 7.8290736890381501e-03;
            v_density[2] = 8.5453943327882340e-03;
            v_density[3] = 8.8197619601121831e-03;
            v_energy[0] = 2.1407169999999998e+09;
            v_energy[1] = 1.1000478000000000e+10;
            v_energy[2] = 1.9860239000000000e+10;
            v_energy[3] = 2.8720000000000000e+10;
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
          1.0132499999999999e+06, 3.4450500000000000e+07, 6.7887750000000000e+07,
          1.0132500000000000e+08};
      constexpr std::array<Real, num> bulkmodulus_true{
          1.4185567142990597e+06, 4.8232210423710555e+07, 9.5046098754186615e+07,
          1.4186000457238689e+08};
      constexpr std::array<Real, num> temperature_true{
          2.9814999999999998e+02, 1.5320999999999999e+03, 2.7660500000000002e+03,
          4.0000000000000000e+03};
      constexpr std::array<Real, num> gruneisen_true{
          4.0000189328753222e-01, 4.0001252676308346e-01, 4.0001367292303214e-01,
          4.0001411193029396e-01};

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

SCENARIO("CarnahanStarling2", "[CarnahanStarling][CarnahanStarling2]") {
  GIVEN("Parameters for a CarnahanStarling EOS") {
    constexpr Real gm1 = 0.4;
    constexpr Real Cv = 7180000.0;
    constexpr Real bb = 1.e-3;
    constexpr Real qq = 42.00e+09;
    constexpr Real qp = -23.0e+7;
    constexpr Real T0 = 200.0;
    constexpr Real P0 = 1000000.0;
    // Create the EOS
    EOS host_eos = CarnahanStarling(gm1, Cv, bb, qq, qp, T0, P0);
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
            v_density[0] = 1.1833012071291069e-03;
            v_density[1] = 7.8290736890381501e-03;
            v_density[2] = 8.5453943327882340e-03;
            v_density[3] = 8.8197619601121831e-03;
            v_energy[0] = 4.4140717000000000e+10;
            v_energy[1] = 5.3000478000000000e+10;
            v_energy[2] = 6.1860239000000000e+10;
            v_energy[3] = 7.0720000000000000e+10;
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
          1.0132499999999999e+06, 3.4450500000000000e+07, 6.7887750000000000e+07,
          1.0132500000000000e+08};
      constexpr std::array<Real, num> bulkmodulus_true{
          1.4185567142990597e+06, 4.8232210423710555e+07, 9.5046098754186615e+07,
          1.4186000457238689e+08};
      constexpr std::array<Real, num> temperature_true{
          2.9814999999999998e+02, 1.5320999999999999e+03, 2.7660500000000002e+03,
          4.0000000000000000e+03};
      constexpr std::array<Real, num> gruneisen_true{
          4.0000189328753222e-01, 4.0001252676308346e-01, 4.0001367292303214e-01,
          4.0001411193029396e-01};

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

SCENARIO("(C-S EoS) Isentropic Bulk Modulus Analytic vs. FD",
         "[CarnahanStarling][CarnahanStarling3]") {
  GIVEN("Parameters for a C-S Gas EOS") {
    constexpr Real gm1 = 0.4;
    constexpr Real Cv = 7180000.0;
    constexpr Real bb = 1.e-3;
    constexpr Real qq = 0.0;
    //  Create the EOS
    EOS host_eos = CarnahanStarling(gm1, Cv, bb, qq);
    EOS eos = host_eos.GetOnDevice();
    GIVEN("Density and energy") {
      constexpr Real density = 1.1833012071291069e-03;
      constexpr Real energy = 2.1407169999999998e+09;
      constexpr Real true_sound_speed = 3.4623877142797290e+04;
      WHEN("A B_S(rho, e) lookup is performed") {
        const Real bulk_modulus =
            eos.BulkModulusFromDensityInternalEnergy(density, energy);
        THEN("The correct sound speed should be computed") {
          const Real sound_speed = std::sqrt(bulk_modulus / density);
          INFO("Density: " << density << "  Energy: " << energy
                           << "  Sound speed: " << sound_speed
                           << " cm/s  True sound speed: " << true_sound_speed << " cm/s");
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
                           << " cm/s  Approximate sound speed: " << ss_approx << " cm/s");
          REQUIRE(isClose(sound_speed, ss_approx, 1e-5));
        }
      }
    }
  }
}

SCENARIO("Recover Ideal Gas from C-S", "[CarnahanStarling][CarnahanStarling4]") {
  GIVEN("Parameters for a CarnahanStarling EOS") {
    constexpr Real gm1 = 0.4;
    constexpr Real Cv = 7180000.0;
    constexpr Real bb = 0;
    constexpr Real qq = 0;
    //  Create the EOS
    EOS host_eos = CarnahanStarling(gm1, Cv, bb, qq);
    EOS eos = host_eos.GetOnDevice();
    EOS ideal_eos = singularity::IdealGas(gm1, Cv);
    GIVEN("Densities and energies") {
      constexpr int num = 1;
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
            v_density[0] = 1.1833068079526625e-03;
            v_energy[0] = 2.1407169999999998e+09;
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create host-side mirrors of the inputs and copy the inputs. These are
      // just used for the comparisons
      auto density = Kokkos::create_mirror_view(v_density);
      auto energy = Kokkos::create_mirror_view(v_energy);
      Kokkos::deep_copy(density, v_density);
      Kokkos::deep_copy(energy, v_energy);
#endif // PORTABILITY_STRATEGY_KOKKOS

      // values for a subset of lookups
      std::array<Real, num> pressure_true;
      std::array<Real, num> bulkmodulus_true;
      std::array<Real, num> temperature_true;
      std::array<Real, num> gruneisen_true;

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
          for (int i = 0; i < num; i++) {
            pressure_true[i] =
                ideal_eos.PressureFromDensityInternalEnergy(density[i], energy[i]);
          }
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
          for (int i = 0; i < num; i++) {
            bulkmodulus_true[i] =
                ideal_eos.BulkModulusFromDensityInternalEnergy(density[i], energy[i]);
          }
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
          for (int i = 0; i < num; i++) {
            temperature_true[i] =
                ideal_eos.TemperatureFromDensityInternalEnergy(density[i], energy[i]);
          }
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
          for (int i = 0; i < num; i++) {
            gruneisen_true[i] =
                ideal_eos.GruneisenParamFromDensityInternalEnergy(density[i], energy[i]);
          }
          array_compare(num, density, energy, h_gruneisen, gruneisen_true, "Density",
                        "Energy");
        }
      }
    }
  }
}

SCENARIO("Test C-S Entropy Calls", "[CarnahanStarling][CarnahanStarling5]") {
  GIVEN("Parameters for a CarnahanStarling EOS") {
    constexpr Real gm1 = 0.4;
    constexpr Real Cv = 7180000.0;
    constexpr Real bb = 1.e-3;
    constexpr Real qq = 42.0e+9;
    constexpr Real qp = -23.0e+9;
    constexpr Real T0 = 200.0;
    constexpr Real P0 = 1000000.0;
    //  Create the EOS
    EOS host_eos = CarnahanStarling(gm1, Cv, bb, qq, qp, T0, P0);
    EOS eos = host_eos.GetOnDevice();
    GIVEN("Densities and energies") {
      constexpr int num = 1;
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
            v_density[0] = 8.8197619601121831e-03;
            v_energy[0] = 7.0720000000000000e+10;
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
      constexpr std::array<Real, num> temperature_true{4.0000000000000000e+03};
      constexpr std::array<Real, num> entropy_true{-2.2983150752058342e+10};

#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create device views for outputs and mirror those views on the host
      Kokkos::View<Real[num]> v_temperature("Temperature");
      Kokkos::View<Real[num]> v_pressure("Pressure");
      Kokkos::View<Real[num]> v_entropy("entr");
      Kokkos::View<Real[num]> v_local_temp("temp");
      auto h_temperature = Kokkos::create_mirror_view(v_temperature);
      auto h_pressure = Kokkos::create_mirror_view(v_pressure);
      auto h_entropy = Kokkos::create_mirror_view(v_entropy);
      auto h_local_temp = Kokkos::create_mirror_view(v_local_temp);
#else
      // Create arrays for the outputs and then pointers to those arrays that
      // will be passed to the functions in place of the Kokkos views
      std::array<Real, num> h_temperature;
      std::array<Real, num> h_pressure;
      std::array<Real, num> h_entropy;
      std::array<Real, num> h_local_temp;
      // Just alias the existing pointers
      auto v_temperature = h_temperature.data();
      auto v_pressure = h_pressure.data();
      auto v_entropy = h_entropy.data();
      auto v_local_temp = h_local_temp.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

      WHEN("A S(rho, e) lookup is performed") {
        eos.EntropyFromDensityInternalEnergy(v_density, v_energy, v_entropy, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_entropy, v_entropy);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned S(rho, e) should be equal to the true value") {
          array_compare(num, density, energy, h_entropy, entropy_true, "Density",
                        "Energy");
        }
      }
      WHEN("A S(rho, T(rho,e)) lookup is performed") {
        eos.TemperatureFromDensityInternalEnergy(v_density, v_energy, v_local_temp, num);
        eos.EntropyFromDensityTemperature(v_density, v_local_temp, v_entropy, num);
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_entropy, v_entropy);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned S(rho, e) should be equal to the true value") {
          array_compare(num, density, energy, h_entropy, entropy_true, "Density",
                        "Energy");
        }
      }
    }
  }
}

SCENARIO("CarnahanStarling6", "[CarnahanStarling][CarnahanStarling6]") {
  GIVEN("Parameters for a CarnahanStarling EOS") {
    constexpr Real gm1 = 0.4;
    constexpr Real Cv = 7180000.0;
    constexpr Real bb = 1.e-3;
    constexpr Real qq = 0.0;
    // Create the EOS
    EOS host_eos = CarnahanStarling(gm1, Cv, bb, qq);
    EOS eos = host_eos.GetOnDevice();
    GIVEN("Pressure and temperature") {
      constexpr int num = 4;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create Kokkos views on device for the input arrays
      Kokkos::View<Real[num]> v_pressure("Pressure");
      Kokkos::View<Real[num]> v_temperature("Temperature");
#else
      // Otherwise just create arrays to contain values and create pointers to
      // be passed to the functions in place of the Kokkos views
      std::array<Real, num> pressure;
      std::array<Real, num> temperature;
      auto v_pressure = pressure.data();
      auto v_temperature = temperature.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

      // Populate the input arrays
      portableFor(
          "Initialize pressure and temperature", 0, 1, PORTABLE_LAMBDA(int i) {
            v_pressure[0] = 1.0132499999999999e+06;
            v_pressure[1] = 3.4450500000000000e+07;
            v_pressure[2] = 6.7887750000000000e+07;
            v_pressure[3] = 1.0132500000000000e+08;
            v_temperature[0] = 2.9814999999999998e+02;
            v_temperature[1] = 1.5320999999999999e+03;
            v_temperature[2] = 2.7660500000000002e+03;
            v_temperature[3] = 4.0000000000000000e+03;
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create host-side mirrors of the inputs and copy the inputs. These are
      // just used for the comparisons
      auto pressure = Kokkos::create_mirror_view(v_pressure);
      auto temperature = Kokkos::create_mirror_view(v_temperature);
      Kokkos::deep_copy(pressure, v_pressure);
      Kokkos::deep_copy(temperature, v_temperature);
#endif // PORTABILITY_STRATEGY_KOKKOS

      // Gold standard values for a subset of lookups
      constexpr std::array<Real, num> density_true{
          1.1833012071291069e-03, 7.8290736890381501e-03, 8.5453943327882340e-03,
          8.8197619601121831e-03};

#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create device views for outputs and mirror those views on the host
      Kokkos::View<Real[num]> v_density("Density");
      Kokkos::View<Real[num]> v_energy("Energy");
      auto h_density = Kokkos::create_mirror_view(v_density);
      auto h_energy = Kokkos::create_mirror_view(v_energy);
#else
      // Create arrays for the outputs and then pointers to those arrays that
      // will be passed to the functions in place of the Kokkos views
      std::array<Real, num> h_density;
      std::array<Real, num> h_energy;
      // Just alias the existing pointers
      auto v_density = h_density.data();
      auto v_energy = h_energy.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

      WHEN("A rho(P, T) lookup is performed") {
        portableFor("rho(P, T) FillEos lookup", 0, num, PORTABLE_LAMBDA(int i) {
          Real cv, bmod;
          static constexpr const unsigned long _output =
              singularity::thermalqs::density |
              singularity::thermalqs::specific_internal_energy;
          eos.FillEos(v_density[i], v_temperature[i], v_energy[i], v_pressure[i], cv,
                      bmod, _output);
        });

#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_density, v_density);
#endif // PORTABILITY_STRATEGY_KOKKOS

        THEN("The returned rho(P, T) should be equal to the true value") {
          array_compare(num, pressure, temperature, h_density, density_true, "Pressure",
                        "Temperature");
        }
      }
    }
  }
}
