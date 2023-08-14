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
#ifndef CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#endif

#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <test/eos_unit_test_helpers.hpp>

using singularity::EOS;
using singularity::NobleAbel;

SCENARIO("NobleAbel1", "[NobleAbel][NobleAbel1]") {
  GIVEN("Parameters for a NobleAbel EOS") {
    constexpr Real gm1 = 0.2110231425091352;
    constexpr Real Cv = 16420000.0;
    constexpr Real bb = 1.435;
    constexpr Real qq = 0.0;
    // Create the EOS
    EOS host_eos = NobleAbel(gm1, Cv, bb, qq);
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
            v_density[0] = 9.7941724217367608e-04;
            v_density[1] = 6.4295356839573397e-03;
            v_density[2] = 7.0119063617845390e-03;
            v_density[3] = 7.2347087589269467e-03;
            v_energy[0] = 4.8956230000000000e+09;
            v_energy[1] = 2.5157082000000000e+10;
            v_energy[2] = 4.5418541000000000e+10;
            v_energy[3] = 6.5680000000000000e+10;
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
          1.0132500000019073e+06, 3.4450500000000000e+07, 6.7887750000001907e+07,
          1.0132500000000000e+08};
      constexpr std::array<Real, num> bulkmodulus_true{
          1.2287962276923475e+06, 4.2108865319896758e+07, 8.3049285363644123e+07,
          1.2399420381644218e+08};
      constexpr std::array<Real, num> temperature_true{
          2.9814999999999998e+02, 1.5320999999999999e+03, 2.7660500000000002e+03,
          4.0000000000000000e+03};
      constexpr std::array<Real, num> gruneisen_true{
          2.1132014531143434e-01, 2.1298825386425976e-01, 2.1316805775971534e-01,
          2.1323692714677225e-01};

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

SCENARIO("NobleAbel2", "[NobleAbel][NobleAbel2]") {
  GIVEN("Parameters for a NobleAbel EOS") {
    constexpr Real gm1 = 0.2110231425091352;
    constexpr Real Cv = 1642.00e+04;
    constexpr Real bb = 1.435;
    constexpr Real qq = 42.00e+09;
    constexpr Real qp = -23.0e+7;
    constexpr Real T0 = 200.0;
    constexpr Real P0 = 1000000.0;
    // Create the EOS
    EOS host_eos = NobleAbel(gm1, Cv, bb, qq, qp, T0, P0);
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
            v_density[0] = 9.7941724217367608e-04;
            v_density[1] = 6.4295356839573397e-03;
            v_density[2] = 7.0119063617845390e-03;
            v_density[3] = 7.2347087589269467e-03;
            v_energy[0] = 4.6895623000000000e+10;
            v_energy[1] = 6.7157082000000000e+10;
            v_energy[2] = 8.7418541000000000e+10;
            v_energy[3] = 1.0768000000000000e+11;
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
          1.0132500000019073e+06, 3.4450500000000000e+07, 6.7887750000001907e+07,
          1.0132500000000000e+08};
      constexpr std::array<Real, num> bulkmodulus_true{
          1.2287962276923475e+06, 4.2108865319896758e+07, 8.3049285363644123e+07,
          1.2399420381644218e+08};
      constexpr std::array<Real, num> temperature_true{
          2.9814999999999998e+02, 1.5320999999999999e+03, 2.7660500000000002e+03,
          4.0000000000000000e+03};
      constexpr std::array<Real, num> gruneisen_true{
          2.1132014531143434e-01, 2.1298825386425976e-01, 2.1316805775971534e-01,
          2.1323692714677225e-01};

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

SCENARIO("NobleAbel Isentropic Bulk Modulus Analytic vs. FD", "[NobleAbel][NobleAbel3]") {
  GIVEN("Parameters for a N-A Gas EOS") {
    constexpr Real gm1 = 2.1102314250913520e-01;
    constexpr Real Cv = 1.6420000000000000e+07;
    constexpr Real bb = 1.4350000000000001e+00;
    constexpr Real qq = 0.0;
    //  Create the EOS
    EOS host_eos = NobleAbel(gm1, Cv, bb, qq);
    EOS eos = host_eos.GetOnDevice();
    GIVEN("Density and energy") {
      constexpr Real density = 9.7941724217367608e-04;
      constexpr Real energy = 4.8956230000000000e+09;
      constexpr Real true_sound_speed = 3.5420612112532610e+04;
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

SCENARIO("Recover Ideal Gas from NobleAbel", "[NobleAbel][NobleAbel4]") {
  GIVEN("Parameters for a NobleAbel EOS") {
    constexpr Real gm1 = 2.1102314250913520e-01;
    constexpr Real Cv = 1.6420000000000000e+07;
    constexpr Real bb = 0;
    constexpr Real qq = 0;
    //  Create the EOS
    EOS host_eos = NobleAbel(gm1, Cv, bb, qq);
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
            v_density[0] = 9.8079571498991256e-04;
            v_energy[0] = 4.8956230000000000e+09;
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

SCENARIO("Test NobleAbel Entropy Calls", "[NobleAbel][NobleAbel5]") {
  GIVEN("Parameters for a NobleAbel EOS") {
    constexpr Real gm1 = 2.1102314250913520e-01;
    constexpr Real Cv = 1.6420000000000000e+07;
    constexpr Real bb = 1.4350000000000001e+00;
    constexpr Real qq = 0.0;
    constexpr Real qp = -23.0e+7;
    constexpr Real T0 = 200.0;
    constexpr Real P0 = 1000000.0;
    //  Create the EOS
    EOS host_eos = NobleAbel(gm1, Cv, bb, qq, qp, T0, P0);
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
            v_density[0] = 6.4295356839573397e-03;
            v_energy[0] = 2.5157082000000000e+10;
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
      constexpr std::array<Real, num> pressure_true{3.4450499999999993e+07};
      constexpr std::array<Real, num> temperature_true{1.5320999999999999e+03};
      constexpr std::array<Real, num> entropy_true{-2.0177705273504910e+08};

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
