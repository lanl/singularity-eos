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

/*
  Profiling tool for EOSPAC EOS
  Authors: Jonah Miller, Chad Meyers, Josh Dolence, Daniel Holiday,
           Richard Berger
  Copyright: LANL
 */

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>

#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_builder.hpp>

using namespace singularity;

using duration = std::chrono::microseconds;
using dvec = std::vector<double>;
using ivec = std::vector<int>;
using Spiner::RegularGrid1D;
using DataBox = Spiner::DataBox<Real>;

#ifdef PORTABILITY_STRATEGY_KOKKOS
using RView = Kokkos::View<Real *>;
using RMirror = typename RView::HostMirror;
#endif

constexpr Real RHO_MIN = 1e-2; // g/cm^3
constexpr Real RHO_MAX = 1e0;
constexpr Real T_MIN = 1.8e2; // Kelvin
constexpr Real T_MAX = 1e3;
constexpr Real E_MIN = 1e9;
constexpr Real E_MAX = 1e13;

inline auto now() { return std::chrono::high_resolution_clock::now(); }

template <typename Function>
inline double get_duration(Function function) {
  auto start = now();
  function();
  auto stop = now();
  return std::chrono::duration_cast<duration>(stop - start).count();
}

// This file can be used to test any of the EOS methods below by defining
// METHOD_UNDER_TEST. CMake will configure multiple targets each defining one
// of these.

// #define METHOD_UNDER_TEST TemperatureFromDensityInternalEnergy
// #define METHOD_UNDER_TEST PressureFromDensityTemperature
// #define METHOD_UNDER_TEST InternalEnergyFromDensityTemperature
// #define METHOD_UNDER_TEST PressureFromDensityInternalEnergy
// #define METHOD_UNDER_TEST SpecificHeatFromDensityTemperature
// #define METHOD_UNDER_TEST SpecificHeatFromDensityInternalEnergy
// #define METHOD_UNDER_TEST BulkModulusFromDensityTemperature
// #define METHOD_UNDER_TEST BulkModulusFromDensityInternalEnergy
// #define METHOD_UNDER_TEST GruneisenParamFromDensityTemperature
// #define METHOD_UNDER_TEST GruneisenParamFromDensityInternalEnergy

#if METHOD_UNDER_TEST == PressureFromDensityTemperature ||                               \
    METHOD_UNDER_TEST == InternalEnergyFromDensityTemperature ||                         \
    METHOD_UNDER_TEST == SpecificHeatFromDensityTemperature ||                           \
    METHOD_UNDER_TEST == BulkModulusFromDensityTemperature ||                            \
    METHOD_UNDER_TEST == GruneisenParamFromDensityTemperature
constexpr Real X_MIN = RHO_MIN;
constexpr Real X_MAX = RHO_MAX;
constexpr Real Y_MIN = T_MIN;
constexpr Real Y_MAX = T_MAX;
#elif METHOD_UNDER_TEST == TemperatureFromDensityInternalEnergy ||                       \
    METHOD_UNDER_TEST == PressureFromDensityInternalEnergy ||                            \
    METHOD_UNDER_TEST == SpecificHeatFromDensityInternalEnergy ||                        \
    METHOD_UNDER_TEST == BulkModulusFromDensityInternalEnergy ||                         \
    METHOD_UNDER_TEST == GruneisenParamFromDensityInternalEnergy
constexpr Real X_MIN = RHO_MIN;
constexpr Real X_MAX = RHO_MAX;
constexpr Real Y_MIN = E_MIN;
constexpr Real Y_MAX = E_MAX;
#endif

#define xstr(s) str(s)
#define str(s) #s

template <typename T>
inline bool get_timing(int ncycles, const ivec &ncells_1d, const double x_min,
                       const double x_max, const double y_min, const double y_max,
                       const std::string &name, T &eos) {
  // use variant for testing
  EOS eos_h = eos;
  bool success = true;
  std::cout << "\t...Profiling EOS: " << name << " (" << xstr(METHOD_UNDER_TEST)
            << ") ..." << std::endl;

  std::ofstream prof_file;
  std::random_device rd; // seed
  std::mt19937 e2(rd()); // generator
  std::uniform_real_distribution<> dist_x(x_min, x_max);
  std::uniform_real_distribution<> dist_y(y_min, y_max);

  ivec ncells(ncells_1d.size());
  dvec zonecycles(ncells_1d.size());
  dvec durations_host(ncells_1d.size());
  dvec durations_host_vec(ncells_1d.size());
  dvec durations_host_vec_scratch(ncells_1d.size());
  dvec durations_device(ncells_1d.size());
  dvec throughputs_host(ncells_1d.size());
  dvec throughputs_host_vec(ncells_1d.size());
  dvec throughputs_host_vec_scratch(ncells_1d.size());
  dvec throughputs_device(ncells_1d.size());

  for (int i = 0; i < ncells_1d.size(); i++) {
    ncells[i] = ncells_1d[i] * ncells_1d[i] * ncells_1d[i];
    zonecycles[i] = static_cast<double>(ncycles * ncells[i]);
  }

  EOS eos_d = eos_h.GetOnDevice();

  for (int n = 0; n < ncells_1d.size(); n++) {
    std::cout << "\t\t...ncells 1d = " << ncells_1d[n] << std::endl;

    const int nc1d = ncells_1d[n];
    const int nc = ncells[n];
    dvec scratch(T::scratch_size(xstr(METHOD_UNDER_TEST), nc) / sizeof(double));
    DataBox F_h(nc1d, nc1d, nc1d);
    DataBox F_hvec(nc1d, nc1d, nc1d);
    DataBox F_hscratch(nc1d, nc1d, nc1d);
#ifdef PORTABILITY_STRATEGY_KOKKOS
    RView x_v("x", nc); // first views
    RView y_v("y", nc);
    RView F_v("F", nc);
    RMirror x_v_hm = Kokkos::create_mirror_view(x_v);
    RMirror y_v_hm = Kokkos::create_mirror_view(y_v);
    RMirror F_v_hm = Kokkos::create_mirror_view(F_v);
    DataBox x_d(x_v.data(), nc1d, nc1d, nc1d); // wrap viws in databoxes
    DataBox y_d(y_v.data(), nc1d, nc1d, nc1d);
    DataBox F_d(F_v.data(), nc1d, nc1d, nc1d);
    DataBox x_h(x_v_hm.data(), nc1d, nc1d, nc1d);
    DataBox y_h(y_v_hm.data(), nc1d, nc1d, nc1d);
    DataBox F_hm(F_v_hm.data(), nc1d, nc1d, nc1d);
#else
    DataBox x_h(nc1d, nc1d, nc1d);
    DataBox y_h(nc1d, nc1d, nc1d);
    DataBox x_d = x_h.slice(3, 0, nc1d);
    DataBox y_d = y_h.slice(3, 0, nc1d);
    DataBox F_d(nc1d, nc1d, nc1d);
    DataBox F_hm = F_d.slice(3, 0, nc1d);
#endif

    for (int k = 0; k < nc1d; k++) {
      for (int j = 0; j < nc1d; j++) {
        for (int i = 0; i < nc1d; i++) {
          x_h(k, j, i) = dist_x(e2);
          y_h(k, j, i) = dist_y(e2);
        }
      }
    }
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::deep_copy(x_v, x_v_hm);
    Kokkos::deep_copy(y_v, y_v_hm);
#endif

    durations_host[n] = get_duration([&]() {
      for (int cycle = 0; cycle < ncycles; cycle++) {
        for (int k = 0; k < nc1d; k++) {
          for (int j = 0; j < nc1d; j++) {
            for (int i = 0; i < nc1d; i++) {
              F_h(k, j, i) = eos_h.METHOD_UNDER_TEST(x_h(k, j, i), y_h(k, j, i));
            }
          }
        }
      }
    });
    durations_host[n] /= static_cast<Real>(nc * ncycles);

    durations_host_vec[n] = get_duration([&]() {
      for (int cycle = 0; cycle < ncycles; cycle++) {
        eos_h.METHOD_UNDER_TEST(x_h.data(), y_h.data(), F_hvec.data(), nc);
      }
    });
    durations_host_vec[n] /= static_cast<Real>(nc * ncycles);

    durations_host_vec_scratch[n] = get_duration([&]() {
      for (int cycle = 0; cycle < ncycles; cycle++) {
        eos_h.METHOD_UNDER_TEST(x_h.data(), y_h.data(), F_hscratch.data(), scratch.data(),
                                nc);
      }
    });
    durations_host_vec_scratch[n] /= static_cast<Real>(nc * ncycles);

    //    durations_device[n] = get_duration([&]() {
    //      portableFor(
    //          "pressure from density and temperature", 0, ncycles, 0, nc1d, 0, nc1d, 0,
    //          nc1d, PORTABLE_LAMBDA(const int &n, const int &k, const int &j, const int
    //          &i) {
    //            F_d(k, j, i) =
    //                eos_d.METHOD_UNDER_TEST(x_d(k, j, i), y_d(k, j, i));
    //          });
    //    });
    //    durations_device[n] /= static_cast<Real>(nc * ncycles);
    //
    throughputs_host[n] = 1e6 / durations_host[n];
    throughputs_host_vec[n] = 1e6 / durations_host_vec[n];
    throughputs_host_vec_scratch[n] = 1e6 / durations_host_vec_scratch[n];
    //    throughputs_device[n] = 1e6 / durations_device[n];
    //
    // #ifdef PORTABILITY_STRATEGY_KOKKOS
    //    Kokkos::deep_copy(F_v_hm, F_v);
    // #endif
    Real F_min = std::numeric_limits<Real>::max();
    Real F_max = std::numeric_limits<Real>::min();
    for (int k = 0; k < nc1d; k++) {
      for (int j = 0; j < nc1d; j++) {
        for (int i = 0; i < nc1d; i++) {
          const auto val = F_h(k, j, i);
          if (val < F_min) F_min = val;
          if (val > F_max) F_max = val;
        }
      }
    }
    std::cout << "\t\t\t...value range = [" << F_min << ", " << F_max << "]" << std::endl;

    F_min = std::numeric_limits<Real>::max();
    F_max = std::numeric_limits<Real>::min();
    for (int k = 0; k < nc1d; k++) {
      for (int j = 0; j < nc1d; j++) {
        for (int i = 0; i < nc1d; i++) {
          const auto val = F_hscratch(k, j, i);
          if (val < F_min) F_min = val;
          if (val > F_max) F_max = val;
        }
      }
    }
    std::cout << "\t\t\t...value range (scratch) = [" << F_min << ", " << F_max << "]"
              << std::endl;

    F_min = std::numeric_limits<Real>::max();
    F_max = std::numeric_limits<Real>::min();
    for (int k = 0; k < nc1d; k++) {
      for (int j = 0; j < nc1d; j++) {
        for (int i = 0; i < nc1d; i++) {
          const auto val = F_hvec(k, j, i);
          if (val < F_min) F_min = val;
          if (val > F_max) F_max = val;
        }
      }
    }
    std::cout << "\t\t\t...value range (vector w scalar calls) = [" << F_min << ", "
              << F_max << "]" << std::endl;

    Real max_diff = 0;
    for (int k = 0; k < nc1d; k++) {
      for (int j = 0; j < nc1d; j++) {
        for (int i = 0; i < nc1d; i++) {
          Real diff = (100. * std::abs(F_h(k, j, i) - F_hscratch(k, j, i)) /
                       (std::abs(F_h(k, j, i)) + 1e-20));
          if (diff > max_diff) max_diff = diff;
        }
      }
    }
    std::cout << "\t\t\t...max difference (scratch) = " << max_diff << "%" << std::endl;
    success = success && (max_diff < 1.0);

    max_diff = 0;
    for (int k = 0; k < nc1d; k++) {
      for (int j = 0; j < nc1d; j++) {
        for (int i = 0; i < nc1d; i++) {
          Real diff = (100. * std::abs(F_h(k, j, i) - F_hvec(k, j, i)) /
                       (std::abs(F_h(k, j, i)) + 1e-20));
          if (diff > max_diff) max_diff = diff;
        }
      }
    }
    std::cout << "\t\t\t...max difference (vector w scalar calls) = " << max_diff << "%"
              << std::endl;

    if (durations_host_vec_scratch[n] >= durations_host_vec[n]) {
      std::cout << "\t\t\t...ERROR: vector calls too slow" << std::endl;
    }

    success = success && (max_diff < 1.0) &&
              durations_host_vec_scratch[n] < durations_host_vec[n];
  }

  prof_file.open(name + "_" + xstr(METHOD_UNDER_TEST) + "_timing.dat",
                 std::ios::out | std::ios::trunc);
  prof_file
      << "# EOS name = " << name << "\n"
      << "# EOS method = " << xstr(METHOD_UNDER_TEST) << "\n"
      << "# Columns:\n"
      << "# [0]: ncells/axis\n"
      << "# [1]: time/point on host   (microseconds)\n"
      << "# [2]: throughput on host   (zone-cycles/node-second)\n"
      << "# [3]: vector functions (w scalar calls) time/point on host   (microseconds)\n"
      << "# [4]: vector functions (w scalar calls) throughput on host   "
         "(zone-cycles/node-second)\n"
      << "# [5]: vectorized (SCRATCH) time/point on host   (microseconds)\n"
      << "# [6]: vectorized (SCRATCH) throughput on host   (zone-cycles/node-second)\n"
      << std::endl;
  for (int n = 0; n < ncells_1d.size(); n++) {
    prof_file << std::setw(10) << ncells_1d[n] << "\t" << durations_host[n] << "\t"
              << throughputs_host[n] << "\t" << durations_host_vec[n] << "\t"
              << throughputs_host_vec[n] << "\t" << durations_host_vec_scratch[n] << "\t"
              << throughputs_host_vec_scratch[n]
              << "\t"
              // << durations_device[n] << "\t"
              // << throughputs_device[n]
              << std::endl;
  }
  prof_file.close();

  eos_d.Finalize();

  return success;
}

int main(int argc, char *argv[]) {
  bool success = true;
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize(argc, argv);
#endif
  {
    ivec ncells_1d;
    if (argc < 3) {
      std::cerr << argv[0] << ": ncycles ncells1, ncells2, ncells3, ..." << std::endl;
      return 1;
    }

    int ncycles = std::atoi(argv[1]);

    ncells_1d.resize(argc - 2);
    for (int i = 2; i < argc; i++) {
      ncells_1d[i - 2] = std::atoi(argv[i]);
    }

    std::cout << "Welcome to the profiler.\n"
              << "Please note that timings are for 1 node and 1 GPU respectively.\n"
              << "Beginning profiling..." << std::endl;

    constexpr int airID = 5030;

    {
      auto eos = EOSPAC(airID);
      success = get_timing(ncycles, ncells_1d, X_MIN, X_MAX, Y_MIN, Y_MAX, "EOSPAC", eos);
    }

    {
      auto eos = ShiftedEOS<EOSPAC>(EOSPAC(airID), 1.0);
      success = get_timing(ncycles, ncells_1d, X_MIN, X_MAX, Y_MIN, Y_MAX,
                           "ShiftedEOS<EOSPAC>", eos) &&
                success;
    }

    {
      auto eos = ScaledEOS<EOSPAC>(EOSPAC(airID), 2.0);
      success = get_timing(ncycles, ncells_1d, 0.5 * X_MIN, 0.5 * X_MAX, Y_MIN, Y_MAX,
                           "ScaledEOS<EOSPAC>", eos) &&
                success;
    }

    {
      auto eos =
          ScaledEOS<ShiftedEOS<EOSPAC>>(ShiftedEOS<EOSPAC>(EOSPAC(airID), 0.3), 2.0);
      success = get_timing(ncycles, ncells_1d, 0.5 * X_MIN, 0.5 * X_MAX, Y_MIN, Y_MAX,
                           "ScaledEOS<ShiftedEOS<EOSPAC>>", eos) &&
                success;
    }

    {
      auto eos = BilinearRampEOS<EOSPAC>(EOSPAC(airID), 1.0, 1.0, 0.0, 0.0);
      success = get_timing(ncycles, ncells_1d, X_MIN, X_MAX, Y_MIN, Y_MAX,
                           "BilinearRampEOS<EOSPAC>", eos) &&
                success;
    }

    {
      auto eos = BilinearRampEOS<ScaledEOS<EOSPAC>>(ScaledEOS<EOSPAC>(EOSPAC(airID), 2.0),
                                                    1.0, 1.0, 0.0, 0.0);
      success = get_timing(ncycles, ncells_1d, 0.5 * X_MIN, 0.5 * X_MAX, Y_MIN, Y_MAX,
                           "BilinearRampEOS<ScaledEOS<EOSPAC>>", eos) &&
                success;
    }

    std::cout << "Done." << std::endl;
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif

  return success ? 0 : 1;
}
