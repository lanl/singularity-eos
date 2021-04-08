//------------------------------------------------------------------------------
// © 2021. Triad National Security, LLC. All rights reserved.  This
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
  Profiling tool for EOS
  Authors: Jonah Miller, Chad Meyers, Josh Dolence, Daniel Holiday
  Copyright: LANL
 */

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

#include "../utils/ports-of-call/portability.hpp"
#include "../utils/ports-of-call/portable_arrays.hpp"

#include "../utils/spiner/databox.hpp"
#include "../utils/spiner/interpolation.hpp"

#include "../eos/eos.hpp"
#include "../eos/eos_builder.hpp"

using namespace singularity;

using duration = std::chrono::microseconds;
using dvec = std::vector<double>;
using ivec = std::vector<int>;
using Spiner::DataBox;
using Spiner::RegularGrid1D;

#ifdef PORTABILITY_STRATEGY_KOKKOS
using RView = Kokkos::View<Real*>;
using RMirror = typename RView::HostMirror;
#endif

constexpr Real RHO_MIN = 1e-2; // g/cm^3
constexpr Real RHO_MAX = 1e2;
constexpr Real T_MIN   = 1e2;  // Kelvin
constexpr Real T_MAX   = 1e4;

inline auto now() {
  return std::chrono::high_resolution_clock::now();
}

template<typename Function>
inline double get_duration(Function function) {
  auto start = now();
  function();
  auto stop = now();
  return std::chrono::duration_cast<duration>(stop-start).count();
}

/*
  TODO(JMM): Only profiles Pressure call
  and does not accept lambdas.
*/
inline void get_timing(int ncycles,
                       const ivec& ncells_1d,
                       const double rho_min, const double rho_max,
                       const double T_min, const double T_max,
                       const std::string& name,
                       EOS& eos_h) {

  std::cout << "\t...Profiling EOS: " << name << "..." << std::endl;

  std::ofstream prof_file;
  std::random_device rd; // seed
  std::mt19937 e2(rd()); // generator
  std::uniform_real_distribution<> dist_rho(rho_min,rho_max);
  std::uniform_real_distribution<> dist_T(T_min,T_max);

  ivec ncells(ncells_1d.size());
  dvec zonecycles(ncells_1d.size());
  dvec durations_host(ncells_1d.size());
  dvec durations_device(ncells_1d.size());
  dvec throughputs_host(ncells_1d.size());
  dvec throughputs_device(ncells_1d.size());

  for (int i = 0; i < ncells_1d.size(); i++) {
    ncells[i] = ncells_1d[i]*ncells_1d[i]*ncells_1d[i];
    zonecycles[i] = static_cast<double>(ncycles*ncells[i]);
  }
  

  EOS eos_d = eos_h.GetOnDevice();

  for (int n = 0; n < ncells_1d.size(); n++) {
    
    std::cout << "\t\t...ncells 1d = " << ncells_1d[n] << std::endl;

    const int nc1d = ncells_1d[n];
    const int nc = ncells[n];
    DataBox P_h(nc1d,nc1d,nc1d);
    #ifdef PORTABILITY_STRATEGY_KOKKOS
    RView rho_v("rho",nc); // first views
    RView T_v("T",nc);
    RView P_v("P",nc);
    RMirror rho_v_hm = Kokkos::create_mirror_view(rho_v);
    RMirror T_v_hm = Kokkos::create_mirror_view(T_v);
    RMirror P_v_hm = Kokkos::create_mirror_view(P_v);
    DataBox rho_d(rho_v.data(),nc1d,nc1d,nc1d); // wrap viws in databoxes
    DataBox T_d(T_v.data(),nc1d,nc1d,nc1d);
    DataBox P_d(P_v.data(),nc1d,nc1d,nc1d);
    DataBox rho_h(rho_v_hm.data(),nc1d,nc1d,nc1d);
    DataBox T_h(T_v_hm.data(),nc1d,nc1d,nc1d);
    DataBox P_hm(P_v_hm.data(),nc1d,nc1d,nc1d);
    #else
    DataBox rho_h(nc1d,nc1d,nc1d);
    DataBox T_h(nc1d,nc1d,nc1d);
    DataBox rho_d = rho_h.slice(3,0,nc1d);
    DataBox T_d = T_h.slice(3,0,nc1d);
    DataBox P_d(nc1d,nc1d,nc1d);
    DataBox P_hm = P_d.slice(3,0,nc1d);
    #endif

    for (int k = 0; k < nc1d; k++) {
      for (int j = 0; j < nc1d; j++) {
        for (int i = 0; i < nc1d; i++) {
          rho_h(k,j,i) = dist_rho(e2);
          T_h(k,j,i) = dist_T(e2);
        }
      }
    }
    #ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::deep_copy(rho_v,rho_v_hm);
    Kokkos::deep_copy(T_v,T_v_hm);
    #endif

    durations_host[n] = get_duration([&]() {
        for (int cycle = 0; cycle < ncycles; cycle++) {
          #pragma omp parallel
          for (int k = 0; k < nc1d; k++) {
            for (int j = 0; j< nc1d; j++) {
              for (int i = 0; i < nc1d; i++) {
                P_h(k,j,i) = eos_h.PressureFromDensityTemperature(rho_h(k,j,i),
                                                                  T_h(k,j,i));
              }
            }
          }
        }
      });
    durations_host[n] /= static_cast<Real>(nc*ncycles);

    durations_device[n] = get_duration([&]() {
        portableFor("pressure from density and temperature",
                    0,ncycles,0,nc1d,0,nc1d,0,nc1d,
                    PORTABLE_LAMBDA(const int& n, const int& k,
                                    const int& j, const int& i) {
                      P_d(k,j,i) = eos_d.PressureFromDensityTemperature(rho_d(k,j,i),
                                                                        T_d(k,j,i));
                    });
      });
    durations_device[n] /= static_cast<Real>(nc*ncycles);

    throughputs_host[n] = 1e6/durations_host[n];
    throughputs_device[n] = 1e6/durations_device[n];

    #ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::deep_copy(P_v_hm, P_v);
    #endif
    Real max_diff = 0;
    for (int k = 0; k < nc1d; k++) {
      for (int j = 0; j < nc1d; j++) {
        for (int i = 0; i < nc1d; i++) {
          Real diff = (100.*std::abs(P_h(k,j,i) - P_hm(k,j,i))
                       /(std::abs(P_h(k,j,i)) + 1e-20));
          if (diff > max_diff) max_diff = diff;
        }
      }
    }
    std::cout << "\t\t\t...max difference = "
              << max_diff << "%"
              << std::endl;
  }

  prof_file.open(name + "_timing.dat",
                 std::ios::out | std::ios::trunc);
  prof_file << "# EOS name = " << name << "\n"
            << "# Columns:\n"
            << "# [0]: ncells/axis\n"
            << "# [1]: time/point on host   (microseconds)\n"
            << "# [2]: throughput on host   (zone-cycles/node-second)\n"
            << "# [3]: time/point on device (microseconds)\n"
            << "# [4]: throughput on device (zone-cycles/GPU-second)"
            << std::endl;
  for (int n = 0; n < ncells_1d.size(); n++) {
    prof_file << std::setw(10)
              << ncells_1d[n] << "\t"
              << durations_host[n] << "\t"
              << throughputs_host[n] << "\t"
              << durations_device[n] << "\t"
              << throughputs_device[n]
              << std::endl;
  }
  prof_file.close();

  eos_d.Finalize();
}

int main(int argc, char* argv[]) {

  #ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize(argc,argv);
  #endif
  {
    ivec ncells_1d;
    if (argc < 3) {
      std::cerr << argv[0]
                << ": ncycles ncells1, ncells2, ncells3, ..."
                << std::endl;
      return 1;
    }
    
    int ncycles = std::atoi(argv[1]);
    
    ncells_1d.resize(argc-2);
    for (int i = 2; i < argc; i++) {
      ncells_1d[i-2] = std::atoi(argv[i]);
    }
    
    std::cout << "Welcome to the profiler.\n"
              << "Please note that timings are for 1 node and 1 GPU respectively.\n"
              << "Beginning profiling..." << std::endl;
    
    EOSBuilder::EOSType type = EOSBuilder::EOSType::IdealGas;
    EOSBuilder::params_t params;
    params["Cv"] = 1e7*0.716; // specific heat in ergs/(g K)
    params["gm1"] = 0.4; // gamma - 1
    EOS eos = EOSBuilder::buildEOS(type,params);
    get_timing(ncycles, ncells_1d,
               RHO_MIN, RHO_MAX, T_MIN, T_MAX,
               "IdealGas", eos);
    
    std::cout << "Done." << std::endl;
  }
  #ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
  #endif

  return 0;
}
