//------------------------------------------------------------------------------
// © 2021-2025. Triad National Security, LLC. All rights reserved.  This
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

#include <chrono>
#include <iostream>
#include <memory>
#include <stdlib.h>
#include <time.h>
#include <utility>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <pte_test_3mat_analytic.hpp>
#include <pte_test_utils.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>
#include <spiner/databox.hpp>

#include <singularity-eos/eos/eos_models.hpp>
#include <singularity-eos/eos/eos_variant.hpp>

#include "eos_unit_test_helpers.hpp"

using DataBox = Spiner::DataBox<Real>;
using singularity::PTESolverPT;
using singularity::PTESolverPTRequiredScratch;
using singularity::PTESolverRhoT;
using singularity::PTESolverRhoTRequiredScratch;
using singularity::Variant;
using EOS = Variant<Gruneisen, DavisReactants, DavisProducts>;

#ifdef PORTABILITY_STRATEGY_KOKKOS
using RView = Kokkos::View<Real *>;
using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
#endif

template <template <typename... Types> class Method_t>
auto TestPTE(const std::string name, const std::size_t nscratch_vars) {
  constexpr Real EPS = 1e-5;
  Real time;
  std::size_t nsuccess = 0;

  // EOS
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::View<EOS *> eos_v("eos", NMAT);
  auto eos_hv = Kokkos::create_mirror_view(eos_v);
#else
  std::vector<EOS> eos_vec(NMAT);
  PortableMDArray<EOS> eos_hv(eos_vec.data(), NMAT);
  PortableMDArray<EOS> eos_v(eos_vec.data(), NMAT);
#endif

  LinearIndexer<decltype(eos_hv)> eos_h(eos_hv);
  set_eos(eos_hv.data());

#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::deep_copy(eos_v, eos_hv);
#endif

  using EOSAccessor = LinearIndexer<decltype(eos_v)>;
  EOSAccessor eos(eos_v);

#ifdef PORTABILITY_STRATEGY_KOKKOS
  RView rho_v("rho", NPTS);
  RView vfrac_v("vfrac", NPTS);
  RView sie_v("sie", NPTS);
  RView temp_v("temp", NPTS);
  RView press_v("press", NPTS);
  RView scratch_v("scratch", NTRIAL * nscratch_vars);
  auto rho_vh = Kokkos::create_mirror_view(rho_v);
  auto vfrac_vh = Kokkos::create_mirror_view(vfrac_v);
  auto sie_vh = Kokkos::create_mirror_view(sie_v);
  auto temp_vh = Kokkos::create_mirror_view(temp_v);
  auto press_vh = Kokkos::create_mirror_view(press_v);
  auto scratch_vh = Kokkos::create_mirror_view(scratch_v);
  DataBox rho_d(rho_v.data(), NTRIAL, NMAT);
  DataBox vfrac_d(vfrac_v.data(), NTRIAL, NMAT);
  DataBox sie_d(sie_v.data(), NTRIAL, NMAT);
  DataBox temp_d(temp_v.data(), NTRIAL, NMAT);
  DataBox press_d(press_v.data(), NTRIAL, NMAT);
  DataBox scratch_d(scratch_v.data(), NTRIAL * nscratch_vars);
  DataBox rho_hm(rho_vh.data(), NTRIAL, NMAT);
  DataBox vfrac_hm(vfrac_vh.data(), NTRIAL, NMAT);
  DataBox sie_hm(sie_vh.data(), NTRIAL, NMAT);
  DataBox temp_hm(temp_vh.data(), NTRIAL, NMAT);
  DataBox press_hm(press_vh.data(), NTRIAL, NMAT);
  DataBox scratch_hm(scratch_vh.data(), NTRIAL * nscratch_vars);

  // histogram
  Kokkos::View<int *, atomic_view> hist_d("histogram", HIST_SIZE);
  auto hist_vh = Kokkos::create_mirror_view(hist_d);
#else
  DataBox rho_d(NTRIAL, NMAT);
  DataBox vfrac_d(NTRIAL, NMAT);
  DataBox sie_d(NTRIAL, NMAT);
  DataBox temp_d(NTRIAL, NMAT);
  DataBox press_d(NTRIAL, NMAT);
  DataBox scratch_d(NTRIAL, nscratch_vars);
  DataBox rho_hm = rho_d.slice(2, 0, NTRIAL);
  DataBox vfrac_hm = vfrac_d.slice(2, 0, NTRIAL);
  DataBox sie_hm = sie_d.slice(2, 0, NTRIAL);
  DataBox temp_hm = temp_d.slice(2, 0, NTRIAL);
  DataBox press_hm = press_d.slice(2, 0, NTRIAL);
  DataBox scratch_hm = scratch_d.slice(2, 0, NTRIAL);

  // histogram
  int hist_vh[HIST_SIZE];
  int *hist_d = hist_vh;
#endif

  // setup state
  for (int n = 0; n < NTRIAL; n++) {
    Indexer2D<decltype(rho_hm)> r(n, rho_hm);
    Indexer2D<decltype(vfrac_hm)> vf(n, vfrac_hm);
    Indexer2D<decltype(sie_hm)> e(n, sie_hm);
    Indexer2D<decltype(temp_hm)> t(n, temp_hm);
    set_state(r, vf, e, t, eos_h);
  }
  for (int i = 0; i < HIST_SIZE; ++i) {
    hist_vh[i] = 0;
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::deep_copy(rho_v, rho_vh);
  Kokkos::deep_copy(vfrac_v, vfrac_vh);
  Kokkos::deep_copy(sie_v, sie_vh);
  Kokkos::deep_copy(temp_v, temp_vh);
  Kokkos::deep_copy(press_v, press_vh);
  Kokkos::deep_copy(hist_d, hist_vh);
#endif

  auto start = std::chrono::high_resolution_clock::now();
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::fence();
#endif
  std::cout << "Starting solver " << name << " with " << NTRIAL << " trials."
            << std::endl;

  portableReduce(
      "PTE!", 0, NTRIAL,
      PORTABLE_LAMBDA(const int &t, std::size_t &ns) {
        Real *lambda[NMAT];
        for (int i = 0; i < NMAT; i++) {
          lambda[i] = nullptr;
        }

        Indexer2D<decltype(rho_d)> rho(t, rho_d);
        Indexer2D<decltype(vfrac_d)> vfrac(t, vfrac_d);
        Indexer2D<decltype(sie_d)> sie(t, sie_d);
        Indexer2D<decltype(temp_d)> temp(t, temp_d);
        Indexer2D<decltype(press_d)> press(t, press_d);

        Real sie_tot = 0.0;
        Real rho_tot = 0.0;
        for (int i = 0; i < NMAT; i++) {
          rho_tot += rho[i] * vfrac[i];
          sie_tot += rho[i] * vfrac[i] * sie[i];
        }
        sie_tot /= rho_tot;

        const Real Tguess =
            ApproxTemperatureFromRhoMatU(NMAT, eos, rho_tot * sie_tot, rho, vfrac);

        auto method = Method_t<EOSAccessor, Indexer2D<decltype(rho_d)>, decltype(lambda)>(
            NMAT, eos, 1.0, sie_tot, rho, vfrac, sie, temp, press, lambda,
            &scratch_d(t * nscratch_vars), Tguess);
        auto status = PTESolver(method);
        if (status.converged) {
          bool in_pte = true;
          for (std::size_t m = 0; m < NMAT - 1; ++m) {
            bool press_close = isClose(press[m], press[m + 1], EPS);
            if (!press_close) {
              printf("Pressures not close! %ld %.14e %.14e %.14e\n", m, press[m],
                     press[m + 1], press[m] - press[m + 1]);
            }
            in_pte = in_pte && press_close;
            bool temp_close = isClose(temp[m], temp[m + 1], EPS);
            if (!temp_close) {
              printf("Temperatures not close! %ld, %.14e %.14e %.14e\n", m, temp[m],
                     temp[m + 1], temp[m] - temp[m + 1]);
            }
            in_pte = in_pte && temp_close;
          }
          ns += in_pte;
        }
        hist_d[method.Niter()] += 1;
      },
      nsuccess);
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::fence();
#endif
  auto stop = std::chrono::high_resolution_clock::now();
  auto sum_time = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::deep_copy(hist_vh, hist_d);
#endif

  Real milliseconds = sum_time.count() / 1e3;
  time = milliseconds;

  std::cout << "Finished " << NTRIAL << " solves in " << milliseconds << " milliseconds"
            << std::endl;
  std::cout << "Solves/ms = " << NTRIAL / milliseconds << std::endl;
  std::cout << "Success: " << nsuccess << "   Failure: " << NTRIAL - nsuccess
            << std::endl;
  std::cout << "Histogram:\n"
            << "iters\tcount\n"
            << "----------------------\n";
  for (int i = 0; i < HIST_SIZE; ++i) {
    std::cout << i << "\t" << hist_vh[i] << "\n";
  }
  std::cout << std::endl;

  return std::make_pair(nsuccess, rho_d);
}

int main(int argc, char *argv[]) {
  int nsuccess = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
#endif
  {
    srand(time(NULL));

    // scratch required for PTE solver
    auto nscratch_vars_rt = PTESolverRhoTRequiredScratch(NMAT);
    auto [ns_rt, rho_rt] = TestPTE<PTESolverRhoT>("PTESolverRhoT", nscratch_vars_rt);

    // // scratch required for PTE solver
    auto nscratch_vars_pt = PTESolverPTRequiredScratch(NMAT);
    auto [ns_pt, rho_pt] = TestPTE<PTESolverPT>("PTESolverPT", nscratch_vars_pt);
    nsuccess += ns_rt;
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif

  // poor-man's ctest integration
  return (nsuccess >= 0.5 * NTRIAL) ? 0 : 1;
}
