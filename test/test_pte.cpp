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

#include <chrono>
#include <iostream>
#include <memory>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <pte_test_3mat_analytic.hpp>
#include <pte_test_utils.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>
#include <spiner/databox.hpp>

#include <singularity-eos/eos/eos_models.hpp>
#include <singularity-eos/eos/eos_variant.hpp>

using DataBox = Spiner::DataBox<Real>;
using singularity::PTESolver;
using singularity::PTESolverPT;
using singularity::PTESolverPTRequiredScratch;
using singularity::PTESolverRhoT;
using singularity::PTESolverRhoTRequiredScratch;
using singularity::Variant;
using EOS = Variant<Gruneisen, DavisReactants, DavisProducts>;

int main(int argc, char *argv[]) {

  int nsuccess = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
#endif
  {
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

    // scratch required for PTE solver
    auto nscratch_vars_rt = PTESolverRhoTRequiredScratch(NMAT);
    auto nscratch_vars_pt = PTESolverPTRequiredScratch(NMAT);

    // state vars
#ifdef PORTABILITY_STRATEGY_KOKKOS
    using RView = Kokkos::View<Real *>;
    using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
    RView rho_v("rho", NPTS);
    RView vfrac_v("vfrac", NPTS);
    RView sie_v("sie", NPTS);
    RView temp_v("temp", NPTS);
    RView press_v("press", NPTS);
    RView scratch_v_rt("scratch", NTRIAL * nscratch_vars_rt);
    RView scratch_v_pt("scratch", NTRIAL * nscratch_vars_pt);
    Kokkos::View<int *, atomic_view> hist_d("histogram", HIST_SIZE);
    auto rho_vh = Kokkos::create_mirror_view(rho_v);
    auto vfrac_vh = Kokkos::create_mirror_view(vfrac_v);
    auto sie_vh = Kokkos::create_mirror_view(sie_v);
    auto temp_vh = Kokkos::create_mirror_view(temp_v);
    auto press_vh = Kokkos::create_mirror_view(press_v);
    auto scratch_vh = Kokkos::create_mirror_view(scratch_v);
    auto hist_vh = Kokkos::create_mirror_view(hist_d);
    DataBox rho_d(rho_v.data(), NTRIAL, NMAT);
    DataBox vfrac_d(vfrac_v.data(), NTRIAL, NMAT);
    DataBox sie_d(sie_v.data(), NTRIAL, NMAT);
    DataBox temp_d(temp_v.data(), NTRIAL, NMAT);
    DataBox press_d(press_v.data(), NTRIAL, NMAT);
    DataBox scratch_d_rt(scratch_v.data(), NTRIAL * nscratch_vars_rt);
    DataBox scratch_d_pt(scratch_v.data(), NTRIAL * nscratch_vars_pt);
    DataBox rho_hm(rho_vh.data(), NTRIAL, NMAT);
    DataBox vfrac_hm(vfrac_vh.data(), NTRIAL, NMAT);
    DataBox sie_hm(sie_vh.data(), NTRIAL, NMAT);
    DataBox temp_hm(temp_vh.data(), NTRIAL, NMAT);
    DataBox press_hm(press_vh.data(), NTRIAL, NMAT);
    DataBox scratch_hm_rt(scratch_vh.data(), NTRIAL * nscratch_vars_rt);
    DataBox scratch_hm_pt(scratch_vh.data(), NTRIAL * nscratch_vars_pt);
#else
    DataBox rho_d(NTRIAL, NMAT);
    DataBox vfrac_d(NTRIAL, NMAT);
    DataBox sie_d(NTRIAL, NMAT);
    DataBox temp_d(NTRIAL, NMAT);
    DataBox press_d(NTRIAL, NMAT);
    DataBox scratch_d_rt(NTRIAL, nscratch_vars_rt);
    DataBox scratch_d_pt(NTRIAL, nscratch_vars_pt);
    DataBox rho_hm = rho_d.slice(2, 0, NTRIAL);
    DataBox vfrac_hm = vfrac_d.slice(2, 0, NTRIAL);
    DataBox sie_hm = sie_d.slice(2, 0, NTRIAL);
    DataBox temp_hm = temp_d.slice(2, 0, NTRIAL);
    DataBox press_hm = press_d.slice(2, 0, NTRIAL);
    DataBox scratch_hm_rt = scratch_d_rt.slice(2, 0, NTRIAL);
    DataBox scratch_hm_pt = scratch_d_pt.slice(2, 0, NTRIAL);
    int hist_vh[HIST_SIZE];
    int *hist_d = hist_vh;
#endif

    // setup state
    srand(time(NULL));
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

#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::View<int, atomic_view> nsuccess_d("n successes");
#else
    PortableMDArray<int> nsuccess_d(&nsuccess, 1);
#endif

    auto start_rt = std::chrono::high_resolution_clock::now();
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif
    std::cout << "Starting PTE RT with " << NTRIAL << " trials." << std::endl;
    portableFor(
        "PTE!", 0, NTRIAL, PORTABLE_LAMBDA(const int &t) {
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

          Real *scratch_rt = &scratch_d_rt(t * nscratch_vars_rt);

          auto method =
              PTESolverRhoT<EOSAccessor, Indexer2D<decltype(rho_d)>, decltype(lambda)>(
                  NMAT, eos, 1.0, sie_tot, rho, vfrac, sie, temp, press, lambda,
                  scratch_rt, Tguess);
          auto status = PTESolver(method);

          if (status.converged) {
            nsuccess_d() += 1;
          }
          hist_d[method.Niter()] += 1;
        });
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif
    auto stop_rt = std::chrono::high_resolution_clock::now();
    auto sum_time_rt =
        std::chrono::duration_cast<std::chrono::microseconds>(stop_rt - start_rt);

    auto start_pt = std::chrono::high_resolution_clock::now();
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif
    std::cout << "Starting PTE PT with " << NTRIAL << " trials." << std::endl;
    portableFor(
        "PTE!", 0, NTRIAL, PORTABLE_LAMBDA(const int &t) {
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

          Real *scratch_pt = &scratch_d_pt(t * nscratch_vars_pt);

          auto method2 =
              PTESolverPT<EOSAccessor, Indexer2D<decltype(rho_d)>, decltype(lambda)>(
                  NMAT, eos, 1.0, sie_tot, rho, vfrac, sie, temp, press, lambda,
                  scratch_pt, Tguess);
          auto status = PTESolver(method2);

          if (status.converged) {
            nsuccess_d() += 1;
          }
          hist_d[method2.Niter()] += 1;
        });
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif
    auto stop_pt = std::chrono::high_resolution_clock::now();
    auto sum_time_pt =
        std::chrono::duration_cast<std::chrono::microseconds>(stop_pt - start_pt);

#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::deep_copy(nsuccess, nsuccess_d);
    Kokkos::deep_copy(hist_vh, hist_d);
#endif

    Real milliseconds_rt = sum_time_rt.count() / 1e3;
    Real milliseconds_pt = sum_time_pt.count() / 1e3;

    std::cout << "Finished " << NTRIAL << " RT solves in " << milliseconds_rt
              << " milliseconds" << std::endl;
    std::cout << "\tSolves/ms = " << NTRIAL / milliseconds_rt << std::endl;
    std::cout << "Finished " << NTRIAL << " PT solves in " << milliseconds_pt
              << " milliseconds" << std::endl;
    std::cout << "\tSolves/ms = " << NTRIAL / milliseconds_pt << std::endl;
    std::cout << "Success: " << nsuccess << "   Failure: " << NTRIAL - nsuccess
              << std::endl;
    std::cout << "Histogram:\n"
              << "iters\tcount\n"
              << "----------------------\n";
    for (int i = 0; i < HIST_SIZE; ++i) {
      std::cout << i << "\t" << hist_vh[i] << "\n";
    }
    std::cout << std::endl;
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif

  // poor-man's ctest integration
  return (nsuccess >= 0.5 * NTRIAL) ? 0 : 1;
}
