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

#include <chrono>
#include <iostream>
#include <memory>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/singularity_eos.hpp>
#include <spiner/databox.hpp>

constexpr int NMAT = 3;
constexpr int NTRIAL = 100;
constexpr int NPTS = NTRIAL * NMAT;
constexpr int HIST_SIZE = 10;

using namespace singularity;
using Spiner::DataBox;

template <typename T>
class LinearIndexer {
 public:
  PORTABLE_FUNCTION LinearIndexer() = default;
  LinearIndexer(const T &t) : data_(t) {}
  PORTABLE_INLINE_FUNCTION
  auto &operator[](const int i) const { return data_(i); }

 private:
  T data_;
};

template <typename T>
class Indexer2D {
 public:
  PORTABLE_FUNCTION Indexer2D() = default;
  PORTABLE_FUNCTION Indexer2D(const int j, const T &t) : j_(j), data_(t) {}
  Indexer2D(const int j, const T &&t) = delete; // prevents r-value binding
  PORTABLE_INLINE_FUNCTION
  auto &operator[](const int i) const { return data_(j_, i); }

 private:
  const int j_;
  const T &data_;
};

template <typename EOSIndexer>
void set_eos(EOSIndexer &&eos) {
  eos[0] = Gruneisen(394000.0, 1.489, 0.0, 0.0, 2.02, 0.47, 8.93, 297.0, 1.0e6, 0.383e7);
  eos[1] = DavisReactants(1.890, 4.115e10, 1.0e6, 297.0, 1.8e5, 4.6, 0.34, 0.56, 0.0,
                          0.4265, 0.001074e10);
  eos[2] =
      DavisProducts(0.798311, 0.58, 1.35, 2.66182, 0.75419, 3.2e10, 0.001072e10, 0.0);
}

template <typename RealIndexer, typename EOSIndexer>
void set_state(RealIndexer &&rho, RealIndexer &&vfrac, RealIndexer &&sie,
               RealIndexer &&temp, EOSIndexer &&eos) {

  rho[0] = 8.93;
  rho[1] = 1.89;
  rho[2] = 2.5;

  Real vsum = 0.;
  for (int i = 0; i < NMAT; i++) {
    temp[i] = 600.0;
    sie[i] = eos[i].InternalEnergyFromDensityTemperature(rho[i], temp[i]);
    vfrac[i] = rand() / (1.0 * RAND_MAX);
    vsum += vfrac[i];
  }

  for (int i = 0; i < NMAT; i++)
    vfrac[i] *= 1.0 / vsum;
}

#ifdef PORTABILITY_STRATEGY_KOKKOS
void run_sg_get_eos_tests() {
  // initialize inputs outputs
  static constexpr const double ev2k = 1.160451930280894026e4;
  EOS eoss[NMAT];
  set_eos(eoss);
  // set volume fractions to equipartition
  // all variables needed to check the pte solve
  Real mfrac[NMAT];
  Real P_true = 5.e10, T_true = 800.0;
  Real T_true_ev = T_true / ev2k;
  Real mass_tot = 0.0, sie_tot_true = 0.0, v_true = 1.0, rho_tot, spvol;
  Real pmax, bmod, dpde, cv;
  mfrac[0] = 0.67;
  mfrac[1] = 0.14;
  mfrac[2] = 0.19;
  for (int i = 0; i < NMAT; ++i) {
    mass_tot += mfrac[i];
  }
  rho_tot = mass_tot / v_true;
  spvol = 1.0 / rho_tot;
  int cell_offset = 1;
  int eos_offset[NMAT];
  for (int m = 0; m < NMAT; ++m) {
    eos_offset[m] = m + 1;
  }
  // do a rho/e(P,T) solve
  Real vfrac_true[NMAT], sie_true[NMAT];
  get_sg_eos(NMAT, 1, 1, -1, eos_offset, eoss, &cell_offset, &P_true, &pmax, &v_true,
             &spvol, &sie_tot_true, &T_true_ev, &bmod, &dpde, &cv, mfrac, vfrac_true,
             sie_true, nullptr, nullptr, nullptr);
  Real sie_tot_check = 0.0;
  for (int m = 0; m < NMAT; ++m) {
    const Real r_m = mfrac[m] / vfrac_true[m];
    const Real p_eos = eoss[m].PressureFromDensityInternalEnergy(r_m, sie_true[m]);
    const Real p_resid = std::abs(P_true - p_eos) / P_true;
    const Real t_eos = eoss[m].TemperatureFromDensityInternalEnergy(r_m, sie_true[m]);
    const Real t_resid = std::abs(T_true - t_eos) / T_true;
    if (p_resid > 1.e-5) printf("p(%i) wrong: %e %e\n", m, P_true, p_eos);
    if (t_resid > 1.e-5) printf("t(%i) wrong: %e %e\n", m, T_true, t_eos);
    sie_tot_check += sie_true[m] * mfrac[m];
  }
  sie_tot_check /= mass_tot;
  printf("sie: %e %e \n", sie_tot_true, sie_tot_check);
  // obtain converged and consistent PTE solution
  // do rho-T input solve
  Real p_check, vfrac_check[NMAT], sie_check[NMAT];
  get_sg_eos(NMAT, 1, 1, -3, eos_offset, eoss, &cell_offset, &p_check, &pmax, &v_true,
             &spvol, &sie_tot_check, &T_true_ev, &bmod, &dpde, &cv, mfrac, vfrac_check,
             sie_check, nullptr, nullptr, nullptr);
  printf("p: %e %e \n", P_true, p_check);
  printf("s: %e %e \n", sie_tot_true, sie_tot_check);
  for (int m = 0; m < NMAT; ++m) {
    printf("m=%i\n", m);
    printf("  rho: %e %e \n", mfrac[m] / vfrac_true[m], mfrac[m] / vfrac_check[m]);
    printf("  sie: %e %e \n", sie_true[m], sie_check[m]);
    // sie_tot_check += sie_true[m] * mfrac[m];
  }
  // ensure outputs are the same
  // do rho-P input solve
  Real t_check;
  get_sg_eos(NMAT, 1, 1, -2, eos_offset, eoss, &cell_offset, &P_true, &pmax, &v_true,
             &spvol, &sie_tot_check, &t_check, &bmod, &dpde, &cv, mfrac, vfrac_check,
             sie_check, nullptr, nullptr, nullptr);
  printf("t: %e %e \n", T_true, t_check * ev2k);
  printf("s: %e %e \n", sie_tot_true, sie_tot_check);
  for (int m = 0; m < NMAT; ++m) {
    printf("m=%i\n", m);
    printf("  rho: %e %e \n", mfrac[m] / vfrac_true[m], mfrac[m] / vfrac_check[m]);
    printf("  sie: %e %e \n", sie_true[m], sie_check[m]);
    // sie_tot_check += sie_true[m] * mfrac[m];
  }

  // ensure outputs are the same
  return;
}
#endif

int main(int argc, char *argv[]) {

  int nsuccess = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
#endif
  {
    // run the various solver tests:
    // since this uses the sg_get_eos, only run
    // if kokkos enable since that function requires it
    // to run the solvers
#ifdef PORTABILITY_STRATEGY_KOKKOS
    run_sg_get_eos_tests();
#endif // PORTABILITY_STRATEGY_KOKKOS
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
    set_eos(eos_h);

#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::deep_copy(eos_v, eos_hv);
#endif

    using EOSAccessor = LinearIndexer<decltype(eos_v)>;
    EOSAccessor eos(eos_v);

    // scratch required for PTE solver
    auto nscratch_vars = PTESolverRhoTRequiredScratch(NMAT);

    // state vars
#ifdef PORTABILITY_STRATEGY_KOKKOS
    using RView = Kokkos::View<Real *>;
    using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
    RView rho_v("rho", NPTS);
    RView vfrac_v("vfrac", NPTS);
    RView sie_v("sie", NPTS);
    RView temp_v("temp", NPTS);
    RView press_v("press", NPTS);
    RView scratch_v("scratch", NTRIAL * nscratch_vars);
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
    DataBox scratch_d(scratch_v.data(), NTRIAL * nscratch_vars);
    DataBox rho_hm(rho_vh.data(), NTRIAL, NMAT);
    DataBox vfrac_hm(vfrac_vh.data(), NTRIAL, NMAT);
    DataBox sie_hm(sie_vh.data(), NTRIAL, NMAT);
    DataBox temp_hm(temp_vh.data(), NTRIAL, NMAT);
    DataBox press_hm(press_vh.data(), NTRIAL, NMAT);
    DataBox scratch_hm(scratch_vh.data(), NTRIAL * nscratch_vars);
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

    auto start = std::chrono::high_resolution_clock::now();
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif
    std::cout << "Starting PTE with " << NTRIAL << " trials." << std::endl;
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

          auto method =
              PTESolverRhoT<EOSAccessor, Indexer2D<decltype(rho_d)>, decltype(lambda)>(
                  NMAT, eos, 1.0, sie_tot, rho, vfrac, sie, temp, press, lambda,
                  &scratch_d(t * nscratch_vars));
          bool success = PTESolver(method);
          //	  if (t == 0) {
          //	    printf("\n\n---\n");
          //	    for (int m = 0; m < NMAT; ++m) {
          //	      printf("r: %e\n", rho[m]);
          //	      printf("v: %e\n", vfrac[m]);
          //	      printf("t: %e\n", temp[m]);
          //	      printf("p: %e\n", press[m]);
          //	    }
          //	  }
          if (success) {
            nsuccess_d() += 1;
          }
          hist_d[method.Niter()] += 1;
        });
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif
    auto stop = std::chrono::high_resolution_clock::now();
    auto sum_time = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::deep_copy(nsuccess, nsuccess_d);
    Kokkos::deep_copy(hist_vh, hist_d);
#endif

    Real milliseconds = sum_time.count() / 1e3;

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
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif

  // poor-man's ctest integration
  if (nsuccess >= 0.5 * NTRIAL) {
    return 0; // exit success
  } else {
    return 1; // exit failure
  }
}
