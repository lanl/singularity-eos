//------------------------------------------------------------------------------
// © 2025. Triad National Security, LLC. All rights reserved.  This
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

// SCENARIO("Implementation of KPT: From out from hydrocode to in to hydrocode") {
// 1:  Massfractions, Internal energy, and density at time t0
// 2:  Calculate full state, at time t0
// 3:  This is given back to hydrocode for updated Internal energy, and density at time
// t1
// 4:  Update mass fractions: Massfractions at t0, full state at t0, Internal energy,
// and density at time t1
//
// 1 must be taken from real run. Use test_pte_3phases for a start.
// 1->2 multiphase EOS Initialized and used: test_pte_3phases
// 1->2 PTE solver parameters given.
// 2 should be given by PTE solver. Compare to what real run gives.
//   and calculate gibbs free energy, and compare to real run.
// 3 Internal energy, and density at time t1 from real run
// 3->4 update model and update method needs to be initialized and used.
// 4 KPT model gives this. Compare to real run

// 1, 1->2, 2 need to be input here. Check in pte-solve testing
// Use pte_test_5phaseSesameSn.hpp and test_pte_3phase.cpp
// This is where we should try to use the fancy skipping phases thing, Indexers.
// Indexers are in pte_test_utils.hpp
// It seems better to actually start from test_pte_3phase.cpp since it has a main.
// Try to figure out if it can be used together with catch2.
// Do separate header files for 1->2 and 3->4 but include both.

#include <chrono>
#include <iostream>
#include <memory>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>
#include <spiner/databox.hpp>
#include <test/pte_test_utils.hpp>

#include <singularity-eos/closure/kinetic_phasetransition_methods.hpp>
#include <singularity-eos/closure/kinetic_phasetransition_models.hpp>
#include <singularity-eos/closure/kinetic_phasetransition_utils.hpp>

#include <test/kpt_full_circle_test.hpp>

using namespace kpt_full_circle_test;
using namespace pte_test_3phaseSesameSn;

using DataBox = Spiner::DataBox<Real>;

using singularity::PTESolverRhoT;
using singularity::PTESolverRhoTRequiredScratch;

using singularity::LogMaxTimeStep;
using singularity::LogRatesCGModel;
using singularity::SmallStepMFUpdate;
using singularity::SmallStepMFUpdateR;
using singularity::SortGibbs;

int main(int argc, char *argv[]) {

  int nsuccess = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
#endif
  // JMM: EOSPAC tests do not work on device.
  if constexpr (PortsOfCall::EXECUTION_IS_HOST) {
    // EOS
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::View<EOS *> eos_v("eos", NMAT);
    Kokkos::View<Real *> CGw_v("mfupdate_CGmodelparamw", NMAT*NMAT);
    Kokkos::View<Real *> CGb_v("mfupdate_CGmodelparamb", NMAT*NMAT);
    auto eos_hv = Kokkos::create_mirror_view(eos_v);
    auto CGw_hv = Kokkos::create_mirror_view(CGw_v);
    auto CGb_hv = Kokkos::create_mirror_view(CGb_v);
#else
    std::vector<EOS> eos_vec(NMAT);
    std::vector<Real> CGw_vec(NMAT*NMAT);
    std::vector<Real> CGb_vec(NMAT*NMAT);
    PortableMDArray<EOS> eos_hv(eos_vec.data(), NMAT);
    PortableMDArray<EOS> eos_v(eos_vec.data(), NMAT);
    auto CGw_hv = CGw_vec.data();
    auto CGb_hv = CGb_vec.data();
    auto CGw_v = CGw_vec.data();
    auto CGb_v = CGb_vec.data();
#endif
    // set-eos and set_mfupdate_params according to material_5phaseSesameSn
    std::cout << "Load EOS tables/data for " << NMAT << " phases" << std::endl;

    set_eos(NMAT, PHASES, eos_hv.data());

    std::cout << "Load mass fraction update model parameters for " << NMAT*(NMAT-1)/2 << " phase transitions" << std::endl;

    set_mfupdate_params(NMAT,PHASES,CGw_hv,CGb_hv);

    std::cout << "Mass fraction update model parameters w and b are in the materials file together with the eos they are determined for. " << std::endl << std::endl;

#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::deep_copy(eos_v, eos_hv);
    Kokkos::deep_copy(CGw_v, CGw_hv);
    Kokkos::deep_copy(CGb_v, CGb_hv);
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
    RView mfrac_v("mfrac", NPTS);
    RView sie_v("sie", NPTS);
    RView temp_v("temp", NPTS);
    RView press_v("press", NPTS);
    RView gibbsrt_v("gibbsRT", NPTS);
    RView gibbsre_v("gibbsRE", NPTS);
    RView scratch_v("scratch", NTRIAL * nscratch_vars);
    Kokkos::View<int *, atomic_view> hist_d("histogram", HIST_SIZE);
    auto rho_vh = Kokkos::create_mirror_view(rho_v);
    auto vfrac_vh = Kokkos::create_mirror_view(vfrac_v);
    auto mfrac_vh = Kokkos::create_mirror_view(mfrac_v);
    auto sie_vh = Kokkos::create_mirror_view(sie_v);
    auto temp_vh = Kokkos::create_mirror_view(temp_v);
    auto press_vh = Kokkos::create_mirror_view(press_v);
    auto gibbsrt_vh = Kokkos::create_mirror_view(gibbsrt_v);
    auto gibbsre_vh = Kokkos::create_mirror_view(gibbsre_v);
    auto scratch_vh = Kokkos::create_mirror_view(scratch_v);
    auto hist_vh = Kokkos::create_mirror_view(hist_d);
    DataBox rho_d(rho_v.data(), NTRIAL, NMAT);
    DataBox vfrac_d(vfrac_v.data(), NTRIAL, NMAT);
    DataBox mfrac_d(mfrac_v.data(), NTRIAL, NMAT);
    DataBox sie_d(sie_v.data(), NTRIAL, NMAT);
    DataBox temp_d(temp_v.data(), NTRIAL, NMAT);
    DataBox press_d(press_v.data(), NTRIAL, NMAT);
    DataBox gibbsrt_d(gibbsrt_v.data(), NTRIAL, NMAT);
    DataBox gibbsre_d(gibbsre_v.data(), NTRIAL, NMAT);
    DataBox scratch_d(scratch_v.data(), NTRIAL * nscratch_vars);
    DataBox rho_hm(rho_vh.data(), NTRIAL, NMAT);
    DataBox vfrac_hm(vfrac_vh.data(), NTRIAL, NMAT);
    DataBox mfrac_hm(mfrac_vh.data(), NTRIAL, NMAT);
    DataBox sie_hm(sie_vh.data(), NTRIAL, NMAT);
    DataBox temp_hm(temp_vh.data(), NTRIAL, NMAT);
    DataBox press_hm(press_vh.data(), NTRIAL, NMAT);
    DataBox gibbsrt_hm(gibbsrt_vh.data(), NTRIAL, NMAT);
    DataBox gibbsre_hm(gibbsre_vh.data(), NTRIAL, NMAT);
    DataBox scratch_hm(scratch_vh.data(), NTRIAL * nscratch_vars);
#else
    DataBox rho_d(NTRIAL, NMAT);
    DataBox vfrac_d(NTRIAL, NMAT);
    DataBox mfrac_d(NTRIAL, NMAT);
    DataBox sie_d(NTRIAL, NMAT);
    DataBox temp_d(NTRIAL, NMAT);
    DataBox press_d(NTRIAL, NMAT);
    DataBox gibbsrt_d(NTRIAL, NMAT);
    DataBox gibbsre_d(NTRIAL, NMAT);
    DataBox scratch_d(NTRIAL, nscratch_vars);
    DataBox rho_hm = rho_d.slice(2, 0, NTRIAL);
    DataBox vfrac_hm = vfrac_d.slice(2, 0, NTRIAL);
    DataBox mfrac_hm = mfrac_d.slice(2, 0, NTRIAL);
    DataBox sie_hm = sie_d.slice(2, 0, NTRIAL);
    DataBox temp_hm = temp_d.slice(2, 0, NTRIAL);
    DataBox press_hm = press_d.slice(2, 0, NTRIAL);
    DataBox gibbsrt_hm = gibbsrt_d.slice(2, 0, NTRIAL);
    DataBox gibbsre_hm = gibbsre_d.slice(2, 0, NTRIAL);
    DataBox scratch_hm = scratch_d.slice(2, 0, NTRIAL);
    int hist_vh[HIST_SIZE];
    int *hist_d = hist_vh;
#endif

    // setup state
    std::cout
        << "Setup Initial states (from Density, Mass fractions, and Internal energy)."
        << std::endl << std::endl;

    srand(time(NULL));
    for (int n = 0; n < NTRIAL; n++) {
      Indexer2D<decltype(rho_hm)> r(n, rho_hm);
      Indexer2D<decltype(vfrac_hm)> vf(n, vfrac_hm);
      Indexer2D<decltype(mfrac_hm)> mf(n, mfrac_hm);
      Indexer2D<decltype(sie_hm)> e(n, sie_hm);
      set_trial_state(n, r, vf, mf, e);
    }
    for (int i = 0; i < HIST_SIZE; ++i) {
      hist_vh[i] = 0;
    }
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::deep_copy(rho_v, rho_vh);
    Kokkos::deep_copy(vfrac_v, vfrac_vh);
    Kokkos::deep_copy(mfrac_v, mfrac_vh);
    Kokkos::deep_copy(sie_v, sie_vh);
    Kokkos::deep_copy(temp_v, temp_vh);
    Kokkos::deep_copy(press_v, press_vh);
    Kokkos::deep_copy(gibbsrt_v, gibbsrt_vh);
    Kokkos::deep_copy(gibbsre_v, gibbsre_vh);
    Kokkos::deep_copy(hist_d, hist_vh);
#endif
    std::cout << "Initial states for " << NTRIAL << " trials set." << std::endl << std::endl;

#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::View<int, atomic_view> nsuccess_d("n successes");
#else
    PortableMDArray<int> nsuccess_d(&nsuccess, 1);
#endif

    auto start = std::chrono::high_resolution_clock::now();
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif
    std::cout << "The trials have initial, time=t0, states: " << std::endl;

    for (int n = 0; n < NTRIAL; n++) {
      Indexer2D<decltype(rho_hm)> r(n, rho_hm);
      Indexer2D<decltype(vfrac_hm)> vf(n, vfrac_hm);
      Indexer2D<decltype(sie_hm)> e(n, sie_hm);
      printinput(n, r, vf);
    }

    std::cout << "Starting PTE solver for " << NTRIAL << " trials." << std::endl
              << std::endl;

    portableFor(
        "PTE in Full Circle test", 0, NTRIAL, PORTABLE_LAMBDA(const int &t) {
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
          if (t == 0) {
            printf("Tguess %.14e\n", Tguess);
          }

          auto method =
              PTESolverRhoT<EOSAccessor, Indexer2D<decltype(rho_d)>, decltype(lambda)>(
                  NMAT, eos, 1.0, sie_tot, rho, vfrac, sie, temp, press, lambda,
                  &scratch_d(t * nscratch_vars), Tguess);
          auto status = PTESolver(method);
          if (status.converged) {
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
    std::cout << "The PTE solver converged: " << nsuccess << ", did NOT converge: " << NTRIAL - nsuccess
              << std::endl;
    std::cout << "Histogram:\n"
              << "iters\tcount\n"
              << "----------------------\n";
    for (int i = 0; i < HIST_SIZE; ++i) {
      std::cout << i << "\t" << hist_vh[i] << "\n";
    }
    std::cout << std::endl;

    std::cout << "Results from the PTE solver are: (time = t0) " << std::endl;

    for (int n = 0; n < NTRIAL; n++) {
      Indexer2D<decltype(rho_hm)> rho(n, rho_hm);
      Indexer2D<decltype(vfrac_hm)> vfrac(n, vfrac_hm);
      Indexer2D<decltype(sie_hm)> sie(n, sie_hm);
      Indexer2D<decltype(temp_hm)> temp(n, temp_hm);
      Indexer2D<decltype(press_hm)> press(n, press_hm);
      printresults(n, rho, vfrac, sie, press, temp);
    }

    std::cout << "t0 results from the PTE solver handed to host code." << std::endl;
    std::cout << "We can also give the host code a suggestion for time step, dt, to t1 = dt + t0, see below." << std::endl;
    std::cout << "Hostcode gives back new Internal energy and density for time t1." << std::endl << std::endl;

#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif
    std::cout << "We still have to update mass/volume fractions to time t1." << std::endl << std::endl;

    std::cout << "Getting t0 Gibbs Free energies, needed for the update model." << std::endl;
    
    for (int n = 0; n < NTRIAL; n++) {
      Indexer2D<decltype(rho_hm)> rho(n, rho_hm);
      Indexer2D<decltype(sie_hm)> sie(n, sie_hm);
      Indexer2D<decltype(temp_hm)> temp(n, temp_hm);
      Indexer2D<decltype(gibbsrt_hm)> gibbsrt(n, gibbsrt_hm);
      Indexer2D<decltype(gibbsre_hm)> gibbsre(n, gibbsre_hm);
      for (int i = 0; i < NMAT; i++) {
        gibbsrt[i] = eos[i].GibbsFreeEnergyFromDensityTemperature(rho[i], temp[i]);
        gibbsre[i] = eos[i].GibbsFreeEnergyFromDensityInternalEnergy(rho[i], sie[i]);
      }
      gibbsprintresults(n, rho, temp, sie, gibbsrt, gibbsre);
    }
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif // PORTABILITY_STRATEGY_KOKKOS

    std::cout << "Updating mass fractions/volume fractions at time = t0 to time = t1." << std::endl << std::endl;

    constexpr int mnum=NMAT*(NMAT-1)/2;
#ifdef PORTABILITY_STRATEGY_KOKKOS
    // Create device views for outputs and mirror those views on the host
    Kokkos::View<int[NMAT]> v_order("Gibbsorder");
    auto h_order = Kokkos::create_mirror_view(v_order);
    Kokkos::View<int[mnum]> v_rateorder("Rateorder");
    auto h_rateorder = Kokkos::create_mirror_view(v_rateorder);
#else
    // Create arrays for the outputs and then pointers to those arrays that
    // will be passed to the functions in place of the Kokkos views
    std::array<int, NMAT> h_order;
    std::array<int, mnum> h_rateorder;
    // Just alias the existing pointers
    auto v_order = h_order.data();
    auto v_rateorder = h_rateorder.data();
#endif // PORTABILITY_STRATEGY_KOKKOS


#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::View<Real[mnum]> v_logrates("LogRates");
    auto h_logrates = Kokkos::create_mirror_view(v_logrates);
    Kokkos::View<Real[mnum]> v_dgibbs("dgibbs");
    auto h_dgibbs = Kokkos::create_mirror_view(v_dgibbs);
    Kokkos::View<int[mnum]> v_fromto("fromto");
    auto h_fromto = Kokkos::create_mirror_view(v_fromto);
#else
    std::array<Real, mnum> h_logrates;
    auto v_logrates = h_logrates.data();
    std::array<Real, mnum> h_dgibbs;
    auto v_dgibbs = h_dgibbs.data();
    std::array<int, mnum> h_fromto;
    auto v_fromto = h_fromto.data();
#endif // PORTABILITY_STRATEGY_KOKKOS
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create Kokkos views on device for the input arrays
      // Create host-side mirrors
    Kokkos::View<Real[NMAT]> v_newmfs("newmassfractions");
    auto newmassfractions = Kokkos::create_mirror_view(v_newmfs);
    Kokkos::View<Real[mnum]> v_deltamfs("massfractiontransfers");
    auto h_deltamfs = Kokkos::create_mirror_view(v_deltamfs);
#else
      // Otherwise just create arrays to contain values and create pointers to
      // be passed to the functions in place of the Kokkos views
    std::array<Real, NMAT> newmassfractions;
    auto v_newmfs = newmassfractions.data();
    std::array<Real, mnum> h_deltamfs;
    auto v_deltamfs = h_deltamfs.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

//  do a verbose full update for each trial.

    Real logmts;

    for (int n = 0; n < NTRIAL; n++) {
      Indexer2D<decltype(gibbsre_hm)> gibbsre(n, gibbsre_hm);
      Indexer2D<decltype(vfrac_hm)> vfrac(n, vfrac_hm);
      Indexer2D<decltype(mfrac_hm)> mfrac(n, mfrac_hm);

      std::cout << "---------------Trial " << n << "---------------" << std::endl << std::endl;

      for (int l = 0; l < NMAT; l++) {
        std::cout << "Phase: " << l
                  << "  initial (t0) mass fractions: " << mfrac[l]
                  << std::endl;
      }
      std::cout << std::endl;

      std::cout << "Using SortGibbs(num,gibbs,order) to give the num phases in order of largest to "
           "smallest Gibbs free energy" << std::endl << std::endl;
      
      SortGibbs(NMAT, gibbsre, v_order);

      std::cout << "Order obtained with SortGibbs, largest to smallest Gibbs: " << std::endl;
      for (int l = 0; l < NMAT; l++) {
        std::cout << "Phase " << v_order[l] << " Gibbs: " << gibbsre[v_order[l]]
                  << std::endl;
      }
      std::cout << std::endl;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
      Kokkos::deep_copy(h_order, v_order);
#endif // PORTABILITY_STRATEGY_KOKKOS

//    array_compare(num, gibbs, phaseorder, order, order_true, "Gibbs", "phase");
//
// create deltagibbs and fromphasetophase
      std::cout << "Rates are calculated from the difference in Gibbs free energy between phase i and k: " << std::endl;
      int ik = 0;
      for (int i = 0; i < NMAT - 1; i++) {
        for (int k = NMAT - 1; k > i; k--) {
          v_dgibbs[ik] = gibbsre[v_order[i]] - gibbsre[v_order[k]];
	  std::cout << "DeltaG between phase " << v_order[i] << " and " << v_order[k] << " is: " << v_dgibbs[ik] << std::endl;
          ik++;
        }
      }
      std::cout << std::endl;

#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
      Kokkos::deep_copy(h_dgibbs, v_dgibbs);
#endif // PORTABILITY_STRATEGY_KOKKOS

      std::cout << "Using LogRij(w,b,num,gibbs,order,fromto) to get the logarithm of the rates from num phases i to j. " << std::endl << std::endl;

      LogRatesCGModel(CGw_v, CGb_v, NMAT, gibbsre, v_order, v_logrates, v_fromto);

      std::cout << "LogRates obtained with LogRatesCGModel: " << std::endl;
      for (int l = 0; l < mnum; l++) {
          std::cout << "From phase i to phase k, ik ( x means 0x): " << v_fromto[l]
                    << "   LogRik: " << v_logrates[l] << std::endl;
      }
      std::cout << std::endl;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
      Kokkos::deep_copy(h_fromto, v_fromto);
      Kokkos::deep_copy(h_logrates, v_logrates);
#endif // PORTABILITY_STRATEGY_KOKKOS
     
//    array_compare(mnum, dgibbs, fromto, h_logrates, logrates_true, "DeltaGibbs",
//                          "FromTo");

      std::cout << "Optional: Give host code a suggestion for time step. " << std::endl << std::endl;
      std::cout << "A LogMaxTimeStep(num,order,logrates) lookup is performed" << std::endl;

      logmts = LogMaxTimeStep(NMAT, mfrac, v_order, v_logrates);
//      isClose(logmts, logmts_true);

      std::cout << "Log(MaxTimeStep) from rates obtained with CGModel: " << logmts
                << std::endl;
      std::cout << "This timestep would give a mass transfer of Rdt*(mass fractions at t0), " << std::endl;
      for (int l = 0; l < mnum; l++) {
        std::cout << "with Rdt = " << std::exp(logmts + v_logrates[l]) 
		  << ", from phase i to phase k, ik ( x means 0x): " << v_fromto[l] << std::endl;
      }
      std::cout << std::endl;

// a time step of 10^(-11) s
      logmts = -25.3284;
      
      std::cout << "Using a 10^(-11)s timestep to update massfractions with SmallStepMFUpdate" << std::endl << std::endl; 

      SmallStepMFUpdate(logmts, NMAT, mfrac, v_order, v_logrates,
                         v_deltamfs, v_newmfs);
      for (int l = 0; l < NMAT; l++) {
        std::cout << "Phase: " << v_order[l]
                  << "  initial mass fractions: " << mfrac[v_order[l]]
                  << std::endl;
        std::cout << "         "
                  << "   final mass fractions: " << v_newmfs[v_order[l]] << "          (" << flag_out_lambda[v_order[l]][n] << ")"
                  << std::endl;
      }
      std::cout << std::endl;
      for (int l = 0; l < mnum; l++) {
        std::cout << "From phase i to phase k, ik ( x means 0x): " << v_fromto[l]
                  << "   Mass fraction transfer: " << v_deltamfs[l] << std::endl;
      }
      std::cout << std::endl;

      std::cout << "Using a 10^(-11)s timestep to update massfractions with SmallStepMFUpdateR" << std::endl << std::endl; 

      SmallStepMFUpdateR(logmts, NMAT, mfrac, v_order, v_logrates,
                         v_deltamfs, v_newmfs);
      for (int l = 0; l < NMAT; l++) {
        std::cout << "Phase: " << v_order[l]
                  << "  initial mass fractions: " << mfrac[v_order[l]]
                  << std::endl;
        std::cout << "         "
                  << "   final mass fractions: " << v_newmfs[v_order[l]] << "          (" << flag_out_lambda[v_order[l]][n] << ")"
                  << std::endl;
      }
      std::cout << std::endl;
      for (int l = 0; l < mnum; l++) {
        std::cout << "From phase i to phase k, ik ( x means 0x): " << v_fromto[l]
                  << "   Mass fraction transfer: " << v_deltamfs[l] << std::endl;
      }
      std::cout << std::endl;

#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
      Kokkos::deep_copy(newmassfractions, v_newmfs);
      Kokkos::deep_copy(h_deltamfs, v_deltamfs);
#endif // PORTABILITY_STRATEGY_KOKKOS
//      array_compare(num, phaseorder, massfractions, newmassfractions,
//                    newmassfractions_true, "Phase", "old massfractions");
//      array_compare(mnum, dgibbs, fromto, h_deltamfs, deltamfs_true,
//                    "DeltaGibbs", "FromTo");
//
      std::cout << "Using SortGibbs(num,rates,order) to give the num phase transitions in order of largest to "
           "smallest rates" << std::endl << std::endl;

      SortGibbs(mnum, v_logrates, v_rateorder);

      std::cout << "Order obtained with SortGibbs, largest to smallest rates: " << std::endl;
      for (int l = 0; l < mnum; l++) {
        std::cout << "Transition " << v_rateorder[l] << "    from/to phases: " << v_fromto[v_rateorder[l]] << " Rate: " << v_logrates[v_rateorder[l]]
                  << std::endl;
      }
      std::cout << std::endl;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
      Kokkos::deep_copy(h_rateorder, v_rateorder);
#endif // PORTABILITY_STRATEGY_KOKKOS

    }

  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif

  // poor-man's ctest integration
  if constexpr (PortsOfCall::EXECUTION_IS_HOST) {
    return (nsuccess >= 0.5 * NTRIAL) ? 0 : 1;
  } else {
    return 0;
  }
}
