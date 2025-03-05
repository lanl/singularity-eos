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

// TODO: Need to make this working on GPUs
//#ifdef PORTABILITY_STRATEGY_NONE
#ifndef PORTABILITY_STRATEGY_KOKKOS

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
#include <singularity-eos/closure/kinetic_phasetransition_methods.hpp>
#include <singularity-eos/closure/kinetic_phasetransition_models.hpp>
#include <singularity-eos/closure/kinetic_phasetransition_utils.hpp>
#include <test/eos_unit_test_helpers.hpp>

using singularity::LogMaxTimeStep;
using singularity::LogRatesCGModel;
using singularity::SmallStepMFUpdate;
using singularity::SortGibbs;

SCENARIO("Testing the KPT framework") {
  GIVEN("Parameters for a 5 phase system") {
    constexpr Real w[25] = {0.,  0.8510E+01, 20., 20., 20., 0.8510E+01, 0.,  20., 20.,
                            20., 20.,        20., 0.,  20., 20.,        20., 20., 20.,
                            0.,  20.,        20., 20., 20., 20.,        0.};

    constexpr Real b[25] = {0.,         0.3060E-04, 0.1000E-05, 0.1000E-05, 0.1000E-05,
                            0.3060E-04, 0.,         0.1000E-05, 0.1000E-05, 0.1000E-05,
                            0.1000E-05, 0.1000E-05, 0.,         0.1000E-05, 0.1000E-05,
                            0.1000E-05, 0.1000E-05, 0.1000E-05, 0.,         0.1000E-05,
                            0.1000E-05, 0.1000E-05, 0.1000E-05, 0.1000E-05, 0.};

    GIVEN("Gibbs free energy") {
      constexpr int num = 5;
      constexpr int mnum = num * (num - 1) / 2;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create Kokkos views on device for the input arrays
      // Create host-side mirrors of the inputs
      Kokkos::View<Real[num]> v_gibbs("gibbs");
      auto gibbs = Kokkos::create_mirror_view(v_gibbs);
      // order phases from high to low gibbs
      Kokkos::View<int[num]> v_order("order");
      auto order = Kokkos::create_mirror_view(v_order);
#else
      // Otherwise just create arrays to contain values and create pointers to
      // be passed to the functions in place of the Kokkos views
      std::array<Real, num> gibbs;
      std::array<int, num> order;
      auto v_gibbs = gibbs.data();
      auto v_order = order.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

      // Populate the input arrays
      portableFor(
          "Initialize gibbs", 0, 1, PORTABLE_LAMBDA(int i) {
            v_gibbs[0] = 0.0292971;
            v_gibbs[1] = 0.0286937;
            v_gibbs[2] = 0.0287409;
            v_gibbs[3] = 1.0e99;
            v_gibbs[4] = 0.02940076;
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // copy the inputs
      Kokkos::deep_copy(gibbs, v_gibbs);
#endif // PORTABILITY_STRATEGY_KOKKOS

      // Gold standard values for a subset of lookups
      constexpr std::array<Real, mnum> logrates_true{
          1.00000000000000007e+210, 1.00000000000000007e+210, 1.00000000000000007e+210,
          1.00000000000000007e+210, 4.99943400447804946e+05,  4.35424707359967171e+05,
          1.07530324485867186e+04,  3.93959978909504514e+02,  3.09367756860215100e+05,
          2.23469012616621285e+03};
      constexpr std::array<int, num> phaseorder{0, 1, 2, 3, 4};
      constexpr std::array<int, num> order_true{3, 4, 0, 2, 1};

      WHEN("Using SortGibbs(num,gibbs,order) to give phases in order of largest to "
           "smallest gibbs") {
        SortGibbs(num, v_gibbs, v_order);

        std::cout << "Order obtained with SortGibbs, largest to smallest Gibbs: "
                  << std::endl;
        for (int l = 0; l < num; l++) {
          std::cout << "Phase " << v_order[l] << " Gibbs: " << v_gibbs[v_order[l]]
                    << std::endl;
        }
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(order, v_order);
#endif // PORTABILITY_STRATEGY_KOKKOS

        THEN("The returned order should be equal to the true order") {
          array_compare(num, gibbs, phaseorder, order, order_true, "Gibbs", "phase");
        }
        // create deltagibbs and fromphasetophase
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::View<Real[mnum]> v_dgibbs("dgibbs");
        auto dgibbs = Kokkos::create_mirror_view(v_dgibbs);
        Kokkos::View<int[mnum]> v_fromto("fromto");
        auto fromto = Kokkos::create_mirror_view(v_fromto);
#else
        std::array<Real, mnum> dgibbs;
        std::array<int, mnum> fromto;
        auto v_dgibbs = dgibbs.data();
        auto v_fromto = fromto.data();
#endif // PORTABILITY_STRATEGY_KOKKOS
        int ik = 0;
        for (int i = 0; i < num - 1; i++) {
          for (int k = num - 1; k > i; k--) {
            v_dgibbs[ik] = v_gibbs[v_order[i]] - v_gibbs[v_order[k]];
            ik++;
          }
        }
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::deep_copy(dgibbs, v_dgibbs);
#endif // PORTABILITY_STRATEGY_KOKKOS

#ifdef PORTABILITY_STRATEGY_KOKKOS
        // Create device views for outputs and mirror those views on the host
        Kokkos::View<Real[mnum]> v_logrates("LogRates");
        auto h_logrates = Kokkos::create_mirror_view(v_logrates);
#else
        // Create arrays for the outputs and then pointers to those arrays that
        // will be passed to the functions in place of the Kokkos views
        std::array<Real, mnum> h_logrates;
        // Just alias the existing pointers
        auto v_logrates = h_logrates.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

        WHEN("A LogRij(w,b,num,gibbs,order,fromto) lookup is performed") {
          LogRatesCGModel(w, b, num, v_gibbs, v_order, v_logrates, v_fromto);

          std::cout << "LogRates obtained with LogRatesCGModel: " << std::endl;
          for (int l = 0; l < mnum; l++) {
            std::cout << "From phase i to phase k, ik ( x means 0x): " << v_fromto[l]
                      << "   LogRik: " << v_logrates[l] << std::endl;
          }
#ifdef PORTABILITY_STRATEGY_KOKKOS
          Kokkos::fence();
          Kokkos::deep_copy(fromto, v_fromto);
          Kokkos::deep_copy(h_logrates, v_logrates);
#endif // PORTABILITY_STRATEGY_KOKKOS

          THEN("The returned LogRij(w,b,gibbs,order) should be equal to the true value") {
            array_compare(mnum, dgibbs, fromto, h_logrates, logrates_true, "DeltaGibbs",
                          "FromTo");
          }

          constexpr Real logmts_true = -2234.69;

#ifdef PORTABILITY_STRATEGY_KOKKOS
          // Create Kokkos views on device for the input arrays
          // Create host-side mirrors of the inputs
          Kokkos::View<Real[num]> v_massfractions("massfractions");
          auto massfractions = Kokkos::create_mirror_view(v_massfractions);
#else
          // Otherwise just create arrays to contain values and create pointers to
          // be passed to the functions in place of the Kokkos views
          std::array<Real, num> massfractions;
          auto v_massfractions = massfractions.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

          // Populate the input arrays
          portableFor(
              "Initialize mass fractions", 0, 1, PORTABLE_LAMBDA(int i) {
                v_massfractions[0] = 0.0;
                v_massfractions[1] = 0.5;
                v_massfractions[2] = 0.5;
                v_massfractions[3] = 0.0;
                v_massfractions[4] = 0.0;
              });
#ifdef PORTABILITY_STRATEGY_KOKKOS
          // copy the inputs
          Kokkos::deep_copy(massfractions, v_massfractions);
#endif // PORTABILITY_STRATEGY_KOKKOS

          Real logmts;

          WHEN("A LogMaxTimeStep(num,order,logrates) lookup is performed") {
            logmts = LogMaxTimeStep(num, v_massfractions, v_order, v_logrates);
            for (int l = 0; l < num; l++) {
              std::cout << "Phase: " << l
                        << "  initial mass fractions: " << v_massfractions[l]
                        << std::endl;
            }
            std::cout << "Log(MaxTimeStep) from Rates obtained with CGModel: " << logmts
                      << std::endl;
            for (int l = 0; l < mnum; l++) {
              std::cout << "From phase i to phase k, ik ( x means 0x): " << v_fromto[l]
                        << "   Max mass fraction transfer: "
                        << std::exp(logmts + v_logrates[l]) << std::endl;
            }

            THEN("The returned Log(MaxTimeStep) should be equal to the true value") {
              CHECK(isClose(logmts, logmts_true));
            }

            // reset the input massfractions arrays
            portableFor(
                "Initialize mass fractions", 0, 1, PORTABLE_LAMBDA(int i) {
                  v_massfractions[0] = 0.2;
                  v_massfractions[1] = 0.0;
                  v_massfractions[2] = 0.1;
                  v_massfractions[3] = 0.3;
                  v_massfractions[4] = 0.4;
                });

#ifdef PORTABILITY_STRATEGY_KOKKOS
            // Create Kokkos views on device for the input arrays
            // Create host-side mirrors
            Kokkos::View<Real[num]> v_newmfs("newmassfractions");
            auto newmassfractions = Kokkos::create_mirror_view(v_newmfs);
            Kokkos::View<Real[mnum]> v_deltamfs("massfractiontransfers");
            auto h_deltamfs = Kokkos::create_mirror_view(v_deltamfs);
#else
            // Otherwise just create arrays to contain values and create pointers to
            // be passed to the functions in place of the Kokkos views
            std::array<Real, num> newmassfractions;
            auto v_newmfs = newmassfractions.data();
            std::array<Real, mnum> h_deltamfs;
            auto v_deltamfs = h_deltamfs.data();
#endif // PORTABILITY_STRATEGY_KOKKOS

            constexpr std::array<Real, num> newmassfractions_true{
                0., 0.800012617417208613, 0.1999873825827914, 0., 0.};
            constexpr std::array<Real, mnum> deltamfs_true{
                0.3, 0., 0., 0., 0.4, 0., 0., 0., 0.2, 0.100012617417208613};

            WHEN("A massfraction update with SmallStepMFUpdate is performed") {
              SmallStepMFUpdate(logmts_true, num, v_massfractions, v_order, v_logrates,
                                v_deltamfs, v_newmfs);
              for (int l = 0; l < num; l++) {
                std::cout << "Phase: " << v_order[l]
                          << "  initial mass fractions: " << v_massfractions[v_order[l]]
                          << std::endl;
                std::cout << "         "
                          << "   final mass fractions: " << v_newmfs[v_order[l]]
                          << std::endl;
              }
              for (int l = 0; l < mnum; l++) {
                std::cout << "From phase i to phase k, ik ( x means 0x): " << v_fromto[l]
                          << "   Mass fraction transfer: " << v_deltamfs[l] << std::endl;
              }
#ifdef PORTABILITY_STRATEGY_KOKKOS
              Kokkos::fence();
              Kokkos::deep_copy(newmassfractions, v_newmfs);
              Kokkos::deep_copy(h_deltamfs, v_deltamfs);
#endif // PORTABILITY_STRATEGY_KOKKOS
              THEN("The new mass fractions should be equal to the true values") {
                array_compare(num, phaseorder, massfractions, newmassfractions,
                              newmassfractions_true, "Phase", "old massfractions");
                array_compare(mnum, dgibbs, fromto, h_deltamfs, deltamfs_true,
                              "DeltaGibbs", "FromTo");
              }
            }
          }
        }
      }
    }
  }
}
//#endif // PORTABILITY_STRATEGY_NONE
#endif // PORTABILITY_STRATEGY_KOKKOS
