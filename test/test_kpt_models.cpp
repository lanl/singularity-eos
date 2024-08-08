//------------------------------------------------------------------------------
// © 2024. Triad National Security, LLC. All rights reserved.  This
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
#include <singularity-eos/closure/kinetic_phasetransition_methods.hpp>
#include <singularity-eos/closure/kinetic_phasetransition_models.hpp>
#include <test/eos_unit_test_helpers.hpp>

using singularity::LogRatesCGModel;
using singularity::SortGibbs;

SCENARIO("First log rate test") {
  GIVEN("Parameters for a 5 phase system") {
    constexpr Real w[25] = {0.,  0.8510E+01, 20., 20., 20., 0.8510E+01, 0.,  20., 20.,
                            20., 20.,        20., 0.,  20., 20.,        20., 20., 20.,
                            0.,  20.,        20., 20., 20., 20.,        0.};

    constexpr Real b[25] = {0.,         0.3060E-04, 0.1000E-05, 0.1000E-05, 0.1000E-05,
                            0.3060E-04, 0.,         0.1000E-05, 0.1000E-05, 0.1000E-05,
                            0.1000E-05, 0.1000E-05, 0.,         0.1000E-05, 0.1000E-05,
                            0.1000E-05, 0.1000E-05, 0.1000E-05, 0.,         0.1000E-05,
                            0.1000E-05, 0.1000E-05, 0.1000E-05, 0.1000E-05, 0.};

    GIVEN("Gibbs free energy and order") {
      constexpr int num = 5;
      constexpr int mnum = num * (num - 1) / 2;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create Kokkos views on device for the input arrays
      Kokkos::View<Real[num]> v_gibbs("gibbs");
      // order phases from high to low gibbs
      Kokkos::View<int[num]> v_order("order");
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
          "Initialize gibbs and order", 0, 1, PORTABLE_LAMBDA(int i) {
            v_gibbs[0] = 0.0292971;
            v_gibbs[1] = 0.0286937;
            v_gibbs[2] = 0.0287409;
            v_gibbs[3] = 1.0e99;
            v_gibbs[4] = 0.02940076;
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      // Create host-side mirrors of the inputs and copy the inputs. These are
      // just used for the comparisons
      auto gibbs = Kokkos::create_mirror_view(v_gibbs);
      auto order = Kokkos::create_mirror_view(v_order);
      Kokkos::deep_copy(gibbs, v_gibbs);
      Kokkos::deep_copy(order, v_order);
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
      }
      // create deltagibbs and fromphasetophase
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<Real[mnum]> v_dgibbs("dgibbs");
      Kokkos::View<int[mnum]> v_fromto("fromto");
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
      auto dgibbs = Kokkos::create_mirror_view(v_dgibbs);
      auto fromto = Kokkos::create_mirror_view(v_fromto);
      Kokkos::deep_copy(dgibbs, v_dgibbs);
      Kokkos::deep_copy(fromto, v_fromto);
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
          std::cout << "From phase i to phase k, ik: " << v_fromto[l]
                    << "   LogRik: " << v_logrates[l] << std::endl;
        }
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        Kokkos::deep_copy(h_logrates, v_logrates);
#endif // PORTABILITY_STRATEGY_KOKKOS
        THEN("The returned LogRij(w,b,gibbs,order) should be equal to the true value") {
          array_compare(mnum, dgibbs, fromto, h_logrates, logrates_true, "DeltaGibbs",
                        "FromTo");
        }
      }
    }
  }
}
