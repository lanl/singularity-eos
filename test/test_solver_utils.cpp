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

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <test/eos_unit_test_helpers.hpp>

SCENARIO("Kahan summation", "[Kahan]") {
  GIVEN("An array that might be subject to catastrophic cancellation") {
    constexpr std::size_t N = 100;
    constexpr Real delta = 1e-12;
    Real *A = (Real *)PORTABLE_MALLOC(N * sizeof(Real));
    portableFor(
        "Set values of array", 0, N, PORTABLE_LAMBDA(const int i) { A[i] = 1 - delta; });
    Real sum = 0;
    WHEN("We use Kahan summation to compute the sum") {
      portableReduce(
          "compute Kahan sum", 0, 1,
          PORTABLE_LAMBDA(const int, Real &s) {
            s = singularity::mix_impl::sum_neumaier(A, N-1, 1); // with offset
          },
          sum);
      THEN("The total is what we expect within machine epsilon") {
        constexpr Real strue = (N - 1) * (1 - delta);
        REQUIRE(isClose(sum, strue, singularity::robust::EPS()));
      }
    }
    PORTABLE_FREE(A);
  }
}

SCENARIO("2x2 matrix solve", "[Matrix][2x2]") {
  GIVEN("A mildly ill-conditioned size 2 linear system") {
    constexpr std::size_t N = 2;
    constexpr Real delta = 1e-11;
    constexpr Real B0_true = 2 * (1 + delta) / delta - 1 / delta;
    constexpr Real B1_true = -1 / delta;
    Real *A = (Real *)PORTABLE_MALLOC(N * N * sizeof(Real));
    Real *B = (Real *)PORTABLE_MALLOC(N * sizeof(Real));
    Real *scr = (Real *)PORTABLE_MALLOC(N * sizeof(Real));
    portableFor(
        "Set values of A and B", 0, 1, PORTABLE_LAMBDA(const int) {
          A[0] = 1;
          A[1] = 1 + delta;
          A[2] = 1;
          A[3] = 1;
          B[0] = 1;
          B[1] = 2;
        });

    WHEN("We solve it with solve_Ax_b_wscr") {
      std::size_t nfail = 0;
      portableReduce(
          "Solve Ax = b", 0, 1,
          PORTABLE_LAMBDA(const int, std::size_t &fails) {
            singularity::mix_impl::solve_Ax_b_wscr(N, A, B, scr);
            fails += !isClose(B[0], B0_true, std::sqrt(delta));
            fails += !isClose(B[1], B1_true, std::sqrt(delta));
          },
          nfail);
      THEN("The error is small") { REQUIRE(nfail == 0); }
    }

    PORTABLE_FREE(A);
    PORTABLE_FREE(B);
    PORTABLE_FREE(scr);
  }
}

SCENARIO("3x3 matrix solve", "[Matrix][3x3]") {
  GIVEN("A mildly ill-conditioned size 3 linear system") {
    constexpr std::size_t N = 3;
    constexpr Real delta = 1e-10;
    constexpr Real denom = delta * (3 + delta);
    constexpr Real B0_true = (delta - 1) / denom;
    constexpr Real B1_true = (delta - 1) / denom;
    constexpr Real B2_true = 2 * (1 + delta) / denom;
    Real *A = (Real *)PORTABLE_MALLOC(N * N * sizeof(Real));
    Real *B = (Real *)PORTABLE_MALLOC(N * sizeof(Real));
    Real *scr = (Real *)PORTABLE_MALLOC(N * sizeof(Real));
    portableFor(
        "Set values of A and B", 0, 1, PORTABLE_LAMBDA(const int) {
          for (std::size_t row = 0; row < N; ++row) {
            for (std::size_t column = 0; column < N; ++column) {
              A[row * N + column] = 1 + (row == column) * delta;
            }
          }
          B[0] = 1;
          B[1] = 1;
          B[2] = 2;
        });

    WHEN("We solve it with solve_Ax_b_wscr") {
      std::size_t nfail = 0;
      portableReduce(
          "Solve Ax = b", 0, 1,
          PORTABLE_LAMBDA(const int, std::size_t &fails) {
            singularity::mix_impl::solve_Ax_b_wscr(N, A, B, scr);
            fails += !isClose(B[0], B0_true, std::sqrt(delta));
            fails += !isClose(B[1], B1_true, std::sqrt(delta));
            fails += !isClose(B[2], B2_true, std::sqrt(delta));
          },
          nfail);
      THEN("The error is small") { REQUIRE(nfail == 0); }
    }

    PORTABLE_FREE(A);
    PORTABLE_FREE(B);
    PORTABLE_FREE(scr);
  }
}
