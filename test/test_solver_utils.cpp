//------------------------------------------------------------------------------
// © 2021-2026. Triad National Security, LLC. All rights reserved.  This
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
            s = singularity::mix_impl::sum_neumaier(A, N - 1, 1); // with offset
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
    // Kokkos QR requires 2 scratch vectors, instead of the 1 required
    // by Cramer's rule or Eigen
    Real *scr = (Real *)PORTABLE_MALLOC(2 * N * sizeof(Real));
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

SCENARIO("We can enforce conservation on a state that might emerge out of a PTE solver"
         "[EnforceMassVolumesSum][EnforceEnergiesSum]") {
  constexpr std::size_t nmat = 3;
  constexpr Real tot_rho = 2.5;
  constexpr Real tot_u = 1.23e12;
  constexpr Real tot_sie = tot_u / tot_rho;

  Real *mu = (Real *)PORTABLE_MALLOC(nmat * sizeof(Real));
  Real *rhobar = (Real *)PORTABLE_MALLOC(nmat * sizeof(Real));
  Real *rho = (Real *)PORTABLE_MALLOC(nmat * sizeof(Real));
  Real *vfrac = (Real *)PORTABLE_MALLOC(nmat * sizeof(Real));
  Real *sie = (Real *)PORTABLE_MALLOC(nmat * sizeof(Real));

  GIVEN("A state that isn't quite right") {
    portableFor(
        "Set thermodyamic state", 0, 1, PORTABLE_LAMBDA(const int) {
          mu[0] = 0.025;
          mu[1] = 0.9; // this one's a metal
          mu[2] = 0.075;

          vfrac[0] = 0.55;
          vfrac[1] = 0.05;
          vfrac[2] = 0.4;
          for (int m = 0; m < nmat; ++m) {
            rhobar[m] = mu[m] * tot_rho;
            rho[m] = rhobar[m] / vfrac[m];
          }

          // perturb and reset
          rho[0] += -0.01;
          rho[1] += 0.004;
          rho[2] -= 0.002;
          for (int m = 0; m < nmat; ++m) {
            vfrac[m] = singularity::robust::ratio(rhobar[m], rho[m]);
          }

          // already perturbed
          sie[0] = 0.5 * tot_u / rhobar[0] + 0.01 * tot_sie;
          sie[1] = 0 - 0.0025 * tot_sie;
          sie[2] = 0.5 * tot_u / rhobar[2] + -0.06 * tot_sie;
        });

    THEN("Things don't quite add up") {
      int nwrong = 0;
      portableReduce(
          "Check it's currently wrong", 0, 1,
          PORTABLE_LAMBDA(const int, int &nw) {
            Real test_rho = 0;
            Real test_u = 0;
            Real test_vfrac = 0;
            for (int m = 0; m < nmat; ++m) {
              test_rho += vfrac[m] * rho[m];
              test_u += vfrac[m] * rho[m] * sie[m];
              test_vfrac += vfrac[m];
            }
            if (isClose(test_vfrac, 1, 1e-12)) {
              nw += 1;
            }
            if (isClose(test_u, tot_u, 1e-12)) {
              nw += 1;
            }
            if (!isClose(test_rho, tot_rho, 1e-12)) {
              nw += 1;
            }
          },
          nwrong);
      REQUIRE(nwrong == 0);

      AND_WHEN("We enforce mass and volume fractions sum") {
        portableFor(
            "Enforce mass and volume fractions sum", 0, 1, PORTABLE_LAMBDA(const int) {
              singularity::MixUtils::EnforceMassVolumesSum(nmat, 1.0, rho, vfrac);
            });
        THEN("They do") {
          int nwrong = 0;
          portableReduce(
              "Check they sum right now", 0, 1,
              PORTABLE_LAMBDA(const int, int &nw) {
                Real test_rho = 0;
                Real test_vfrac = 0;
                for (int m = 0; m < nmat; ++m) {
                  test_rho += vfrac[m] * rho[m];
                  test_vfrac += vfrac[m];
                }
                if (!isClose(test_vfrac, 1, 1e-12)) {
                  nw += 1;
                }
                if (!isClose(test_rho, tot_rho, 1e-12)) {
                  nw += 1;
                }
              },
              nwrong);
          REQUIRE(nwrong == 0);
        }
        AND_WHEN("We enforce energies sum") {
          portableFor(
              "Enforce energies sum", 0, 1, PORTABLE_LAMBDA(const int) {
                singularity::MixUtils::EnforceEnergiesSum(nmat, tot_rho, tot_sie, rho,
                                                          vfrac, sie);
              });
          THEN("They do") {
            int nwrong = 0;
            portableReduce(
                "Check energies sum", 0, 1,
                PORTABLE_LAMBDA(const int, int &nw) {
                  Real test_u = 0;
                  for (int m = 0; m < nmat; ++m) {
                    test_u += vfrac[m] * rho[m] * sie[m];
                  }
                  if (!isClose(test_u, tot_u, 1e-12)) {
                    nw += 1;
                  }
                },
                nwrong);
            REQUIRE(nwrong == 0);
          }
        }
      }
    }
  }
  PORTABLE_FREE(mu);
  PORTABLE_FREE(rhobar);
  PORTABLE_FREE(rho);
  PORTABLE_FREE(vfrac);
  PORTABLE_FREE(sie);
}
