//------------------------------------------------------------------------------
// © 2021-2024. Triad National Security, LLC. All rights reserved.  This
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
#include <iostream> // debug
#include <limits>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/eos/eos.hpp>

#ifdef SINGULARITY_BUILD_CLOSURE
#include <singularity-eos/eos/singularity_eos.hpp>
#endif

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <test/eos_unit_test_helpers.hpp>

#ifdef SPINER_USE_HDF
#ifdef SINGULARITY_TEST_STELLAR_COLLAPSE
SCENARIO("Test 3D reinterpolation to fast log grid", "[StellarCollapse]") {
  WHEN("We generate a 3D DataBox of a 2D power law times a line") {
    using singularity::StellarCollapse;
    constexpr int N2 = 100;
    constexpr int N1 = 101;
    constexpr int N0 = 102;
    StellarCollapse::Grid_t g2(1.0 / N2, 1, N2);
    StellarCollapse::Grid_t g1(1, 3, N1);
    StellarCollapse::Grid_t g0(2, 4, N0);
    StellarCollapse::DataBox db(N2, N1, N0);

    db.setRange(2, g2);
    db.setRange(1, g1);
    db.setRange(0, g0);

    for (int i2 = 0; i2 < N2; ++i2) {
      Real x2 = g2.x(i2);
      for (int i1 = 0; i1 < N1; ++i1) {
        Real lx1 = g1.x(i1);
        for (int i0 = 0; i0 < N0; ++i0) {
          Real lx0 = g0.x(i0);
          db(i2, i1, i0) = std::log10(x2) + 2 * lx1 + 2 * lx0;
        }
      }
    }

    THEN("The databox should interpolate values correctly") {
      const Real x2 = 0.5;
      const Real lx1 = 2;
      const Real x1 = std::pow(10., lx1);
      const Real lx0 = 3;
      const Real x0 = std::pow(10., lx0);
      REQUIRE(isClose(db.interpToReal(x2, lx1, lx0), std::log10(x2 * x1 * x1 * x0 * x0),
                      1e-5));
    }

    THEN("We can re-interpolate to fast log space") {
      StellarCollapse::DataBox scratch(N2, N1, N0);
      StellarCollapse::dataBoxToFastLogs(db, scratch, true);

      AND_THEN("The fast-log gridded table contains correct ranges") {
        REQUIRE(db.range(2) == g2);
        REQUIRE(db.range(1).nPoints() == N1);
        REQUIRE(isClose(db.range(1).min(),
                        singularity::FastMath::log10(std::pow(10, g1.min())), 1e-12));
        REQUIRE(isClose(db.range(1).max(),
                        singularity::FastMath::log10(std::pow(10, g1.max())), 1e-12));
        REQUIRE(db.range(0).nPoints() == N0);
        REQUIRE(isClose(db.range(0).min(),
                        singularity::FastMath::log10(std::pow(10, g0.min())), 1e-12));
        REQUIRE(isClose(db.range(0).max(),
                        singularity::FastMath::log10(std::pow(10, g0.max())), 1e-12));
      }

      AND_THEN("The re-interpolated fast log is a sane number") {}

      AND_THEN("The fast-log table approximately interpolates the power law") {
        const Real x2 = 0.5;
        const Real x1 = 100;
        const Real x0 = 1000;
        const Real lx1 = singularity::FastMath::log10(x1);
        const Real lx0 = singularity::FastMath::log10(x0);
        const Real lval_interp = db.interpToReal(x2, lx1, lx0);
        const Real val_interp = singularity::FastMath::pow10(lval_interp);
        const Real val_true = x2 * x1 * x1 * x0 * x0;
        const Real rel_diff =
            0.5 * std::abs(val_interp - val_true) / (val_true + val_interp);
        REQUIRE(rel_diff <= 1e-3);
      }
      scratch.finalize();
    }
    db.finalize();
  }
}

SCENARIO("Stellar Collapse EOS", "[StellarCollapse]") {
  using singularity::IdealGas;
  using singularity::StellarCollapse;
  const std::string savename = "stellar_collapse_ideal_2.sp5";
  GIVEN("A stellar collapse EOS") {
    const std::string filename = "./goldfiles/stellar_collapse_ideal.h5";
    THEN("We can load the file") { // don't bother filtering bmod here.
      StellarCollapse sc(filename, false, false);
      AND_THEN("Some properties we expect for ideal gas hold") {
        std::vector<Real> lambda(2);
        Real rho, t, sie, p, cv, b, dpde, dvdt;
        sc.ValuesAtReferenceState(rho, t, sie, p, cv, b, dpde, dvdt, lambda);
        Real yemin = sc.YeMin();
        Real yemax = sc.YeMax();
        int N = 123;
        Real dY = (yemax - yemin) / (N + 1);
        for (int i = 0; i < N; ++i) {
          Real Ye = yemin + i * dY;
          lambda[0] = Ye;
          REQUIRE(isClose(sie, sc.InternalEnergyFromDensityTemperature(rho, t, lambda)));
          Real Xa, Xh, Xn, Xp, Abar, Zbar;
          sc.MassFractionsFromDensityTemperature(rho, t, Xa, Xh, Xn, Xp, Abar, Zbar,
                                                 lambda);
          REQUIRE(isClose(Ye, Xp));
        }
        Real rhomin = sc.rhoMin();
        Real rhomax = sc.rhoMax();
        Real drho = (rhomax - rhomin) / (N + 1);
        for (int i = 0; i < N; ++i) {
          Real rho = rhomin + i * drho;
          REQUIRE(isClose(sie, sc.InternalEnergyFromDensityTemperature(rho, t, lambda)));
        }
      }
      GIVEN("An Ideal Gas equation of state") {
        constexpr Real gamma = 1.4;
        constexpr Real mp = 1.67262171e-24;
        constexpr Real kb = 1.3806505e-16;
        constexpr Real Cv = kb / (mp * (gamma - 1)); // mean molecular weight = mp
        IdealGas ig(gamma - 1, Cv);
        auto ig_d = ig.GetOnDevice();
        THEN("The tabulated gamma Stellar Collapse and the gamma agree roughly") {
          Real yemin = sc.YeMin();
          Real yemax = sc.YeMax();
          Real tmin = sc.TMin();
          Real tmax = sc.TMax();
          Real ltmin = std::log10(tmin);
          Real ltmax = std::log10(tmax);
          Real lrhomin = std::log10(sc.rhoMin());
          Real lrhomax = std::log10(sc.rhoMax());
          auto sc_d = sc.GetOnDevice();

          int nwrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
          using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
          Kokkos::View<int, atomic_view> nwrong_d("wrong");
#else
          PortableMDArray<int> nwrong_d(&nwrong_h, 1);
#endif

          const int N = 123;
          const Real dY = (yemax - yemin) / (N + 1);
          const Real dlT = (ltmax - ltmin) / (N + 1);
          const Real dlR = (lrhomax - lrhomin) / (N + 1);
          portableFor(
              "fill eos", 0, N, 0, N, 0, N,
              PORTABLE_LAMBDA(const int &k, const int &j, const int &i) {
                Real lambda[2];
                Real Ye = yemin + k * dY;
                Real lT = ltmin + j * dlT;
                Real lR = lrhomin + i * dlR;
                Real T = std::pow(10., lT);
                Real R = std::pow(10., lR);
                Real e1, e2, p1, p2, cv1, cv2, b1, b2;
                unsigned long output = (singularity::thermalqs::pressure |
                                        singularity::thermalqs::specific_internal_energy |
                                        singularity::thermalqs::specific_heat |
                                        singularity::thermalqs::bulk_modulus);
                lambda[0] = Ye;

                sc_d.FillEos(R, T, e1, p1, cv1, b1, output, lambda);
                ig_d.FillEos(R, T, e2, p2, cv2, b2, output, lambda);
                if (!isClose(e1, e2)) {
                  nwrong_d() += 1;
                }
                if (!isClose(p1, p2)) {
                  nwrong_d() += 1;
                }
                if (!isClose(cv1, cv2)) {
                  nwrong_d() += 1;
                }
                if (!isClose(b1, b2)) {
                  nwrong_d() += 1;
                }
              });
#ifdef PORTABILITY_STRATEGY_KOKKOS
          Kokkos::deep_copy(nwrong_h, nwrong_d);
#endif
          REQUIRE(nwrong_h == 0);
          sc_d.Finalize();
        }
        ig_d.Finalize();
        ig.Finalize();
      }
      AND_THEN("We can save the file to SP5") {
        sc.Save(savename);
        AND_THEN("We can load the sp5 file") {
          StellarCollapse sc2(savename, true);
          AND_THEN("The two stellar collapse EOS's agree") {

            Real yemin = sc.YeMin();
            Real yemax = sc.YeMax();
            Real tmin = sc.TMin();
            Real tmax = sc.TMax();
            Real ltmin = std::log10(tmin);
            Real ltmax = std::log10(tmax);
            Real lrhomin = std::log10(sc.rhoMin());
            Real lrhomax = std::log10(sc.rhoMax());
            REQUIRE(yemin == sc2.YeMin());
            REQUIRE(yemax == sc2.YeMax());
            REQUIRE(sc.TMin() == sc2.TMin());
            REQUIRE(sc.TMax() == sc2.TMax());
            REQUIRE(isClose(lrhomin, std::log10(sc2.rhoMin()), 1e-12));
            REQUIRE(isClose(lrhomax, std::log10(sc2.rhoMax()), 1e-12));

            auto sc1_d = sc.GetOnDevice();
            auto sc2_d = sc2.GetOnDevice();

            int nwrong_h = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
            using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
            Kokkos::View<int, atomic_view> nwrong_d("wrong");
#else
            PortableMDArray<int> nwrong_d(&nwrong_h, 1);
#endif

            const int N = 123;
            constexpr Real gamma = 1.4;
            const Real dY = (yemax - yemin) / (N + 1);
            const Real dlT = (ltmax - ltmin) / (N + 1);
            const Real dlR = (lrhomax - lrhomin) / (N + 1);
            portableFor(
                "fill eos", 0, N, 0, N, 0, N,
                PORTABLE_LAMBDA(const int &k, const int &j, const int &i) {
                  Real lambda[2];
                  Real Ye = yemin + k * dY;
                  Real lT = ltmin + j * dlT;
                  Real lR = lrhomin + i * dlR;
                  Real T = std::pow(10., lT);
                  Real R = std::pow(10., lR);
                  Real e1, e2, p1, p2, cv1, cv2, b1, b2, s1, s2;
                  unsigned long output =
                      (singularity::thermalqs::pressure |
                       singularity::thermalqs::specific_internal_energy |
                       singularity::thermalqs::specific_heat |
                       singularity::thermalqs::bulk_modulus);
                  lambda[0] = Ye;

                  sc1_d.FillEos(R, T, e1, p1, cv1, b1, output, lambda);
                  sc2_d.FillEos(R, T, e2, p2, cv2, b2, output, lambda);
                  // Fill entropy. Will need to change later.
                  s1 = sc1_d.EntropyFromDensityTemperature(R, T, lambda);
                  s2 = p2 * std::pow(R, -gamma); // ideal
                  if (!isClose(e1, e2)) nwrong_d() += 1;
                  if (!isClose(p1, p2)) nwrong_d() += 1;
                  if (!isClose(cv1, cv2)) nwrong_d() += 1;
                  if (!isClose(b1, b2)) nwrong_d() += 1;
                  if (!isClose(s1, s2)) nwrong_d() += 1;
                });
#ifdef PORTABILITY_STRATEGY_KOKKOS
            Kokkos::deep_copy(nwrong_h, nwrong_d);
#endif
            REQUIRE(nwrong_h == 0);

            sc1_d.Finalize();
            sc2_d.Finalize();
          }
          sc2.Finalize();
        }
      }
      sc.Finalize();
    }
  }
}
#endif // SINGULARITY_TEST_STELLAR_COLLAPSE
#endif // USE_HDF5
