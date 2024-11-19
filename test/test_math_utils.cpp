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

#include <iostream>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/fast-math/logs.hpp>
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>
#include <singularity-eos/eos/eos.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <test/eos_unit_test_helpers.hpp>

SCENARIO("EOS Variant Type", "[Variant][EOS]") {
  using singularity::EOS;
  // print out the eos type
  std::cout << demangle(typeid(EOS).name()) << std::endl;
}

SCENARIO("Test that we can either throw an error on host or do nothing on device",
         "[RequireMaybe]") {
  // TODO(JMM): For whatever reason, the preprocessor does not like it if I
  // call `PORTABLE_ALWAYS_THROW_OR_ABORT
  REQUIRE_MAYBE_THROWS(PORTABLE_ALWAYS_THROW_OR_ABORT("Error message"));
}

#ifndef SINGULARITY_USE_TRUE_LOG_GRIDDING
#ifndef SINGULARITY_NQT_PORTABLE
SCENARIO("Test that the fast math magic numbers are all correct", "[FastMath]") {
  REQUIRE(singularity::FastMath::FP64LE::one_as_int == 4607182418800017408);
  REQUIRE(singularity::FastMath::FP64LE::scale_down == 2.22044604925031e-16);
  REQUIRE(singularity::FastMath::FP64LE::scale_up == 4503599627370496.0);
  REQUIRE(singularity::FastMath::FP64LE::mantissa_mask == 4503599627370495);
  REQUIRE(singularity::FastMath::FP64LE::low_mask == 67108863);
}
#endif
#endif

SCENARIO("Test that fast logs are invertible and run on device", "[FastMath]") {
  GIVEN("A set of values to invert over a large dynamic range") {
    constexpr Real LXMIN = -20;
    constexpr Real LXMAX = 32;
    constexpr int NX = 1000;
    constexpr Real DLX = (LXMAX - LXMIN) / (NX - 1);
    Real *x = (Real *)PORTABLE_MALLOC(NX * sizeof(Real));
    portableFor(
        "Set x values", 0, NX, PORTABLE_LAMBDA(const int i) {
          const Real lx = LXMIN + i * DLX;
          x[i] = std::pow(10., lx);
        });
    THEN("The fast exp of the fast log returns the original") {
      int nw_ie = 0;
#ifdef PORTABILITY_STRATEGY_KOKKOS
      using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
      Kokkos::View<int, atomic_view> n_wrong_ie("wrong_ie");
#else
      PortableMDArray<int> n_wrong_ie(&nw_ie, 1);
#endif
      portableFor(
          "try out the fast math", 0, NX, PORTABLE_LAMBDA(const int i) {
            constexpr Real machine_eps = std::numeric_limits<Real>::epsilon();
            constexpr Real acceptable_err = 1000 * machine_eps;
            const Real lx = singularity::FastMath::log10(x[i]);
            const Real elx = singularity::FastMath::pow10(lx);
            const Real rel_err = 2.0 * std::abs(x[i] - elx) /
                                 (std::abs(x[i]) + std::abs(elx) + machine_eps);
            n_wrong_ie() += (rel_err > acceptable_err);
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::deep_copy(nw_ie, n_wrong_ie);
#endif
      REQUIRE(nw_ie == 0);
    }
    PORTABLE_FREE(x);
  }
}

SCENARIO("Rudimentary test of the root finder", "[RootFinding1D]") {

  GIVEN("Root finding") {
    using namespace RootFinding1D;

    THEN("A root can be found for shift = 1, scale = 2, offset = 0.5") {
      int ntimes = 100;
      Real guess = 0;
      Real root;
      Status status;
      Real shift = 1;
      Real scale = 2;
      Real offset = 0.5;

      auto f = PORTABLE_LAMBDA(const Real x) { return myAtan(x, shift, scale, offset); };
      Status *statusesp = (Status *)PORTABLE_MALLOC(ntimes * sizeof(Status));
      Real *rootsp = (Real *)PORTABLE_MALLOC(ntimes * sizeof(Real));
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<Status *, Kokkos::MemoryTraits<Kokkos::Unmanaged>> statuses(statusesp,
                                                                               ntimes);
      Kokkos::View<Real *, Kokkos::MemoryTraits<Kokkos::Unmanaged>> roots(rootsp, ntimes);
#else
      PortableMDArray<Status> statuses;
      PortableMDArray<Real> roots;
      statuses.NewPortableMDArray(statusesp, ntimes);
      roots.NewPortableMDArray(rootsp, ntimes);
#endif
      portableFor(
          "find roots", 0, ntimes, PORTABLE_LAMBDA(const int i) {
            statuses(i) = regula_falsi(f, 0, guess, -1, 3, 1e-10, 1e-10, roots(i));
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::View<Status> s_copy(statuses, 0);
      Kokkos::View<Real> r_copy(roots, 0);
      Kokkos::deep_copy(root, r_copy);
      Kokkos::deep_copy(status, s_copy);
#else
      status = statuses(ntimes - 1);
      root = roots(ntimes - 1);
#endif

      PORTABLE_FREE(statusesp);
      PORTABLE_FREE(rootsp);
      REQUIRE(status == Status::SUCCESS);
      REQUIRE(isClose(root, 0.744658));
    }
  }
}
