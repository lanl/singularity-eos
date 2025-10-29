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

#include <array>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/generic_indexer.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif
#include <catch2/matchers/catch_matchers_floating_point.hpp>

SCENARIO("Generic Indexer", "[GenericIndexer]") {
  GIVEN("An array of arays populated with data") {
    constexpr std::size_t n = 5;
    constexpr std::size_t m = 7;
    auto data = std::array<std::array<Real, m>, n>{};
    for (std::size_t i = 0; i < n; i++) {
      for (std::size_t j = 0; j < m; j++) {
        data[i][j] = (i + 1) * (j + 1);
      }
    }
    WHEN("A flattened array is created with the same size and populated with "
         "the same data") {
      auto flat_data = std::array<Real, n * m>{};
      for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = 0; j < m; j++) {
          flat_data[i * m + j] = data[i][j];
        }
      }
      WHEN("We use the generic indexer to index into material data") {
        struct constant_offset {
          std::size_t offset;
          constexpr constant_offset(std::size_t offset_) : offset{offset_} {}
          constexpr std::size_t operator[](std::size_t i) {
            return i * offset;
          }
        };
        auto mat_indexer = singularity::GenericIndexer(flat_data.data(), constant_offset{m});
        THEN("The data returned should be same as if we accessed the array of "
             "arrays") {
          for (std::size_t i = 0; i < n; i++) {
            auto* mat_array = &mat_indexer[i];
            for (std::size_t j = 0; j < m; j++) {
              INFO("i: " << i << " j: " << j);
              CHECK_THAT(data[i][j], Catch::Matchers::WithinRel(mat_array[j], 1.0e-12));
            }
          }
        }
      }
    }
  }
}

