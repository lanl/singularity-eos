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

#ifdef SINGULARITY_USE_SPINER

#include <cstdio>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>

#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/spiner_table_utils.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

using namespace singularity;
class DummyEOS;
using STricks = table_utils::SpinerTricks<DummyEOS>;

class DummyEOS {
 public:
  friend class table_utils::SpinerTricks<DummyEOS>;
  using Grid_t = Spiner::RegularGrid1D<std::size_t>;
  using DataBox = Spiner::DataBox<std::size_t, Grid_t>;

  DummyEOS() = default;
  DummyEOS(int N) : N_(N), a_(N), b_(N), memoryStatus_(DataStatus::OnHost) {
    for (int i = 0; i < N; ++i) {
      a_(i) = i;
      b_(i) = i + N;
    }
  }
  void Finalize() { STricks::Finalize(this); }
  DummyEOS GetOnDevice() { return STricks::GetOnDevice(this); }
  std::size_t DynamicMemorySizeInBytes() const {
    return STricks::DynamicMemorySizeInBytes(this);
  }
  std::size_t DumpDynamicMemory(char *dst) const {
    return STricks::DumpDynamicMemory(dst, this);
  }
  std::size_t SetDynamicMemory(char *src) { return STricks::SetDynamicMemory(src, this); }

  // Public so we can inspect in unit tests
  int N_ = -1;
  DataBox a_, b_;
  DataStatus memoryStatus_ = DataStatus::Deallocated;

  // Private so we can ensure class friendship is working
 private:
#define DBLIST &a_, &b_
  std::vector<const DataBox *> GetDataBoxPointers_() const {
    return std::vector<const DataBox *>{DBLIST};
  }
  std::vector<DataBox *> GetDataBoxPointers_() { return std::vector<DataBox *>{DBLIST}; }
#undef DBLIST
};

constexpr int N = 5;
SCENARIO("The SpinerTricks tool can copy databoxes to device and finalize them",
         "[SpinerTricks]") {
  WHEN("We initialize a host-side DummyEOS") {
    DummyEOS eos_h(N);
    REQUIRE(eos_h.memoryStatus_ == DataStatus::OnHost);
    THEN("We can copy it to device") {
      DummyEOS eos_d = eos_h.GetOnDevice();
      REQUIRE(eos_d.memoryStatus_ == DataStatus::OnDevice);
      AND_THEN("The device-side memory is correct") {
        int nwrong = 0;
        portableReduce(
            "Check Dummy EOS", 0, N,
            PORTABLE_LAMBDA(const int i, int &nw) {
              nw += (eos_d.a_(i) != i) + (eos_d.b_(i) != i + N);
            },
            nwrong);
        REQUIRE(nwrong == 0);
      }
      eos_d.Finalize();
      REQUIRE(eos_d.memoryStatus_ == DataStatus::Deallocated);
    }

    eos_h.Finalize();
    REQUIRE(eos_h.memoryStatus_ == DataStatus::Deallocated);
  }
}

SCENARIO("The SpinerTricks tool can serialize dynamic memory", "[SpinerTricks]") {
  GIVEN("A DummyEOS") {
    DummyEOS eos1(N);
    std::size_t size = eos1.DynamicMemorySizeInBytes();
    REQUIRE(eos1.memoryStatus_ == DataStatus::OnHost);
    REQUIRE(size == 2 * N * sizeof(std::size_t));
    THEN("We can serialize it") {
      char *data = (char *)malloc(size);
      eos1.DumpDynamicMemory(data);
      WHEN("We can de-serialize it twice") {
        DummyEOS eos2(N), eos3(N);
        eos2.SetDynamicMemory(data);
        eos3.SetDynamicMemory(data);
        THEN("The data is different from the original object") {
          REQUIRE(STricks::DataBoxesPointToDifferentMemory(eos1, eos2));
          REQUIRE(STricks::DataBoxesPointToDifferentMemory(eos1, eos3));
          AND_THEN("The serialized objects point to the same underlying memory as each "
                   "other") {
            REQUIRE(STricks::DataBoxesPointToSameMemory(eos2, eos3));
          }
        }
      }
      free(data);
    }
    eos1.Finalize();
  }
}

#endif // SINGULARITY_USE_SPINER
