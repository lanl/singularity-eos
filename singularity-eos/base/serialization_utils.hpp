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

#ifndef SINGULARITY_EOS_BASE_SERIALIZATION_UTILS_
#define SINGULARITY_EOS_BASE_SERIALIZATION_UTILS_

#include <numeric>
#include <utility>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>

#include <singularity-eos/base/constants.hpp>

namespace singularity {

struct MemberResizer {
  template <typename Collection_t>
  void operator()(Collection_t &collection, std::size_t count) {
    collection.resize(count);
  }
};

template <typename Container_t, typename Resizer_t = MemberResizer>
struct BulkSerializer {
  Container_t eos_objects;

  BulkSerializer() = default;
  BulkSerializer(const Container_t &eos_objects_) : eos_objects(eos_objects_) {}
  void Finalize() {
    for (auto &eos : eos_objects) {
      eos.Finalize();
    }
  }
  const auto Size() const { return eos_objects.size(); }
  std::size_t SerializedSizeInBytes() const {
    return sizeof(std::size_t) + accumulate_([](std::size_t tot, const auto &eos) {
             return tot + eos.SerializedSizeInBytes();
           });
  }
  std::size_t SharedMemorySizeInBytes() const {
    return accumulate_([](std::size_t tot, const auto &eos) {
      return tot + eos.SharedMemorySizeInBytes();
    });
  }
  std::size_t Serialize(char *dst) {
    std::size_t vsize = eos_objects.size();
    memcpy(dst, &vsize, sizeof(std::size_t));
    std::size_t offst = sizeof(std::size_t);
    for (auto &eos : eos_objects) {
      offst += eos.Serialize(dst + offst);
    }
    PORTABLE_ALWAYS_REQUIRE(offst == SerializedSizeInBytes(), "Serialization failed!");
    return offst;
  }
  auto Serialize() {
    std::size_t size = SerializedSizeInBytes();
    char *dst = (char *)malloc(size);
    std::size_t new_size = Serialize(dst);
    PORTABLE_ALWAYS_REQUIRE(size == new_size, "Serialization failed!");
    return std::make_pair(size, dst);
  }
  std::size_t DeSerialize(char *src,
                          const SharedMemSettings &stngs = DEFAULT_SHMEM_STNGS) {
    std::size_t vsize;
    memcpy(&vsize, src, sizeof(std::size_t));
    Resizer_t()(eos_objects, vsize);
    SharedMemSettings stngs_loc = stngs; // non-const
    std::size_t offst = sizeof(std::size_t);
    std::size_t shared_offst = 0;
    for (auto &eos : eos_objects) {
      std::size_t shared_bytes = eos.SharedMemorySizeInBytes();
      offst += eos.DeSerialize(src + offst, stngs_loc);
      stngs_loc.data += shared_bytes;
      shared_offst += shared_bytes;
    }
    PORTABLE_ALWAYS_REQUIRE(offst == SerializedSizeInBytes(), "De-Serialization failed!");
    PORTABLE_ALWAYS_REQUIRE(shared_offst == SharedMemorySizeInBytes(),
                            "De-Serialization into shared memory failed!");
    return offst;
  }

 private:
  template <typename F>
  std::size_t accumulate_(const F &f) const {
    std::size_t init = 0;
    return std::accumulate(eos_objects.cbegin(), eos_objects.cend(), init, f);
  }
};

template <typename EOS>
using VectorSerializer = BulkSerializer<std::vector<EOS>, MemberResizer>;

} // namespace singularity
#endif // SINGULARITY_EOS_BASE_SERIALIZATION_UTILS_
