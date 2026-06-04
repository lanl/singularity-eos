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

#ifndef SINGULARITY_EOS_BASE_TABLE_UTILS_HPP_
#define SINGULARITY_EOS_BASE_TABLE_UTILS_HPP_
#ifdef SINGULARITY_USE_SPINER

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/constants.hpp>

#include <singularity-utils/bounds.hpp> // Bounds
#include <spiner/databox.hpp>

namespace singularity {
namespace table_utils {

// JMM: Making this a struct with static methods, rather than a
// namespace, saves a few "friend" declarations. Broadly these methods
// rely on the EOS providing a `GetDataBoxPointers_` method and
// declaring this class a friend. Doing so lets us loop over each data
// box member of a spiner-based EOS ond operate on it. For example, by
// calling `finalize` or `getOnDevice`. The intent here is that the
// static methods of `SpinerTricks` are called from class methods of
// spiner-based EOS's.
template <typename EOS>
class SpinerTricks { // class instead of struct for -Wmismatched-tags
 public:
  static auto GetOnDevice(EOS *peos_h) {
    // trivially copy all but dynamic memory
    EOS eos_d = *peos_h;
    auto pdbs_d = eos_d.GetDataBoxPointers_();
    auto pdbs_h = peos_h->GetDataBoxPointers_();
    int idb = 0;
    for (auto *pdb_d : pdbs_d) {
      auto *pdb_h = pdbs_h[idb++];
      *pdb_d = pdb_h->getOnDevice();
    }
    // set memory status
    eos_d.memoryStatus_ = DataStatus::OnDevice;
    return eos_d;
  }
  static void Finalize(EOS *peos) {
    if (peos->memoryStatus_ != DataStatus::UnManaged) {
      for (auto *pdb : peos->GetDataBoxPointers_()) {
        // JMM: on host, some databoxes might be slices.
        if (pdb->ownsAllocatedMemory()) {
          pdb->finalize();
        }
      }
    }
    peos->memoryStatus_ = DataStatus::Deallocated;
  }
  static std::size_t DynamicMemorySizeInBytes(const EOS *peos) {
    std::size_t out = 0;
    for (const auto *pdb : peos->GetDataBoxPointers_()) {
      out += pdb->sizeBytes();
    }
    return out;
  }
  static std::size_t DumpDynamicMemory(char *dst, const EOS *peos) {
    std::size_t offst = 0;
    for (const auto *pdb : peos->GetDataBoxPointers_()) {
      std::size_t size = pdb->sizeBytes();
      memcpy(dst + offst, pdb->data(), size);
      offst += size;
    }
    return offst;
  }
  static std::size_t SetDynamicMemory(char *src, EOS *peos) {
    std::size_t offst = 0;
    for (auto *pdb : peos->GetDataBoxPointers_()) {
      offst += pdb->setPointer(src + offst);
    }
    peos->memoryStatus_ = DataStatus::UnManaged;
    return offst;
  }
  static bool DataBoxesPointToSameMemory(const EOS &eos_a, const EOS &eos_b) {
    const auto pdbs_a = eos_a.GetDataBoxPointers_();
    const auto pdbs_b = eos_b.GetDataBoxPointers_();
    int idb = 0;
    for (const auto *pdb_a : pdbs_a) {
      const auto &db_a = *pdb_a;
      const auto &db_b = *(pdbs_b[idb++]);
      if (&(db_a(0)) != &(db_b(0))) return false;
    }
    return true;
  }
  static bool DataBoxesPointToDifferentMemory(const EOS &eos_a, const EOS &eos_b) {
    return !DataBoxesPointToSameMemory(eos_a, eos_b);
  }
};
} // namespace table_utils
} // namespace singularity

#endif // SINGULARITY_USE_SPINER
#endif // SINGULARITY_EOS_BASE_TABLE_UTILS_HPP_
