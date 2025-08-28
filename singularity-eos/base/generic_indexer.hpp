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

#include <utility>

namespace singularity {

// Index into an arbitrary array with an arbitrary map
template<typename arrT, typename mapT>
struct GenericIndexer {
  arrT arr_;
  mapT map_;

  template<typename arrT_, typename mapT_>
  constexpr
  GenericIndexer(arrT_&& arr_in, mapT_&& map_in)
      : arr_{std::forward<arrT_>(arr_in)}
      , map_{std::forward<mapT_>(map_in)}
  {}

  template<typename idxT>
  constexpr
  auto &operator[](const idxT i) { return arr_[map_[i]]; }

  template<typename idxT>
  constexpr
  const auto &operator[](const idxT i) const { return arr_[map_[i]]; }
};
// CTAD
template<typename arrT_, typename mapT_>
GenericIndexer(arrT_&& arr_in, mapT_&& map_in) -> GenericIndexer<arrT_, mapT_>;

} // namespace singularity
