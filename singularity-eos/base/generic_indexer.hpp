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

#ifndef SINGULARITY_EOS_BASE_GENERIC_INDEXER_HPP_
#define SINGULARITY_EOS_BASE_GENERIC_INDEXER_HPP_

#include <utility>

namespace singularity {

template <typename Arr, typename Map>
struct GenericIndexer {
  Arr arr_;
  Map map_;

  template <typename ArrT_, typename MapT_>
  constexpr GenericIndexer(ArrT_ &&arr_in, MapT_ &&map_in)
      : arr_(std::forward<ArrT_>(arr_in)), map_(std::forward<MapT_>(map_in)) {}

  // & : non-const lvalue
  template <class I>
  constexpr decltype(auto) operator[](I i) & {
    return arr_[map_[i]];
  }

  // const& : const lvalue
  template <class I>
  constexpr decltype(auto) operator[](I i) const & {
    return arr_[map_[i]];
  }

  // && : rvalue (indexer is a temporary) — forward arr_’s value category
  template <class I>
  constexpr decltype(auto) operator[](I i) && {
    return std::forward<Arr>(arr_)[map_[i]]; // move only if Arr is a value type
  }

  // const rvalue indexer
  template <class I>
  constexpr decltype(auto) operator[](I i) const && {
    return std::forward<const Arr>(arr_)[map_[i]]; // preserves const
  }
};

// CTAD: preserve references for lvalues, decay for rvalues
template <class ArrT_, class MapT_>
GenericIndexer(ArrT_ &&, MapT_ &&) -> GenericIndexer<ArrT_, MapT_>;

} // namespace singularity

#endif // #ifndef SINGULARITY_EOS_BASE_GENERIC_INDEXER_HPP_
