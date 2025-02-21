//------------------------------------------------------------------------------
// © 2025. Triad National Security, LLC. All rights reserved.  This
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

#ifndef SINGULARITY_EOS_BASE_INDEXABLE_TYPES_
#define SINGULARITY_EOS_BASE_INDEXABLE_TYPES_

#include <array>
#include <type_traits>
#include <utility>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/variadic_utils.hpp>

#define SINGULARITY_DECLARE_INDEXABLE_TYPE(TYPE_NAME)                                    \
  struct TYPE_NAME : public singularity::IndexableTypes::IndexableTypesBase<TYPE_NAME> {}

namespace singularity {
namespace IndexableTypes {
// Base class for indexable types that uses CRTP
template <typename CRTP>
struct IndexableTypesBase {
  template <typename, typename = void>
  struct IsOverloadedIn : std::false_type {};
  template <typename T>
  struct IsOverloadedIn<T, std::void_t<decltype(std::declval<T>()[std::declval<CRTP>()])>>
      : std::true_type {};
};
} // namespace IndexableTypes

namespace IndexerUtils {
// Convenience function for accessing an indexer by either type or
// natural number index depending on what is available
template <typename T, typename Indexer_t>
inline auto &Get(Indexer_t &&lambda, std::size_t idx = 0) {
  if constexpr (T::template IsOverloadedIn<Indexer_t>::value) {
    return lambda[T()];
  } else {
    return lambda[idx];
  }
}

// This is a convenience struct to easily build a small indexer with
// a set of indexable types. Note however that accessing with
// indexable types is an O(N) operation with this implementation.
// You probably don't want to use this in real code.
template <typename... Ts>
class VariadicIndexer {
 public:
  VariadicIndexer() = default;
  template <typename T,
            typename = std::enable_if_t<variadic_utils::contains<T, Ts...>::value>>
  PORTABLE_FORCEINLINE_FUNCTION Real &operator[](const T &t) {
    // get the index of T in the variadic type list
    constexpr std::size_t N = sizeof...(Ts);
    std::size_t idx = 0;
    std::size_t i = 0;
    ((std::is_same_v<T, Ts> ? idx = i : ++i), ...);

    return data_[idx];
  }
  PORTABLE_FORCEINLINE_FUNCTION
  Real &operator[](const std::size_t idx) { return data_[idx]; }
  PORTABLE_FORCEINLINE_FUNCTION
  Real &operator[](const int idx) { return data_[idx]; }
  static inline constexpr std::size_t size() { return sizeof...(Ts); }

 private:
  std::array<Real, sizeof...(Ts)> data_;
};

} // namespace IndexerUtils

namespace IndexableTypes {
SINGULARITY_DECLARE_INDEXABLE_TYPE(MeanIonizationState);
SINGULARITY_DECLARE_INDEXABLE_TYPE(LogDensity);
SINGULARITY_DECLARE_INDEXABLE_TYPE(LogTemperature);
SINGULARITY_DECLARE_INDEXABLE_TYPE(MeanAtomicMass);
SINGULARITY_DECLARE_INDEXABLE_TYPE(MeanAtomicNumber);
SINGULARITY_DECLARE_INDEXABLE_TYPE(ElectronFraction);
} // namespace IndexableTypes
} // namespace singularity
#endif // SINGULARITY_EOS_BASE_INDEXABLE_TYPES_
