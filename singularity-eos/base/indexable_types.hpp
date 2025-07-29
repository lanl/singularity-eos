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

namespace singularity {
namespace IndexerUtils {
// Convenience function for accessing an indexer by either type or
// natural number index depending on what is available
template <typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION auto &Get(Indexer_t &&lambda, std::size_t idx = 0) {
  if constexpr (variadic_utils::is_indexable_v<Indexer_t, T>) {
    return lambda[T()];
  } else {
    return lambda[idx];
  }
}

// This is a convenience struct to easily build a small indexer with
// a set of indexable types.
template <typename Data_t, typename... Ts>
class VariadicIndexerBase {
 public:
  VariadicIndexerBase() = default;
  PORTABLE_FORCEINLINE_FUNCTION
  VariadicIndexerBase(const Data_t &data) : data_(data) {}
  template <typename T,
            typename = std::enable_if_t<variadic_utils::contains<T, Ts...>::value>>
  PORTABLE_FORCEINLINE_FUNCTION Real &operator[](const T &t) {
    constexpr std::size_t idx = variadic_utils::GetIndexInTL<T, Ts...>();
    return data_[idx];
  }
  PORTABLE_FORCEINLINE_FUNCTION
  Real &operator[](const std::size_t idx) { return data_[idx]; }
  template <typename T,
            typename = std::enable_if_t<variadic_utils::contains<T, Ts...>::value>>
  PORTABLE_FORCEINLINE_FUNCTION const Real &operator[](const T &t) const {
    constexpr std::size_t idx = variadic_utils::GetIndexInTL<T, Ts...>();
    return data_[idx];
  }
  PORTABLE_FORCEINLINE_FUNCTION
  const Real &operator[](const std::size_t idx) const { return data_[idx]; }
  static inline constexpr std::size_t size() { return sizeof...(Ts); }

 private:
  Data_t data_;
};
// uses an array
template <typename... Ts>
using VariadicIndexer = VariadicIndexerBase<std::array<Real, sizeof...(Ts)>, Ts...>;
// uses a Real*
template <typename... Ts>
using VariadicPointerIndexer = VariadicIndexerBase<Real *, Ts...>;
} // namespace IndexerUtils

namespace IndexableTypes {
struct MeanIonizationState {};
struct LogDensity {};
struct LogTemperature {};
struct MeanAtomicMass {};
struct MeanAtomicNumber {};
struct ElectronFraction {};
struct RootStatus {};
struct TableStatus {};
template<int mat_idx>
struct MassFraction {};
} // namespace IndexableTypes
} // namespace singularity
#endif // SINGULARITY_EOS_BASE_INDEXABLE_TYPES_
