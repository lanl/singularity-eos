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

// Identifies an indexer as a type-based indexer
template <class, class = void>
struct is_type_indexer : std::false_type {};
template <class Indexer_t>
struct is_type_indexer<Indexer_t,
                       std::void_t<decltype(std::decay_t<Indexer_t>::is_type_indexable)>>
    : std::bool_constant<std::decay_t<Indexer_t>::is_type_indexable> {};
template <class Indexer_t>
constexpr bool is_type_indexer_v = is_type_indexer<Indexer_t>::value;

// The "safe" version of Get(). This function will ONLY return a value IF that
// type-based index is present in the Indexer OR if the Indexer doesn't support
// type-based indexing.
template <typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION
bool safeGet(Indexer_t const &lambda, std::size_t const idx, Real &out) {
  // If null then nothing happens
  if (variadic_utils::is_nullptr(lambda)) {
    return false;
  }

  // Return value if type index is available
  if constexpr (variadic_utils::is_indexable_v<Indexer_t, T>) {
    out = lambda[T{}];
    return true;
  }

  // Do nothing if lambda has type indexing BUT doesn't have this type index
  if constexpr (is_type_indexer_v<Indexer_t>) {
    return false;
  }

  // Fall back to numerical indexing if no type indexing
  if constexpr (variadic_utils::has_whole_num_index<Indexer_t>::value) {
    out = lambda[idx];
    return true;
  }

  // Something else...
  return false;
}

// Overload when no index is provided
template <typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION
bool safeGet(Indexer_t const &lambda, Real &out) {
  // If null then nothing happens
  if (variadic_utils::is_nullptr(lambda)) {
    return false;
  }

  // Return value if type index is available
  if constexpr (variadic_utils::is_indexable_v<Indexer_t, T>) {
    out = lambda[T{}];
    return true;
  }

  // Do nothing if lambda has type indexing BUT doesn't have this type index
  if constexpr (is_type_indexer_v<Indexer_t>) {
    return false;
  }

  // Something else...
  return false;
}

// Same as above but causes an error condition (static or dynamic) if the value
// can't be obtained
template <typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION
Real safeMustGet(Indexer_t const &lambda, std::size_t const idx) {
  // Error on null pointer
  PORTABLE_ALWAYS_REQUIRE(!variadic_utils::is_nullptr(lambda),
                          "Indexer can't be nullptr");

  // Return type-based index. Static assert that type MUST exist in indexer
  if constexpr (is_type_indexer_v<Indexer_t>) {
    static_assert(variadic_utils::is_indexable_v<Indexer_t, T>);
    return lambda[T{}];
  }

  // Fall back to numerical indexing if no type indexing
  if constexpr (variadic_utils::has_whole_num_index<Indexer_t>::value) {
    return lambda[idx];
  }

  // Something else...
  PORTABLE_ALWAYS_THROW_OR_ABORT("Cannot obtain value from unknown indexer");
}

// Break out "Set" functionality from "Get". The original "Get()" did both, but
// the "safe" version needs to separate that functionality for setting the
// values in a lambda
template <typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION
bool safeSet(Indexer_t &lambda, std::size_t const idx, Real const in) {
  // If null then nothing happens
  if (variadic_utils::is_nullptr(lambda)) {
    return false;
  }

  // Return value if type index is available
  if constexpr (variadic_utils::is_indexable_v<Indexer_t, T>) {
    lambda[T{}] = in;
    return true;
  }

  // Do nothing if lambda has type indexing BUT doesn't have this type index
  if constexpr (is_type_indexer_v<Indexer_t>) {
    return false;
  }

  // Fall back to numerical indexing if no type indexing
  if constexpr (variadic_utils::has_whole_num_index<Indexer_t>::value) {
    lambda[idx] = in;
    return true;
  }

  // Something else...
  return false;
}

// Overload without numeric index
template <typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION
bool safeSet(Indexer_t &lambda, Real const in) {
  // If null then nothing happens
  if (variadic_utils::is_nullptr(lambda)) {
    return false;
  }

  // Return value if type index is available
  if constexpr (variadic_utils::is_indexable_v<Indexer_t, T>) {
    lambda[T{}] = in;
    return true;
  }

  // Do nothing if lambda has type indexing BUT doesn't have this type index
  if constexpr (is_type_indexer_v<Indexer_t>) {
    return false;
  }

  // Something else...
  return false;
}

// Same as above but causes an error condition (static or dynamic) if the value
// can't be obtained
template <typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION
void safeMustSet(Indexer_t &lambda, std::size_t const idx, Real const in) {
  // Error on null pointer
  PORTABLE_ALWAYS_REQUIRE(!variadic_utils::is_nullptr(lambda),
                          "Indexer can't be nullptr");

  // Return type-based index. Static assert that type MUST exist in indexer
  if constexpr (is_type_indexer_v<Indexer_t>) {
    static_assert(variadic_utils::is_indexable_v<Indexer_t, T>);
    lambda[T{}] = in;
  }

  // Fall back to numerical indexing if no type indexing
  if constexpr (variadic_utils::has_whole_num_index<Indexer_t>::value) {
    lambda[idx] = in;
  }

  // Something else...
  PORTABLE_ALWAYS_THROW_OR_ABORT("Cannot obtain value from unknown indexer");
}

// NOTE: this Get is "unsafe" because it can allow you to overwrite a type-based
// index since it automatically falls back to numeric indexing if the type
// index isn't present.

// Convenience function for accessing an indexer by either type or natural
// number index depending on what is available.
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
  // Any class that wants to be recognized as indexable (so that we don't
  // accidentally fall back to integer indexing when we don't want to) needs to
  // include this.
  constexpr static bool is_type_indexable = true;

  // JHP: another option for the `is_type_indexable` flag is to take the ADL
  // route. Essentially this would involve defining a friend function that
  // could be defined in an appropriate namesapce so that theoretically a host
  // code could use a TPL container with type-based indexing and allow that
  // container to be flagged in our code as acceptable. This seems like a bit
  // of a heavy hammer for what we need here though. We can easily change this
  // if a TPL provides a type that is being used for this purpose.

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
} // namespace IndexableTypes
} // namespace singularity
#endif // SINGULARITY_EOS_BASE_INDEXABLE_TYPES_
