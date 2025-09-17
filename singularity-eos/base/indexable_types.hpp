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
#include <ports-of-call/portable_errors.hpp>
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

namespace impl {

// Simple way to switch between pure type indexing or also allowing intergers
enum class AllowedIndexing { Numeric, TypeOnly };

// The "safe" version of Get(). This function will ONLY return a value IF that
// type-based index is present in the Indexer OR if the Indexer doesn't support
// type-based indexing.
template <AllowedIndexing AI, typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION bool SafeGet(Indexer_t const &lambda, std::size_t const idx,
                                           Real &out) {
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

  // Fall back to numeric indexing if allowed
  if constexpr (AI == AllowedIndexing::Numeric) {
    if constexpr (variadic_utils::has_int_index_v<Indexer_t>) {
      out = lambda[idx];
      return true;
    }
  }

  // Something else...
  return false;
}

// Break out "Set" functionality from "Get". The original "Get()" did both, but
// the "safe" version needs to separate that functionality for setting the
// values in a lambda
template <AllowedIndexing AI, typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION bool SafeSet(Indexer_t &lambda, std::size_t const idx,
                                           Real const in) {
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

  // Fall back to numeric indexing if allowed
  if constexpr (AI == AllowedIndexing::Numeric) {
    if constexpr (variadic_utils::has_int_index_v<Indexer_t>) {
      lambda[idx] = in;
      return true;
    }
  }

  // Something else...
  return false;
}

// Same as above but causes an error condition (static or dynamic) if the value
// can't be obtained. Note that the `decltype(auto)` is intended to preserve the
// value category of the square bracket operator of the `Indexer_t` type. This
// allows references to be returned since there is also no possibility of the
// call doing nothing (i.e. like SafeGet and SafeSet), and thus it can be used
// for either setting or getting values. This should also allow for `const`
// correctness downstream in the wrappers where `lambda` is `const`
template <AllowedIndexing AI, typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION decltype(auto) SafeMustGetSet(Indexer_t &&lambda,
                                                            std::size_t const idx) {
  // Error on null pointer
  PORTABLE_ALWAYS_REQUIRE(!variadic_utils::is_nullptr(lambda),
                          "Indexer can't be nullptr");

  // Return type-based index. Static assert that type MUST exist in indexer
  if constexpr (is_type_indexer_v<Indexer_t>) {
    static_assert(variadic_utils::is_indexable_v<Indexer_t, T>);
    // Use std::forward to maintain value category for lambda, and use
    // parentheses to do the same for the output of the lambda[] operation
    return (std::forward<Indexer_t>(lambda)[T{}]);
  }

  // Fall back to numerical indexing if allowed
  if constexpr (AI == AllowedIndexing::Numeric) {
    if constexpr (variadic_utils::has_int_index_v<Indexer_t>) {
      // Use std::forward to maintain value category for lambda, and use
      // parentheses to do the same for the output of the lambda[] operation
      return (std::forward<Indexer_t>(lambda)[idx]);
    }
  }

  // Something else...
  PORTABLE_ALWAYS_THROW_OR_ABORT("Cannot obtain value from unknown indexer");

  // This code is unreachable, but we need a concrete return type so that
  // everything downstream is happy
  Real ret = 0;
  return ret;
}

} // namespace impl

// Overload when numerical index is provided
template <typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION bool SafeGet(Indexer_t const &lambda, std::size_t const idx,
                                           Real &out) {
  return impl::SafeGet<impl::AllowedIndexing::Numeric, T>(lambda, idx, out);
}

// Overload when numerical index isn't provided
template <typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION bool SafeGet(Indexer_t const &lambda, Real &out) {
  std::size_t idx = 0;
  return impl::SafeGet<impl::AllowedIndexing::TypeOnly, T>(lambda, idx, out);
}

// Overload when numerical index is provided
template <typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION bool SafeSet(Indexer_t &lambda, std::size_t const idx,
                                           Real const in) {
  return impl::SafeSet<impl::AllowedIndexing::Numeric, T>(lambda, idx, in);
}

// Overload when numerical index isn't provided
template <typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION bool SafeSet(Indexer_t &lambda, Real const in) {
  std::size_t idx = 0;
  return impl::SafeSet<impl::AllowedIndexing::TypeOnly, T>(lambda, idx, in);
}

// Overload when numerical index is provided
template <typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION Real SafeMustGet(Indexer_t const &lambda,
                                               std::size_t const idx) {
  return impl::SafeMustGetSet<impl::AllowedIndexing::Numeric, T>(lambda, idx);
}

// Overload when numerical index isn't provided
template <typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION Real SafeMustGet(Indexer_t const &lambda) {
  std::size_t idx = 0;
  return impl::SafeMustGetSet<impl::AllowedIndexing::TypeOnly, T>(lambda, idx);
}

// Overload when numerical index is provided
template <typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION void SafeMustSet(Indexer_t &lambda, std::size_t const idx,
                                               Real const in) {
  impl::SafeMustGetSet<impl::AllowedIndexing::Numeric, T>(lambda, idx) = in;
}

// Overload when numerical index isn't provided
template <typename T, typename Indexer_t>
PORTABLE_FORCEINLINE_FUNCTION void SafeMustSet(Indexer_t &lambda, Real const in) {
  std::size_t idx = 0;
  impl::SafeMustGetSet<impl::AllowedIndexing::TypeOnly, T>(lambda, idx) = in;
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
