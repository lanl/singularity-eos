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

#ifndef SINGULARITY_EOS_BASE_TUPLE_UTILS_HPP_
#define SINGULARITY_EOS_BASE_TUPLE_UTILS_HPP_

#include <tuple>

namespace singularity {
namespace tuple_utils {

namespace impl {

// Helper that copies a tuple into a new tuple while also applying a
// transformation function to each element. Wrapped in `impl` namespace due to
// index sequence and called below with a more normal function signature
template <typename TupleT, typename FunctionT, std::size_t... I>
constexpr auto tuple_transform_iteration(TupleT &&tup, FunctionT &&function,
                                         std::index_sequence<I...>) {
  return std::make_tuple(function(std::get<I>(std::forward<TupleT>(tup)))...);
}
} // namespace impl

// Returns a copy of the tuple with a transform function called on each element
template <typename TupleT, typename FunctionT>
constexpr auto tuple_transform(TupleT &&tup, FunctionT &&function) {
  constexpr std::size_t N = std::tuple_size_v<std::remove_reference_t<TupleT>>;
  return impl::tuple_transform_iteration(std::forward<TupleT>(tup),
                                         std::forward<FunctionT>(function),
                                         std::make_index_sequence<N>{});
}

// SFINAE helper to dicriminate between packs and tuples
template <class>
struct is_std_tuple : std::false_type {};
template <class... Ts>
struct is_std_tuple<std::tuple<Ts...>> : std::true_type {};
template <class T>
inline constexpr bool is_std_tuple_v = is_std_tuple<std::decay_t<T>>::value;

} // namespace tuple_utils
} // namespace singularity

#endif // #ifndef SINGULARITY_EOS_BASE_TUPLE_UTILS_HPP_
