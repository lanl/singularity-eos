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

#include <tuple>

namespace singularity {
namespace tuple_utils {

// Should this actually live in variadic utils?
namespace impl {

// Helper to transform a tuple into a new tuple by applying a function on each element.
template <typename TupleT, typename FunctionT, std::size_t... I>
constexpr auto tuple_transform_iteration(TupleT&& tup, FunctionT function,
                                         std::index_sequence<I...>) {
    return std::make_tuple(function(std::get<I>(std::forward<TupleT>(tup)))...);
}
} // namespace impl

// Returns a copy of the tuple with a transform function called on each element
template <typename TupleT, typename FunctionT>
constexpr auto tuple_transform(TupleT&& tup, FunctionT function) {
    constexpr std::size_t N = std::tuple_size_v<std::remove_reference_t<TupleT>>;
    return impl::tuple_transform_iteration(std::forward<TupleT>(tup), function,
                                           std::make_index_sequence<N>{});
}

} // namespace tuple_utils
} // namespace singularity
