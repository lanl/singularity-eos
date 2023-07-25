//------------------------------------------------------------------------------
// © 2021-2023. Triad National Security, LLC. All rights reserved.  This
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

#ifndef SINGULARITY_EOS_BASE_VARIADIC_UTILS_HPP_
#define SINGULARITY_EOS_BASE_VARIADIC_UTILS_HPP_

#include <type_traits>

namespace singularity {
namespace detail {

// variadic list
template <typename... Ts>
struct type_list {};

// variadic list of modifiers
template <template <typename> class... Ts>
struct adapt_list {};

// this flattens a typelist of typelists to a single typelist

// first parameter - accumulator
// second parameter - input list
template <class T, class U>
struct flatten_helper;

// first case - the head of the type_list is type_list too
// expand this type_list and continue
template <class... Ts, class... Heads, class... Tail>
struct flatten_helper<type_list<Ts...>, type_list<type_list<Heads...>, Tail...>> {
  using type =
      typename flatten_helper<type_list<Ts...>, type_list<Heads..., Tail...>>::type;
};

// second case - the head of the type_list is not a type_list
// append it to our new, flattened list
template <class... Ts, class Head, class... Tail>
struct flatten_helper<type_list<Ts...>, type_list<Head, Tail...>> {
  using type = typename flatten_helper<type_list<Ts..., Head>, type_list<Tail...>>::type;
};

// base case - input type_list is empty
// return our flattened list
template <class... Ts>
struct flatten_helper<type_list<Ts...>, type_list<>> {
  using type = type_list<Ts...>;
};

// wrapper around flatten_helper
template <class T>
struct flatten;

// start with an empty accumulator
template <class... Ts>
struct flatten<type_list<Ts...>> {
  using type = typename flatten_helper<type_list<>, type_list<Ts...>>::type;
};

// filter nested variadic templates
template <template <typename> class FIRST, template <typename> class... REST>
constexpr auto remove_first(adapt_list<FIRST, REST...>) {
  return adapt_list<REST...>{};
}

// nested specialization filter type
template <template <typename> class Template, typename T>
struct is_not_duplicate_nested : std::true_type {};

template <template <typename> class Template, typename T>
struct is_not_duplicate_nested<Template, Template<Template<T>>> : std::false_type {};

// filter
template <template <typename> class, template <template <typename> class, class> class,
          template <class...> class, class...>
struct filter;

template <template <typename> class ADAPTER,
          template <template <typename> class, class> class Pred,
          template <class...> class Variadic>
struct filter<ADAPTER, Pred, Variadic> {
  using type = Variadic<>;
};

template <template <typename> class ADAPTER,
          template <template <typename> class, class> class Pred,
          template <class...> class Variadic, class T, class... Ts>
struct filter<ADAPTER, Pred, Variadic, T, Ts...> {
  template <class, class>
  struct Cons;
  template <class Head, class... Tail>
  struct Cons<Head, Variadic<Tail...>> {
    using type = Variadic<Head, Tail...>;
  };

  using type = typename std::conditional<
      Pred<ADAPTER, T>::value,
      typename Cons<T, typename filter<ADAPTER, Pred, Variadic, Ts...>::type>::type,
      typename filter<ADAPTER, Pred, Variadic, Ts...>::type>::type;
};

template <template <typename> class U, typename... Ts>
constexpr auto filter_nested(type_list<Ts...>) {
  using f_list = typename filter<U, is_not_duplicate_nested, type_list, Ts...>::type;
  return f_list{};
}

template <template <typename> class FIRST, template <typename> class... Us,
          typename... Ts>
constexpr auto filter_nested_variadic(adapt_list<FIRST, Us...> m, type_list<Ts...> l) {
  using t1 = type_list<decltype(filter_nested<FIRST>(l))>;
  constexpr typename flatten<t1>::type l1{};
  return filter_nested_variadic(remove_first(m), l1);
}

template <typename... Ts>
constexpr auto filter_nested_variadic(adapt_list<>, type_list<Ts...> l) {
  return l;
}

// apply class template to typelist
template <template <typename> class T, typename... Us>
constexpr auto transform_list(type_list<Us...>) {
  return type_list<T<Us>...>{};
}

template <template <typename> class... Ts, typename... Us>
constexpr auto transform_variadic_list(type_list<Us...> list,
                                       adapt_list<Ts...> mod_list) {
  using t1 = type_list<decltype(transform_list<Ts>(list))...>;
  constexpr typename flatten<t1>::type l1{};
  constexpr auto l2 = filter_nested_variadic(mod_list, l1);
  return l2;
}

// concat typelists
template <typename... S, typename... T>
constexpr auto concat(type_list<S...>, type_list<T...>) {
  return type_list<S..., T...>{};
}

template <typename T, typename... S>
constexpr auto concat(T t, S... ss) {
  return concat(t, concat(ss...));
}

// get number of items in variadic list
template <template <typename> class... Ts>
constexpr auto pack_size(adapt_list<Ts...>) {
  return sizeof...(Ts);
}

template <typename... Ts>
constexpr auto pack_size(type_list<Ts...>) {
  return sizeof...(Ts);
}

} // namespace detail
} // namespace singularity

#endif // SINGULARITY_EOS_BASE_VARIADIC_UTILS_HPP_
