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

#ifndef SINGULARITY_EOS_BASE_VARIADIC_UTILS_HPP_
#define SINGULARITY_EOS_BASE_VARIADIC_UTILS_HPP_

#include <type_traits>
#include <utility>

namespace singularity {
namespace variadic_utils {

// Some generic variatic utilities
// ======================================================================

// Useful for generating nullptr of a specific pointer type
template <typename T>
inline constexpr T *np() {
  return nullptr;
}

// C++14 implementation of std::remove_cvref (available since C++20)
// credit to CJ + Diego
template <typename T>
struct remove_cvref {
  typedef std::remove_cv_t<std::remove_reference_t<T>> type;
};

// Helper types
template <typename T>
using remove_cvref_t = typename remove_cvref<T>::type;

// SFINAE to check if a value is a null ptr
template <typename T, typename = typename std::enable_if<
                          std::is_pointer<remove_cvref_t<T>>::value>::type>
constexpr inline bool is_nullptr(T &&t) {
  return std::forward<T>(t) == nullptr;
}
template <typename T, typename std::enable_if<
                          !std::is_pointer<remove_cvref_t<T>>::value>::type * = nullptr>
constexpr inline bool is_nullptr(T &&) {
  return false;
}

// Backport of C++17 bool_constant.
// With C++17, can be replaced with
// using std::bool_constant
template <bool B>
using bool_constant = std::integral_constant<bool, B>;

// Implementation of std::conjunction/std::disjunction without C++17
// With C++17, can be replaced with
// using std::disjunction
template <bool...>
struct bool_pack {};
template <bool... Bs>
using conjunction = std::is_same<bool_pack<true, Bs...>, bool_pack<Bs..., true>>;
template <bool... Bs>
struct disjunction : bool_constant<!conjunction<!Bs...>::value> {};

// Checks if T is contained in the pack Ts
template <typename T, typename... Ts>
using contains = disjunction<std::is_same<T, Ts>::value...>;

template <typename T, typename... Ts>
constexpr bool contains_v() {
  return contains<T, Ts...>::value;
}

// variadic list
template <typename... Ts>
struct type_list {};

// variadic list of modifiers
template <template <typename> class... Ts>
struct adapt_list {};

// provide index of a type in a type list
template <typename T, typename Head, typename... Ts>
constexpr std::size_t GetIndexInTL(std::size_t current_index) {
  if constexpr (std::is_same_v<T, Head>) {
    return current_index;
  } else {
    static_assert(sizeof...(Ts) > 0, "Type T must be in type list!");
    return GetIndexInTL<T, Ts...>(current_index + 1);
  }
}
template <typename T, typename Head, typename... Ts>
constexpr std::size_t GetIndexInTL() {
  return GetIndexInTL<T, Head, Ts...>(0);
}

// is_indexable similar to is_invokable
template <typename, typename, typename = void>
struct is_indexable : std::false_type {};
template <typename T, typename Index>
struct is_indexable<T, Index,
                    std::void_t<decltype(std::declval<T>()[std::declval<Index>()])>>
    : std::true_type {};
template <typename T, typename Index>
constexpr bool is_indexable_v = is_indexable<T, Index>::value;

// Check if a type can accept an int index
template <class T, class = void>
struct has_int_index : std::false_type {};
template <class T>
struct has_int_index<T, std::void_t<decltype(std::declval<T>()[std::declval<int>()])>>
    : std::true_type {};
template <typename T>
constexpr bool has_int_index_v = has_int_index<T>::value;

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

template <template <typename> class FIRST, template <typename> class... Us,
          typename... Ts>
constexpr auto filter_nested_variadic(adapt_list<FIRST, Us...> m, type_list<Ts...> l) {
  using t1 =
      type_list<typename filter<FIRST, is_not_duplicate_nested, type_list, Ts...>::type>;
  constexpr typename flatten<t1>::type l1{};
  return filter_nested_variadic(adapt_list<Us...>{}, l1);
}

template <typename... Ts>
constexpr auto filter_nested_variadic(adapt_list<>, type_list<Ts...> l) {
  return l;
}

// apply class template to typelist
template <template <typename> class T, typename... Us>
struct transform_list_struct {
  using type = type_list<T<Us>...>;
};

template <template <typename> class... Ts, typename... Us>
constexpr auto transform_variadic_list(type_list<Us...> list,
                                       adapt_list<Ts...> mod_list) {
  using t1 = type_list<typename transform_list_struct<Ts, Us...>::type...>;
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

} // namespace variadic_utils
} // namespace singularity

#endif // SINGULARITY_EOS_BASE_VARIADIC_UTILS_HPP_
