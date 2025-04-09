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

#ifndef _SINGULARITY_EOS_EOS_EOS_BASE_
#define _SINGULARITY_EOS_EOS_EOS_BASE_

#include <cstring>
#include <limits>
#include <string>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>
#include <singularity-eos/base/variadic_utils.hpp>

namespace singularity {
namespace mfuncname {
static inline auto member_func_name(const char *type_name, const char *func_name) {
  return std::string(type_name) + std::string("::") + std::string(func_name);
}
} // namespace mfuncname
} // namespace singularity
#define SG_MEMBER_FUNC_NAME()                                                            \
  singularity::mfuncname::member_func_name(typeid(CRTP).name(), __func__);

namespace singularity {
namespace eos_base {

namespace impl {
constexpr std::size_t MAX_NUM_CHARS = 121;
// Cuda doesn't have strcat, so we implement it ourselves
PORTABLE_FORCEINLINE_FUNCTION
char *StrCat(char *destination, const char *source) {
  std::size_t i, j; // not in loops because they're re-used.

  // specifically avoid strlen, which isn't on GPU
  for (i = 0; destination[i] != '\0'; i++) {
  }
  // assumes destination has enough memory allocated
  for (j = 0; source[j] != '\0'; j++) {
    // MAX_NUM_CHARS-1 to leave room for null terminator
    std::size_t ipj = i + j;
    PORTABLE_REQUIRE(ipj < MAX_NUM_CHARS - 1,
                     "Concat string must be within allowed size");
    destination[ipj] = source[j];
  }
  // null terminate destination string
  destination[i + j] = '\0';

  // the destination is returned by standard `strcat()`
  return destination;
}
} // namespace impl

// This Macro adds the `using` statements that allow for the base class
// VECTOR functionality to overload the scalar implementations in the derived
// classes. Do not add functions here that are not overloads of derived class features.
// TODO(JMM): Should we have more macros that capture just some of these?
// JMM: Use VA_ARGS to capture more complex template types
#define SG_ADD_BASE_CLASS_USINGS(...)                                                    \
  using EosBase<__VA_ARGS__>::TemperatureFromDensityInternalEnergy;                      \
  using EosBase<__VA_ARGS__>::InternalEnergyFromDensityTemperature;                      \
  using EosBase<__VA_ARGS__>::PressureFromDensityTemperature;                            \
  using EosBase<__VA_ARGS__>::PressureFromDensityInternalEnergy;                         \
  using EosBase<__VA_ARGS__>::MinInternalEnergyFromDensity;                              \
  using EosBase<__VA_ARGS__>::SpecificHeatFromDensityTemperature;                        \
  using EosBase<__VA_ARGS__>::SpecificHeatFromDensityInternalEnergy;                     \
  using EosBase<__VA_ARGS__>::BulkModulusFromDensityTemperature;                         \
  using EosBase<__VA_ARGS__>::BulkModulusFromDensityInternalEnergy;                      \
  using EosBase<__VA_ARGS__>::GruneisenParamFromDensityTemperature;                      \
  using EosBase<__VA_ARGS__>::GruneisenParamFromDensityInternalEnergy;                   \
  using EosBase<__VA_ARGS__>::FillEos;                                                   \
  using EosBase<__VA_ARGS__>::EntropyFromDensityTemperature;                             \
  using EosBase<__VA_ARGS__>::EntropyFromDensityInternalEnergy;                          \
  using EosBase<__VA_ARGS__>::GibbsFreeEnergyFromDensityTemperature;                     \
  using EosBase<__VA_ARGS__>::GibbsFreeEnergyFromDensityInternalEnergy;                  \
  using EosBase<__VA_ARGS__>::scratch_size;                                              \
  using EosBase<__VA_ARGS__>::max_scratch_size;

// This macro adds these methods to a derived class. Due to scope,
// these can't be implemented in the base class, unless we make
// _AZbar public. Not all EOS's may want these default functions
// TODO(JMM): Should we go the alternate route and make _AZbar public?
#define SG_ADD_DEFAULT_MEAN_ATOMIC_FUNCTIONS(_AZbar)                                     \
  PORTABLE_INLINE_FUNCTION                                                               \
  Real MeanAtomicMass() const { return _AZbar.Abar; }                                    \
  PORTABLE_INLINE_FUNCTION                                                               \
  Real MeanAtomicNumber() const { return _AZbar.Zbar; }

// This macro adds several methods that most modifiers will
// want. Not ALL modifiers will want these methods as written here,
// so use this macro with care.
// TODO(JMM): Find a better solution. Multiple inheritence and mixins
// dont' seem to work as desired here.
#define SG_ADD_MODIFIER_METHODS(T, t_)                                                   \
  static inline constexpr bool IsModified() { return true; }                             \
  inline constexpr T UnmodifyOnce() { return t_; }                                       \
  std::size_t DynamicMemorySizeInBytes() const { return t_.DynamicMemorySizeInBytes(); } \
  std::size_t SharedMemorySizeInBytes() const { return t_.SharedMemorySizeInBytes(); }   \
  std::size_t DumpDynamicMemory(char *dst) { return t_.DumpDynamicMemory(dst); }         \
  std::size_t SetDynamicMemory(char *src,                                                \
                               const SharedMemSettings &stngs = DEFAULT_SHMEM_STNGS) {   \
    return t_.SetDynamicMemory(src, stngs);                                              \
  }                                                                                      \
  constexpr bool AllDynamicMemoryIsShareable() const {                                   \
    return t_.AllDynamicMemoryIsShareable();                                             \
  }

#define SG_ADD_MODIFIER_MEAN_METHODS(t_)                                                 \
  PORTABLE_INLINE_FUNCTION                                                               \
  Real MeanAtomicMass() const { return t_.MeanAtomicMass(); }                            \
  PORTABLE_INLINE_FUNCTION                                                               \
  Real MeanAtomicNumber() const { return t_.MeanAtomicNumber(); }

class Factor {
  Real value_ = 1.0;
  bool is_set_ = false;

 public:
  bool is_set() const { return is_set_; }

  Real get() const { return value_; }

  void set(Real v) {
    is_set_ = true;
    value_ = v;
  }

  void apply(Real v) {
    is_set_ = true;
    value_ *= v;
  }

  void clear() {
    is_set_ = false;
    value_ = 1.0;
  }
};

struct Transform {
  Factor x, y, f;
};

/*
  This is a utility struct used to bundle mean atomic
  mass/number. Used in the default implementations of MeanAtomicMass
  and MeanAtomicNumber provided by the base class.
 */
struct MeanAtomicProperties {
  Real Abar, Zbar;

  // default is hydrogen
  static constexpr Real DEFAULT_ABAR = 1.0;
  static constexpr Real DEFAULT_ZBAR = 1.0;

  PORTABLE_INLINE_FUNCTION
  MeanAtomicProperties(Real Abar_, Real Zbar_) : Abar(Abar_), Zbar(Zbar_) {}
  PORTABLE_INLINE_FUNCTION
  MeanAtomicProperties() : Abar(DEFAULT_ABAR), Zbar(DEFAULT_ZBAR) {}
  PORTABLE_INLINE_FUNCTION
  void CheckParams() const {
    PORTABLE_ALWAYS_REQUIRE(Abar > 0, "Positive mean atomic mass");
    PORTABLE_ALWAYS_REQUIRE(Zbar > 0, "Positive mean atomic number");
  }
  PORTABLE_INLINE_FUNCTION
  void PrintParams() const {
    printf("      Abar  = %g\n", Abar);
    printf("      Zbar  = %g\n", Zbar);
  }
};

/*
This is a CRTP that allows for static inheritance so that default behavior for
various member functions can be defined.
https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern

In particular, the default behavior for the vector version of the EOS lookup
member functions is to perform a `portableFor` loop over all of the input states
*/
template <typename CRTP>
class EosBase {
 public:
  template <typename T, typename R>
  struct is_raw_pointer
      : std::is_same<std::remove_reference_t<std::remove_cv_t<T>>, R *> {};

  // Generic evaluator
  template <typename Functor_t>
  PORTABLE_INLINE_FUNCTION void EvaluateDevice(const Functor_t f) const {
    const CRTP copy = *(static_cast<CRTP const *>(this));
    f(copy);
  }
  template <typename Functor_t>
  void EvaluateHost(Functor_t &f) const {
    const CRTP copy = *(static_cast<CRTP const *>(this));
    f(copy);
  }

  // EOS builder helpers
  // Checks if an EOS can be modified
  template <template <class> typename Mod, typename... Ts>
  constexpr bool ModifiedInList() const {
    return variadic_utils::contains_v<Mod<CRTP>, Ts...>();
  }
  // Modifies an EOS object by pulling out underlying type and modifying it
  // This one returns Mod<CRT>
  template <template <class> typename Mod, typename... Args>
  constexpr Mod<CRTP> Modify(Args &&...args) const {
    CRTP unmodified = *(static_cast<CRTP const *>(this));
    return Mod<CRTP>(std::move(unmodified), std::forward<Args>(args)...);
  }
  // These are overloads needed for the variant, as std::visit must be
  // able to return a variant of the same type every time. This lets
  // us do so, even though sometimes we don't modify the object.
  template <template <class> typename Mod, typename... Args>
  constexpr Mod<CRTP> ConditionallyModify(std::true_type, Args &&...args) const {
    return Modify<Mod>(std::forward<Args>(args)...);
  }
  template <template <class> typename Mod, typename... Args>
  constexpr CRTP ConditionallyModify(std::false_type, Args &&...args) const {
    CRTP unmodified = *(static_cast<CRTP const *>(this));
    return unmodified;
  }
  template <template <class> typename Mod, typename... Ts, typename... Args>
  constexpr auto ConditionallyModify(const variadic_utils::type_list<Ts...> &tl,
                                     Args &&...args) const {
    constexpr bool do_mod = variadic_utils::contains_v<Mod<CRTP>, Ts...>();
    return ConditionallyModify<Mod>(variadic_utils::bool_constant<do_mod>(),
                                    std::forward<Args>(args)...);
  }

  // Scalar member functions that get shared
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GibbsFreeEnergyFromDensityTemperature(
      const Real rho, const Real T,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const CRTP copy = *(static_cast<CRTP const *>(this));
    Real sie = copy.InternalEnergyFromDensityTemperature(rho, T, lambda);
    Real P = copy.PressureFromDensityTemperature(rho, T, lambda);
    Real S = copy.EntropyFromDensityTemperature(rho, T, lambda);
    return sie + (P / rho) - T * S;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GibbsFreeEnergyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const CRTP copy = *(static_cast<CRTP const *>(this));
    Real T = copy.TemperatureFromDensityInternalEnergy(rho, sie, lambda);
    Real P = copy.PressureFromDensityTemperature(rho, T, lambda);
    Real S = copy.EntropyFromDensityTemperature(rho, T, lambda);
    return sie + (P / rho) - T * S;
  }

  // Vector member functions
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                       RealIndexer &&temperatures, const int num,
                                       LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          temperatures[i] =
              copy.TemperatureFromDensityInternalEnergy(rhos[i], sies[i], lambdas[i]);
        });
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void
  TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                       RealIndexer &&temperatures, Real * /*scratch*/,
                                       const int num, LambdaIndexer &&lambdas) const {
    TemperatureFromDensityInternalEnergy(std::forward<ConstRealIndexer>(rhos),
                                         std::forward<ConstRealIndexer>(sies),
                                         std::forward<RealIndexer>(temperatures), num,
                                         std::forward<LambdaIndexer>(lambdas));
  }

  template <typename LambdaIndexer>
  inline void TemperatureFromDensityInternalEnergy(const Real *rhos, const Real *sies,
                                                   Real *temperatures, Real * /*scratch*/,
                                                   const int num, LambdaIndexer &&lambdas,
                                                   Transform && = Transform()) const {
    TemperatureFromDensityInternalEnergy(rhos, sies, temperatures, num,
                                         std::forward<LambdaIndexer>(lambdas));
  }

  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&sies, const int num,
                                                   LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          sies[i] = copy.InternalEnergyFromDensityTemperature(rhos[i], temperatures[i],
                                                              lambdas[i]);
        });
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&sies, Real * /*scratch*/,
                                                   const int num,
                                                   LambdaIndexer &&lambdas) const {
    InternalEnergyFromDensityTemperature(std::forward<ConstRealIndexer>(rhos),
                                         std::forward<ConstRealIndexer>(temperatures),
                                         std::forward<RealIndexer>(sies), num,
                                         std::forward<LambdaIndexer>(lambdas));
  }
  template <typename LambdaIndexer>
  inline void InternalEnergyFromDensityTemperature(const Real *rhos,
                                                   const Real *temperatures, Real *sies,
                                                   Real * /*scratch*/, const int num,
                                                   LambdaIndexer &&lambdas,
                                                   Transform && = Transform()) const {
    InternalEnergyFromDensityTemperature(rhos, temperatures, sies, num,
                                         std::forward<LambdaIndexer>(lambdas));
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void PressureFromDensityTemperature(ConstRealIndexer &&rhos,
                                             ConstRealIndexer &&temperatures,
                                             RealIndexer &&pressures, const int num,
                                             LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          pressures[i] =
              copy.PressureFromDensityTemperature(rhos[i], temperatures[i], lambdas[i]);
        });
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void
  PressureFromDensityTemperature(ConstRealIndexer &&rhos, ConstRealIndexer &&temperatures,
                                 RealIndexer &&pressures, Real * /*scratch*/,
                                 const int num, LambdaIndexer &&lambdas) const {
    PressureFromDensityTemperature(std::forward<ConstRealIndexer>(rhos),
                                   std::forward<ConstRealIndexer>(temperatures),
                                   std::forward<RealIndexer>(pressures), num,
                                   std::forward<LambdaIndexer>(lambdas));
  }
  template <typename LambdaIndexer>
  inline void PressureFromDensityTemperature(const Real *rhos, const Real *temperatures,
                                             Real *pressures, Real * /*scratch*/,
                                             const int num, LambdaIndexer &&lambdas,
                                             Transform && = Transform()) const {
    PressureFromDensityTemperature(rhos, temperatures, pressures, num,
                                   std::forward<LambdaIndexer>(lambdas));
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                ConstRealIndexer &&sies,
                                                RealIndexer &&pressures, const int num,
                                                LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          pressures[i] =
              copy.PressureFromDensityInternalEnergy(rhos[i], sies[i], lambdas[i]);
        });
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void
  PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                    RealIndexer &&pressures, Real * /*scratch*/,
                                    const int num, LambdaIndexer &&lambdas) const {
    PressureFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(pressures), num, std::forward<LambdaIndexer>(lambdas));
  }
  template <typename LambdaIndexer>
  inline void PressureFromDensityInternalEnergy(const Real *rhos, const Real *sies,
                                                Real *pressures, Real * /*scratch*/,
                                                const int num, LambdaIndexer &&lambdas,
                                                Transform && = Transform()) const {
    PressureFromDensityInternalEnergy(rhos, sies, pressures, num,
                                      std::forward<LambdaIndexer>(lambdas));
  }
  ///
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void MinInternalEnergyFromDensity(ConstRealIndexer &&rhos, RealIndexer &&sies,
                                           const int num, LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          sies[i] = copy.MinInternalEnergyFromDensity(rhos[i], lambdas[i]);
        });
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void MinInternalEnergyFromDensity(ConstRealIndexer &&rhos, RealIndexer &&sies,
                                           Real * /*scratch*/, const int num,
                                           LambdaIndexer &&lambdas) const {
    MinInternalEnergyFromDensity(std::forward<ConstRealIndexer>(rhos),
                                 std::forward<RealIndexer>(sies), num,
                                 std::forward<LambdaIndexer>(lambdas));
  }
  template <typename LambdaIndexer>
  inline void MinInternalEnergyFromDensity(const Real *rhos, Real *sies,
                                           Real * /*scratch*/, const int num,
                                           LambdaIndexer &&lambdas,
                                           Transform && = Transform()) const {
    MinInternalEnergyFromDensity(rhos, num, std::forward<LambdaIndexer>(lambdas));
  }
  ///
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void EntropyFromDensityTemperature(ConstRealIndexer &&rhos,
                                            ConstRealIndexer &&temperatures,
                                            RealIndexer &&entropies, const int num,
                                            LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          entropies[i] =
              copy.EntropyFromDensityTemperature(rhos[i], temperatures[i], lambdas[i]);
        });
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void
  EntropyFromDensityTemperature(ConstRealIndexer &&rhos, ConstRealIndexer &&temperatures,
                                RealIndexer &&entropies, Real * /*scratch*/,
                                const int num, LambdaIndexer &&lambdas) const {
    EntropyFromDensityTemperature(std::forward<ConstRealIndexer>(rhos),
                                  std::forward<ConstRealIndexer>(temperatures),
                                  std::forward<RealIndexer>(entropies), num,
                                  std::forward<LambdaIndexer>(lambdas));
  }
  template <typename LambdaIndexer>
  inline void EntropyFromDensityTemperature(const Real *rhos, const Real *temperatures,
                                            Real *entropies, Real * /*scratch*/,
                                            const int num, LambdaIndexer &&lambdas,
                                            Transform && = Transform()) const {
    EntropyFromDensityTemperature(rhos, temperatures, entropies, num,
                                  std::forward<LambdaIndexer>(lambdas));
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void EntropyFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                               ConstRealIndexer &&sies,
                                               RealIndexer &&entropies, const int num,
                                               LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          entropies[i] =
              copy.EntropyFromDensityInternalEnergy(rhos[i], sies[i], lambdas[i]);
        });
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void
  EntropyFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                   RealIndexer &&entropies, Real * /*scratch*/,
                                   const int num, LambdaIndexer &&lambdas) const {
    EntropyFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(entropies), num, std::forward<LambdaIndexer>(lambdas));
  }
  template <typename LambdaIndexer>
  inline void EntropyFromDensityInternalEnergy(const Real *rhos, const Real *sies,
                                               Real *entropies, Real * /*scratch*/,
                                               const int num, LambdaIndexer &&lambdas,
                                               Transform && = Transform()) const {
    EntropyFromDensityInternalEnergy(rhos, sies, entropies, num,
                                     std::forward<LambdaIndexer>(lambdas));
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                                 ConstRealIndexer &&temperatures,
                                                 RealIndexer &&cvs, const int num,
                                                 LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          cvs[i] = copy.SpecificHeatFromDensityTemperature(rhos[i], temperatures[i],
                                                           lambdas[i]);
        });
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                                 ConstRealIndexer &&temperatures,
                                                 RealIndexer &&cvs, Real * /*scratch*/,
                                                 const int num,
                                                 LambdaIndexer &&lambdas) const {
    SpecificHeatFromDensityTemperature(std::forward<ConstRealIndexer>(rhos),
                                       std::forward<ConstRealIndexer>(temperatures),
                                       std::forward<RealIndexer>(cvs), num,
                                       std::forward<LambdaIndexer>(lambdas));
  }
  template <typename LambdaIndexer>
  inline void SpecificHeatFromDensityTemperature(const Real *rhos,
                                                 const Real *temperatures, Real *cvs,
                                                 Real * /*scratch*/, const int num,
                                                 LambdaIndexer &&lambdas,
                                                 Transform && = Transform()) const {
    SpecificHeatFromDensityTemperature(rhos, temperatures, cvs, num,
                                       std::forward<LambdaIndexer>(lambdas));
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                    ConstRealIndexer &&sies,
                                                    RealIndexer &&cvs, const int num,
                                                    LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          cvs[i] =
              copy.SpecificHeatFromDensityInternalEnergy(rhos[i], sies[i], lambdas[i]);
        });
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void
  SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                        RealIndexer &&cvs, Real * /*scratch*/,
                                        const int num, LambdaIndexer &&lambdas) const {
    SpecificHeatFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(cvs), num, std::forward<LambdaIndexer>(lambdas));
  }
  template <typename LambdaIndexer>
  inline void SpecificHeatFromDensityInternalEnergy(const Real *rhos, const Real *sies,
                                                    Real *cvs, Real * /*scratch*/,
                                                    const int num,
                                                    LambdaIndexer &&lambdas,
                                                    Transform && = Transform()) const {
    SpecificHeatFromDensityInternalEnergy(rhos, sies, cvs, num,
                                          std::forward<LambdaIndexer>(lambdas));
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                                ConstRealIndexer &&temperatures,
                                                RealIndexer &&bmods, const int num,
                                                LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          bmods[i] = copy.BulkModulusFromDensityTemperature(rhos[i], temperatures[i],
                                                            lambdas[i]);
        });
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                                ConstRealIndexer &&temperatures,
                                                RealIndexer &&bmods, Real * /*scratch*/,
                                                const int num,
                                                LambdaIndexer &&lambdas) const {
    BulkModulusFromDensityTemperature(std::forward<ConstRealIndexer>(rhos),
                                      std::forward<ConstRealIndexer>(temperatures),
                                      std::forward<RealIndexer>(bmods), num,
                                      std::forward<LambdaIndexer>(lambdas));
  }
  template <typename LambdaIndexer>
  inline void BulkModulusFromDensityTemperature(const Real *rhos,
                                                const Real *temperatures, Real *bmods,
                                                Real * /*scratch*/, const int num,
                                                LambdaIndexer &&lambdas,
                                                Transform && = Transform()) const {
    BulkModulusFromDensityTemperature(rhos, temperatures, bmods, num,
                                      std::forward<LambdaIndexer>(lambdas));
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&sies,
                                                   RealIndexer &&bmods, const int num,
                                                   LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          bmods[i] =
              copy.BulkModulusFromDensityInternalEnergy(rhos[i], sies[i], lambdas[i]);
        });
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void
  BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                       RealIndexer &&bmods, Real * /*scratch*/,
                                       const int num, LambdaIndexer &&lambdas) const {
    BulkModulusFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(bmods), num, std::forward<LambdaIndexer>(lambdas));
  }
  template <typename LambdaIndexer>
  inline void BulkModulusFromDensityInternalEnergy(const Real *rhos, const Real *sies,
                                                   Real *bmods, Real * /*scratch*/,
                                                   const int num, LambdaIndexer &&lambdas,
                                                   Transform && = Transform()) const {
    BulkModulusFromDensityInternalEnergy(rhos, sies, bmods, num,
                                         std::forward<LambdaIndexer>(lambdas));
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&gm1s, const int num,
                                                   LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          gm1s[i] = copy.GruneisenParamFromDensityTemperature(rhos[i], temperatures[i],
                                                              lambdas[i]);
        });
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&gm1s, Real * /*scratch*/,
                                                   const int num,
                                                   LambdaIndexer &&lambdas) const {
    GruneisenParamFromDensityTemperature(std::forward<ConstRealIndexer>(rhos),
                                         std::forward<ConstRealIndexer>(temperatures),
                                         std::forward<RealIndexer>(gm1s), num,
                                         std::forward<LambdaIndexer>(lambdas));
  }
  template <typename LambdaIndexer>
  inline void GruneisenParamFromDensityTemperature(const Real *rhos,
                                                   const Real *temperatures, Real *gm1s,
                                                   Real * /*scratch*/, const int num,
                                                   LambdaIndexer &&lambdas,
                                                   Transform && = Transform()) const {
    GruneisenParamFromDensityTemperature(rhos, temperatures, gm1s, num,
                                         std::forward<LambdaIndexer>(lambdas));
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                      ConstRealIndexer &&sies,
                                                      RealIndexer &&gm1s, const int num,
                                                      LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          gm1s[i] =
              copy.GruneisenParamFromDensityInternalEnergy(rhos[i], sies[i], lambdas[i]);
        });
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                      ConstRealIndexer &&sies,
                                                      RealIndexer &&gm1s,
                                                      Real * /*scratch*/, const int num,
                                                      LambdaIndexer &&lambdas) const {
    GruneisenParamFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(gm1s), num, std::forward<LambdaIndexer>(lambdas));
  }
  template <typename LambdaIndexer>
  inline void GruneisenParamFromDensityInternalEnergy(const Real *rhos, const Real *sies,
                                                      Real *gm1s, Real * /*scratch*/,
                                                      const int num,
                                                      LambdaIndexer &&lambdas,
                                                      Transform && = Transform()) const {
    GruneisenParamFromDensityInternalEnergy(rhos, sies, gm1s, num,
                                            std::forward<LambdaIndexer>(lambdas));
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void GibbsFreeEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                                    ConstRealIndexer &&Ts,
                                                    RealIndexer &&Gs, const int num,
                                                    LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          Gs[i] = copy.GibbsFreeEnergyFromDensityTemperature(rhos[i], Ts[i], lambdas[i]);
        });
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void
  GibbsFreeEnergyFromDensityTemperature(ConstRealIndexer &&rhos, ConstRealIndexer &&Ts,
                                        RealIndexer &&Gs, Real * /*scratch*/,
                                        const int num, LambdaIndexer &&lambdas) const {
    GibbsFreeEnergyFromDensityTemperature(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(Ts),
        std::forward<RealIndexer>(Gs), num, std::forward<LambdaIndexer>(lambdas));
  }
  template <typename LambdaIndexer>
  inline void GibbsFreeEnergyFromDensityTemperature(const Real *rhos, const Real *Ts,
                                                    Real *Gs, Real * /*scratch*/,
                                                    const int num,
                                                    LambdaIndexer &&lambdas,
                                                    Transform && = Transform()) const {
    GibbsFreeEnergyFromDensityTemperature(rhos, Ts, Gs, num,
                                          std::forward<LambdaIndexer>(lambdas));
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void GibbsFreeEnergyFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                       ConstRealIndexer &&sies,
                                                       RealIndexer &&Gs, const int num,
                                                       LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          Gs[i] =
              copy.GibbsFreeEnergyFromDensityInternalEnergy(rhos[i], sies[i], lambdas[i]);
        });
  }
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer,
            typename = std::enable_if_t<!is_raw_pointer<RealIndexer, Real>::value>>
  inline void GibbsFreeEnergyFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                       ConstRealIndexer &&sies,
                                                       RealIndexer &&Gs,
                                                       Real * /*scratch*/, const int num,
                                                       LambdaIndexer &&lambdas) const {
    GibbsFreeEnergyFromDensityInternalEnergy(
        std::forward<ConstRealIndexer>(rhos), std::forward<ConstRealIndexer>(sies),
        std::forward<RealIndexer>(Gs), num, std::forward<LambdaIndexer>(lambdas));
  }
  template <typename LambdaIndexer>
  inline void GibbsFreeEnergyFromDensityInternalEnergy(const Real *rhos, const Real *sies,
                                                       Real *Gs, Real * /*scratch*/,
                                                       const int num,
                                                       LambdaIndexer &&lambdas,
                                                       Transform && = Transform()) const {
    GibbsFreeEnergyFromDensityInternalEnergy(rhos, sies, Gs, num,
                                             std::forward<LambdaIndexer>(lambdas));
  }

  template <typename RealIndexer, typename LambdaIndexer>
  inline void FillEos(RealIndexer &&rhos, RealIndexer &&temps, RealIndexer &&energies,
                      RealIndexer &&presses, RealIndexer &&cvs, RealIndexer &&bmods,
                      const int num, const unsigned long output,
                      LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          copy.FillEos(rhos[i], temps[i], energies[i], presses[i], cvs[i], bmods[i],
                       output, lambdas[i]);
        });
  }

  // Report minimum values of density and temperature
  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumDensity() const { return 0; }
  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumTemperature() const { return 0; }

  // Report maximum value of density. Default is unbounded.
  // JMM: Should we use actual infinity, the largest real, or just a
  // big number?  For comparisons, actual infinity is better. It also
  // has the advantage of being projective with modifiers that modify
  // the max. On the other hand, it's more fraught if someone tries to
  // put it into a formula without guarding against it.
  PORTABLE_FORCEINLINE_FUNCTION
  Real MaximumDensity() const { return 1e100; }

  // These are for the PT space PTE solver to bound the iterations in
  // a safe range.
  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumPressure() const { return 0; }
  // Gruneisen EOS's often have a maximum density, which implies a maximum pressure.
  PORTABLE_FORCEINLINE_FUNCTION
  Real MaximumPressureAtTemperature([[maybe_unused]] const Real T) const { return 1e100; }

  PORTABLE_INLINE_FUNCTION
  Real RhoPmin(const Real temp) const { return 0.0; }

  static inline unsigned long scratch_size(const std::string method,
                                           const unsigned int nelements) {
    return 0;
  }
  static inline unsigned long max_scratch_size(const unsigned int nelements) { return 0; }
  PORTABLE_INLINE_FUNCTION
  static int nlambda() noexcept { return 0.; }

  // JMM: EOS's which encapsulate a mix or reactions may wish to vary
  // this.  For example, Helmholtz and StellarCollapse. This isn't the
  // default, so by default the base class provides a specialization.
  // for models where density and temperature are required, the EOS
  // developer is in charge of either throwing an error or choosing
  // reasonable defaults.
  // TODO(JMM): Should we provide vector implementations if we depend
  // on rho, T, etc?
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicMassFromDensityTemperature(
      const Real rho, const Real T,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    CRTP copy = *(static_cast<CRTP const *>(this));
    return copy.MeanAtomicMass();
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicNumberFromDensityTemperature(
      const Real rho, const Real T,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    CRTP copy = *(static_cast<CRTP const *>(this));
    return copy.MeanAtomicNumber();
  }

  // Default entropy behavior is to cause an error
  PORTABLE_FORCEINLINE_FUNCTION
  void EntropyIsNotEnabled(const char *eosname) const {
    // Construct the error message using char* so it works on device
    // WARNING: This needs to be updated if EOS names get longer
    // base msg length 32 + 5 chars = 37 chars
    // + 1 char for null terminator
    // maximum allowed EOS length = 44 chars
    char msg[impl::MAX_NUM_CHARS] = "Singularity-EOS: Entropy is not enabled for the '";
    impl::StrCat(msg, eosname);
    impl::StrCat(msg, "' EOS");
    PORTABLE_ALWAYS_THROW_OR_ABORT(msg);
  }

  // Default MinInternalEnergyFromDensity behavior is to just return the zero-K isotherm.
  // This should be fine for all thermodynamically consistent EOS, but could cause issues
  // with EOS that aren't thermodynamically consistent.
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MinInternalEnergyFromDensity(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    CRTP copy = *(static_cast<CRTP const *>(this));
    return copy.InternalEnergyFromDensityTemperature(rho, 0.);
  }

  // This error is useful for EOS where the zero-K approximation is invalid for whatever
  // reason
  PORTABLE_FORCEINLINE_FUNCTION
  void MinInternalEnergyIsNotEnabled(const char *eosname) const {
    // Construct the error message using char* so it works on device
    // WARNING: This needs to be updated if EOS names get longer
    // base msg length 32 + 5 chars = 37 chars
    // + 1 char for null terminator
    // maximum allowed EOS length = 44 chars
    char msg[impl::MAX_NUM_CHARS] =
        "Singularity-EOS: MinInternalEnergyFromDensity() is not enabled for the '";
    impl::StrCat(msg, eosname);
    impl::StrCat(msg, "' EOS");
    PORTABLE_ALWAYS_THROW_OR_ABORT(msg);
  }

  // JMM: This method is often going to be overloaded for special cases.
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const {
    using RootFinding1D::findRoot; // more robust but slower. Better default.
    using RootFinding1D::Status;

    // Pressure is not monotone in density at low densities, which can
    // prevent convergence. We want to approach tension from above,
    // not below. Choose close to, but above, normal density for a
    // metal like copper.
    constexpr Real DEFAULT_RHO_GUESS = 12;

    CRTP copy = *(static_cast<CRTP const *>(this));

    // P(rho) not monotone. When relevant, bound rhopmin.
    Real rhomin = std::max(copy.RhoPmin(temp), copy.MinimumDensity());
    Real rhomax = copy.MaximumDensity();
    PORTABLE_REQUIRE(rhomax > rhomin, "max bound > min bound");

    auto PofRT = [&](const Real r) {
      return copy.PressureFromDensityTemperature(r, temp, lambda);
    };
    Real rhoguess = rho;                                // use input density
    if ((rhoguess <= rhomin) || (rhoguess >= rhomax)) { // avoid edge effects
      if ((rhomin < DEFAULT_RHO_GUESS) && (DEFAULT_RHO_GUESS < rhomax)) {
        rhoguess = DEFAULT_RHO_GUESS;
      } else {
        rhoguess = 0.5 * (rhomin + rhomax);
      }
    }
    auto status = findRoot(PofRT, press, rhoguess, rhomin, rhomax, robust::EPS(),
                           robust::EPS(), rho);
    // JMM: This needs to not fail and instead return something sane.
    // If root find failed to converge, density will at least be
    // within bounds.
    if (status != Status::SUCCESS) {
      PORTABLE_WARN("DensityEnergyFromPressureTemperature failed to find root\n");
    }
    sie = copy.InternalEnergyFromDensityTemperature(rho, temp, lambda);
    return;
  }
  PORTABLE_INLINE_FUNCTION void DensityEnergyFromPressureTemperature(const Real press,
                                                                     const Real temp,
                                                                     Real &rho,
                                                                     Real &sie) const {
    CRTP copy = *(static_cast<CRTP const *>(this));
    copy.DensityEnergyFromPressureTemperature(press, temp, static_cast<Real *>(nullptr),
                                              rho, sie);
  }

  // Serialization
  /*
    The methodology here is there are *three* size methods all EOS's provide:
    - `SharedMemorySizeInBytes()` which is the amount of memory a class can share
    - `DynamicMemorySizeInBytes()` which is the amount of memory not covered by
    `sizeof(this)`
    - `SerializedSizeInBytes()` which is the total size of the object.

    I wanted serialization machinery to work if you use a standalone
    class or if you use the variant. To make that possible, each class
    provides its own implementation of `SharedMemorySizeInBytes` and
    `DynamicMemorySizeInBytes()`. But then there is a separate
    implementation for the variant and for the base class for
    `SerializedSizeInBytes`, `Serialize`, and `DeSerialize`.
   */

  // JMM: These must frequently be special-cased.
  std::size_t DynamicMemorySizeInBytes() const { return 0; }
  std::size_t DumpDynamicMemory(char *dst) { return 0; }
  std::size_t SetDynamicMemory(char *src,
                               const SharedMemSettings &stngs = DEFAULT_SHMEM_STNGS) {
    return 0;
  }
  // JMM: These usually don't need to be special cased.
  std::size_t SharedMemorySizeInBytes() const {
    const CRTP *pcrtp = static_cast<const CRTP *>(this);
    return pcrtp->DynamicMemorySizeInBytes();
  }
  constexpr bool AllDynamicMemoryIsShareable() const { return true; }
  // JMM: These are generic and never need to be special-cased.
  // However, there must be a separate implementation for these
  // separately in the base class and in the variant.
  std::size_t SerializedSizeInBytes() const {
    // sizeof(*this) returns the size of JUST the base class.
    const CRTP *pcrtp = static_cast<const CRTP *>(this);
    std::size_t dyn_size = pcrtp->DynamicMemorySizeInBytes();
    return dyn_size + sizeof(CRTP);
  }
  std::size_t Serialize(char *dst) {
    CRTP *pcrtp = static_cast<CRTP *>(this);
    memcpy(dst, pcrtp, sizeof(CRTP));
    std::size_t offst = sizeof(CRTP);
    std::size_t dyn_size = pcrtp->DynamicMemorySizeInBytes();
    if (dyn_size > 0) {
      offst += pcrtp->DumpDynamicMemory(dst + sizeof(CRTP));
    }
    const std::size_t tot_size = pcrtp->SerializedSizeInBytes();
    PORTABLE_ALWAYS_REQUIRE(offst == tot_size, "Serialization failed!");
    return offst;
  }
  auto Serialize() {
    CRTP *pcrtp = static_cast<CRTP *>(this);
    std::size_t size = pcrtp->SerializedSizeInBytes();
    char *dst = (char *)malloc(size);
    std::size_t size_new = Serialize(dst);
    PORTABLE_ALWAYS_REQUIRE(size_new == size, "Serialization failed!");
    return std::make_pair(size, dst);
  }
  std::size_t DeSerialize(char *src,
                          const SharedMemSettings &stngs = DEFAULT_SHMEM_STNGS) {
    CRTP *pcrtp = static_cast<CRTP *>(this);
    memcpy(pcrtp, src, sizeof(CRTP));
    std::size_t offst = sizeof(CRTP);
    std::size_t dyn_size = pcrtp->DynamicMemorySizeInBytes();
    if (dyn_size > 0) {
      const bool sizes_same = pcrtp->AllDynamicMemoryIsShareable();
      if (stngs.CopyNeeded() && sizes_same) {
        memcpy(stngs.data, src + offst, dyn_size);
      }
      offst += pcrtp->SetDynamicMemory(src + offst, stngs);
    }
    const std::size_t tot_size = pcrtp->SerializedSizeInBytes();
    PORTABLE_ALWAYS_REQUIRE(offst == tot_size, "Deserialization failed!");
    return offst;
  }

  // Tooling for indexers
  template <typename T>
  static inline constexpr bool NeedsLambda() {
    return false;
  }
  template <typename T>
  static inline constexpr bool NeedsLambda(const T &t) {
    return NeedsLambda<T>();
  }

  // Tooling for modifiers
  static inline constexpr bool IsModified() { return false; }

  inline constexpr decltype(auto) UnmodifyOnce() { return *static_cast<CRTP *>(this); }

  inline constexpr decltype(auto) GetUnmodifiedObject() {
    if constexpr (CRTP::IsModified()) {
      auto unmodified =
          ((static_cast<CRTP *>(this))->UnmodifyOnce()).GetUnmodifiedObject();
      return unmodified;
    } else {
      return *static_cast<CRTP *>(this);
    }
  }
};
} // namespace eos_base
} // namespace singularity

#undef SG_MEMBER_FUNC_NAME
#endif
