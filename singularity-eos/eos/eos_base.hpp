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

#ifndef _SINGULARITY_EOS_EOS_EOS_BASE_
#define _SINGULARITY_EOS_EOS_EOS_BASE_

#include <cstring>
#include <string>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>

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
constexpr std::size_t MAX_NUM_CHARS = 81;
// Cuda doesn't have strcat, so we implement it ourselves
PORTABLE_FORCEINLINE_FUNCTION
char *StrCat(char *destination, const char *source) {
  int i, j; // not in loops because they're re-used.

  // specifically avoid strlen, which isn't on GPU
  for (i = 0; destination[i] != '\0'; i++) {
  }
  // assumes destination has enough memory allocated
  for (j = 0; source[j] != '\0'; j++) {
    // MAX_NUM_CHARS-1 to leave room for null terminator
    PORTABLE_REQUIRE((i + j) < MAX_NUM_CHARS - 1,
                     "Concat string must be within allowed size");
    destination[i + j] = source[j];
  }
  // null terminate destination string
  destination[i + j] = '\0';

  // the destination is returned by standard `strcat()`
  return destination;
}
} // namespace impl

// This Macro adds the `using` statements that allow for the base class
// vector functionality to overload the scalar implementations in the derived
// classes
// TODO(JMM): Should we have more macros that capture just some of these?
#define SG_ADD_BASE_CLASS_USINGS(EOSDERIVED)                                             \
  using EosBase<EOSDERIVED>::TemperatureFromDensityInternalEnergy;                       \
  using EosBase<EOSDERIVED>::InternalEnergyFromDensityTemperature;                       \
  using EosBase<EOSDERIVED>::PressureFromDensityTemperature;                             \
  using EosBase<EOSDERIVED>::PressureFromDensityInternalEnergy;                          \
  using EosBase<EOSDERIVED>::SpecificHeatFromDensityTemperature;                         \
  using EosBase<EOSDERIVED>::SpecificHeatFromDensityInternalEnergy;                      \
  using EosBase<EOSDERIVED>::BulkModulusFromDensityTemperature;                          \
  using EosBase<EOSDERIVED>::BulkModulusFromDensityInternalEnergy;                       \
  using EosBase<EOSDERIVED>::GruneisenParamFromDensityTemperature;                       \
  using EosBase<EOSDERIVED>::GruneisenParamFromDensityInternalEnergy;                    \
  using EosBase<EOSDERIVED>::MinimumDensity;                                             \
  using EosBase<EOSDERIVED>::MinimumTemperature;                                         \
  using EosBase<EOSDERIVED>::PTofRE;                                                     \
  using EosBase<EOSDERIVED>::FillEos;                                                    \
  using EosBase<EOSDERIVED>::EntropyFromDensityTemperature;                              \
  using EosBase<EOSDERIVED>::EntropyFromDensityInternalEnergy;                           \
  using EosBase<EOSDERIVED>::EntropyIsNotEnabled;                                        \
  using EosBase<EOSDERIVED>::IsModified;                                                 \
  using EosBase<EOSDERIVED>::UnmodifyOnce;                                               \
  using EosBase<EOSDERIVED>::GetUnmodifiedObject;

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
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  TemperatureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                       RealIndexer &&temperatures, Real * /*scratch*/,
                                       const int num, LambdaIndexer &&lambdas) const {
    TemperatureFromDensityInternalEnergy(rhos, sies, temperatures, num, lambdas);
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
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void InternalEnergyFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&sies, Real * /*scratch*/,
                                                   const int num,
                                                   LambdaIndexer &&lambdas) const {
    InternalEnergyFromDensityTemperature(rhos, temperatures, sies, num, lambdas);
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
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  PressureFromDensityTemperature(ConstRealIndexer &&rhos, ConstRealIndexer &&temperatures,
                                 RealIndexer &&pressures, Real * /*scratch*/,
                                 const int num, LambdaIndexer &&lambdas) const {
    PressureFromDensityTemperature(rhos, temperatures, pressures, num, lambdas);
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
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  PressureFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                    RealIndexer &&pressures, Real * /*scratch*/,
                                    const int num, LambdaIndexer &&lambdas) const {
    PressureFromDensityInternalEnergy(rhos, sies, pressures, num, lambdas);
  }
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
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  EntropyFromDensityTemperature(ConstRealIndexer &&rhos, ConstRealIndexer &&temperatures,
                                RealIndexer &&entropies, Real * /*scratch*/,
                                const int num, LambdaIndexer &&lambdas) const {
    EntropyFromDensityTemperature(rhos, temperatures, entropies, num, lambdas);
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
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  EntropyFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                   RealIndexer &&entropies, Real * /*scratch*/,
                                   const int num, LambdaIndexer &&lambdas) const {
    EntropyFromDensityInternalEnergy(rhos, sies, entropies, num, lambdas);
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
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void SpecificHeatFromDensityTemperature(ConstRealIndexer &&rhos,
                                                 ConstRealIndexer &&temperatures,
                                                 RealIndexer &&cvs, Real * /*scratch*/,
                                                 const int num,
                                                 LambdaIndexer &&lambdas) const {
    SpecificHeatFromDensityTemperature(rhos, temperatures, cvs, num, lambdas);
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
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  SpecificHeatFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                        RealIndexer &&cvs, Real * /*scratch*/,
                                        const int num, LambdaIndexer &&lambdas) const {
    SpecificHeatFromDensityInternalEnergy(rhos, sies, cvs, num, lambdas);
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
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void BulkModulusFromDensityTemperature(ConstRealIndexer &&rhos,
                                                ConstRealIndexer &&temperatures,
                                                RealIndexer &&bmods, Real * /*scratch*/,
                                                const int num,
                                                LambdaIndexer &&lambdas) const {
    BulkModulusFromDensityTemperature(rhos, temperatures, bmods, num, lambdas);
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
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void
  BulkModulusFromDensityInternalEnergy(ConstRealIndexer &&rhos, ConstRealIndexer &&sies,
                                       RealIndexer &&bmods, Real * /*scratch*/,
                                       const int num, LambdaIndexer &&lambdas) const {
    BulkModulusFromDensityInternalEnergy(rhos, sies, bmods, num, lambdas);
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
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void GruneisenParamFromDensityTemperature(ConstRealIndexer &&rhos,
                                                   ConstRealIndexer &&temperatures,
                                                   RealIndexer &&gm1s, Real * /*scratch*/,
                                                   const int num,
                                                   LambdaIndexer &&lambdas) const {
    GruneisenParamFromDensityTemperature(rhos, temperatures, gm1s, num, lambdas);
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
  template <typename RealIndexer, typename ConstRealIndexer, typename LambdaIndexer>
  inline void GruneisenParamFromDensityInternalEnergy(ConstRealIndexer &&rhos,
                                                      ConstRealIndexer &&sies,
                                                      RealIndexer &&gm1s,
                                                      Real * /*scratch*/, const int num,
                                                      LambdaIndexer &&lambdas) const {
    GruneisenParamFromDensityInternalEnergy(rhos, sies, gm1s, num, lambdas);
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

  template <typename RealIndexer, typename LambdaIndexer>
  inline void PTofRE(RealIndexer &&rhos, RealIndexer &&sies, RealIndexer &&presses,
                     RealIndexer &&temps, RealIndexer &&dpdrs, RealIndexer &&dpdes,
                     RealIndexer &&dtdrs, RealIndexer &&dtdes, const int num,
                     LambdaIndexer &&lambdas) const {
    static auto const name = SG_MEMBER_FUNC_NAME();
    static auto const cname = name.c_str();
    CRTP copy = *(static_cast<CRTP const *>(this));
    portableFor(
        cname, 0, num, PORTABLE_LAMBDA(const int i) {
          copy.PTofRE(rhos[i], sies[i], lambdas[i], presses[i], temps[i], dpdrs[i],
                      dpdes[i], dtdrs[i], dtdes[i]);
        });
  }
  // Scalar version of PTofRE
  PORTABLE_INLINE_FUNCTION
  void PTofRE(Real &rho, Real &sie, Real *lambda, Real &press, Real &temp, Real &dpdr,
              Real &dpde, Real &dtdr, Real &dtde) const {
    // Get the dervived class
    auto eos = static_cast<CRTP const &>(*this);

    press = eos.PressureFromDensityInternalEnergy(rho, sie, lambda);
    temp = eos.TemperatureFromDensityInternalEnergy(rho, sie, lambda);
    const Real drho = rho * 1.0e-6;
    const Real de = sie * 1.0e-6;
    const Real Pr = eos.PressureFromDensityInternalEnergy(rho + drho, sie, lambda);
    const Real Pe = eos.PressureFromDensityInternalEnergy(rho, sie + de, lambda);
    const Real Tr = eos.TemperatureFromDensityInternalEnergy(rho + drho, sie, lambda);
    const Real Te = eos.TemperatureFromDensityInternalEnergy(rho, sie + de, lambda);
    dpdr = (Pr - press) / drho;
    dpde = (Pe - press) / de;
    dtdr = (Tr - temp) / drho;
    dtde = (Te - temp) /
           de; // Would it be better to skip the calculation of Te and return 1/cv?
    return;
  }

  PORTABLE_INLINE_FUNCTION
  Real RhoPmin(const Real temp) const { return 0.0; }

  // Default entropy behavior is to return an error
  PORTABLE_FORCEINLINE_FUNCTION
  void EntropyIsNotEnabled(const char *eosname) const {
    // Construct the error message using char* so it works on device
    // WARNING: This needs to be updated if EOS names get longer
    // base msg length 32 + 5 chars = 37 chars
    // + 1 char for null terminator
    // maximum allowed EOS length = 44 chars
    char msg[impl::MAX_NUM_CHARS] = "Entropy is not enabled for the '";
    impl::StrCat(msg, eosname);
    impl::StrCat(msg, "' EOS");
    PORTABLE_ALWAYS_THROW_OR_ABORT(msg);
  }

  // Tooling for modifiers
  inline constexpr bool IsModified() const { return false; }

  inline constexpr decltype(auto) UnmodifyOnce() { return *static_cast<CRTP *>(this); }

  inline constexpr decltype(auto) GetUnmodifiedObject() { return *static_cast<CRTP *>(this); }
};
} // namespace eos_base
} // namespace singularity

#undef SG_MEMBER_FUNC_NAME
#endif
