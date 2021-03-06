//------------------------------------------------------------------------------
// © 2021-2022. Triad National Security, LLC. All rights reserved.  This
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

#include <string>

#include <ports-of-call/portability.hpp>

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
  using EosBase<EOSDERIVED>::FillEos;

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

  // Specialzied vector version of PTofRE maybe more suited for EOSPAC
  // }
  // template<typename RealIndexer, typename LambdaIndexer>
  // inline
  // void PTofRE(RealIndexer &&rhos, RealIndexer &&sies,
  //             RealIndexer &&presses, RealIndexer &&temps,
  //             RealIndexer &&dpdrs, RealIndexer &&dpdes,
  //             RealIndexer &&dtdrs, RealIndexer &&dtdes,
  //             const int num, LambdaIndexer &&lambdas) const {
  //   // Get the dervived class
  //   auto eos = static_cast<CRTP const&>(*this);

  //   // Lookup at density and energy
  //   eos.PressureFromDensityInternalEnergy(rhos, sies, num, presses,
  //                                               lambdas);
  //   eos.TemperatureFromDensityInternalEnergy(rhos, sies, num, temps,
  //                                                  lambdas);
  //   // Peturbation factor
  //   Real factor = 1. + 1.0e-06;

  //   // Perturb densities and do lookups
  //   portableFor(
  //       'PerturbDensities', 0, num, PORTABLE_LAMBDA(const int i) {
  //         rhos[i] *= factor;
  //       }
  //   );
  //   eos.PressureFromDensityInternalEnergy(rhos, sies, num, dpdrs,
  //                                               lambdas);
  //   eos.TemperatureFromDensityInternalEnergy(rhos, sies, num, dtdrs,
  //                                                  lambdas);

  //   // Reset densities, perturb energies, and do lookups
  //   portableFor(
  //       'PerturbEnergiesResetDensities', 0, num,
  //       PORTABLE_LAMBDA(const int i) {
  //         sies[i] *= factor;
  //         rhos[i] /= factor;
  //       }
  //   );
  //   eos.PressureFromDensityInternalEnergy(rhos, sies, dpdes, num,
  //                                               lambdas);
  //   eos.TemperatureFromDensityInternalEnergy(rhos, sies, dtdes, num,
  //                                                  lambdas);

  //   // Reset the energies to their original values
  //   portableFor(
  //       'ResetEnergies', 0, num, PORTABLE_LAMBDA(const int i) {
  //         sies[i] /= factor;
  //       }
  //   );

  //   // Calculate the derivatives
  //   portableFor(
  //       'CalculateDerivatives', 0., num, PORTABLE_LAMBDA(const int i) {
  //         dpdrs[i] = (dpdrs[i] - presses[i]) / (rhos[i] * (1. - factor));
  //         dpdes[i] = (dpdes[i] - presses[i]) / (sies[i] * (1. - factor));
  //         dtdrs[i] = (dtdrs[i] - temps[i]) / (rhos[i] * (1. - factor));
  //         dtdes[i] = (dtdes[i] - temps[i]) / (sies[i] * (1. - factor));
  //       }
  //   );

  //   return;
  // }
};
} // namespace eos_base
} // namespace singularity

#undef SG_MEMBER_FUNC_NAME
#endif
