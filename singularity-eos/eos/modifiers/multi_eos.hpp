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

#ifndef _SINGULARITY_EOS_EOS_MULTI_EOS_
#define _SINGULARITY_EOS_EOS_MULTI_EOS_

#include "stdio.h"
#include <cstdlib>
#include <iostream>
#include <tuple>
#include <algorithm>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/tuple_utils.hpp>
#include <singularity-eos/eos/eos_base.hpp>

// SFINAE helper macro to make sure the lambda has mass fraction information for
// all EOS in the MultiEOS
// TODO: Does this trigger an issue where all EOS in the Variant must be passed
// mass fraction information? Do we want this behavior? Or do we want a runtime
// error when this EOS is in the variant but _isn't_ passed mass fraction
// information?
#define LAMBDA_HAS_MASS_FRACTION_INDICES                                                 \
  typename std::enable_if<                                                               \
        variadic_utils::is_indexable_v<                                                  \
            Indexer_t, IndexableTypes::MassFraction<nmat_>                               \
        >                                                                                \
      >::type

namespace singularity {

using namespace eos_base;

/*
This is a modifer that takes three closure model template parameters with an
arbitrary number of EOS and combines them so that they can essentially be
grouped as a single material. Mass fractions are then provided through the
lambda parameter and a closure rule determines the combined material response.

Note that three closure models are needed because the the different known
mixture quantities, i.e. pressure, density, temperature, and internal energy,
can be combined in different ways to specify the system. The PTE solver for
each case _must_ be different unless an iterative wrapper is employed. Each
closure model must mirror the API for singularity-eos closure models, although
it may be preferable to write analytic closures (or analytic matrix inversions)
in certain cases.

From a design perspective, this code involves a LOT of fold expressions to
expand "loops" over the EOS at compile time.
*/

template<
      typename ClosureRE
    , typename ClosureRT
    , typename ClosureRP
    , typename... EOSModelsT
>
class MultiEOS : public EosBase<MultiEOS<ClosureRE, ClosureRT, ClosureRP, EOSModelsT...>> {
 public:
  SG_ADD_BASE_CLASS_USINGS(MultiEOS)

  static std::string EosType() {
    std::string all_models;
    // Use a fold to evaluate the expression in parentheses with all elements of
    // the pack (C++17)
    ((all_models += (all_models.empty() ? "" : ", ") + EOSModelsT::EosType()), ...);
    return "MultiEOS<"ClosureRE::MethodType() + ", " + ClosureRT::MethodType() + ", "
      + ClosureRP::MethodType() + ", " + all_models + ">";
  }

  static std::string EosPyType() {
    // Use a left fold expression to sum the model names (C++17)
    std::string all_models = (std::string() + ... + EOSModelsT::EosType());
    return "MultiEOS" + ClosureRE::MethodType() + ClosureRT::MethodType()
      + ClosureRP::MethodType() + all_models;
  }

  MultiEOS() = default;

  template<typename... EOSModelsT_>
  MultiEOS(EOSModelsT_... eos_models)
      : models_(std::forward<EOSModelsT_>(eos_models))
  {
    CheckParams();
  }

  // Overload that makes copying the models easier (e.g. GetOnDevice()). Make
  // explicit to avoid accidentally storing a tuple of a tuple
  explicit MultiEOS(std::tuple<EOSModelsT...> model_tuple)
      : models_(std::forward<std::tuple<EOSModelsT...>>(model_tuple))
  {
    CheckParams();
  }

  PORTABLE_INLINE_FUNCTION void CheckParams() const {
    // Again we're using a fold expression now with std::apply to call
    // CheckParams() on each model (C++17)
    std::apply([](const auto&... eos_models) {
      (eos_models.CheckParams(), ...);
    }, models_);
  }

  auto GetOnDevice() {
    // Because GetOnDevice() returns a copy of the EOS, we can't really just
    // use std::apply, and instead we need to also populate a new copy of the
    // model tuple
    auto models_on_device = tuple_transform(
      models_, [](const auto& eos) { return eos.GetOnDevice(); }
    );
    return MultiEOS<ClosureRE, ClosureRT, ClosureRP, EOSModelsT...>(models_on_device);
  }

  inline void Finalize() {
    // All the fold expressions!
    std::apply([](const auto&... eos_models) {
      (eos_models.Finalize(), ...);
    }, models_);
  }

  /* Density-energy inputs

  PTESolverRhoT(const std::size_t nmat, EOS_t &&eos, const Real vfrac_tot,
                const Real sie_tot, Real_t &&rho, Real_t &&vfrac, Real_t &&sie,
                Real_t &&temp, Real_t &&press, Lambda_t &&lambda, Real *scratch,
                const Real Tnorm = 0.0, const MixParams &params = MixParams())

  PTESolverPT(const std::size_t nmat, EOS_t &&eos, const Real vfrac_tot,
              const Real sie_tot, Real_t &&rho, Real_t &&vfrac, Real_t &&sie,
              Real_t &&temp, Real_t &&press, Lambda_t &&lambda, Real *scratch,
              const Real Tnorm = 0.0, const MixParams &params = MixParams())

  */

  /* Density-temperature inputs

  PTESolverFixedT(const std::size_t nmat, EOS_t &&eos, const Real vfrac_tot,
                  const Real T_true, Real_t &&rho, Real_t &&vfrac, Real_t &&sie,
                  CReal_t &&temp, Real_t &&press, Lambda_t &lambda, Real *scratch,
                  const MixParams &params = MixParams())

  */

  /* Density-pressure inputs

  PTESolverFixedP(const std::size_t nmat, EOS_t &&eos, const Real vfrac_tot,
                  const Real P, Real_t &&rho, Real_t &&vfrac, Real_t &&sie,
                  Real_t &&temp, CReal_t &&press, Lambda_t &lambda, Real *scratch,
                  const MixParams &params = MixParams())

  */

  template<typename RealIndexer, LambdaIndexer>
  inline SolverStatus SolvePTEFromDensityEnergy(
      const Real sie_tot, const Real density_tot, Real pressure, Real temperature,
      RealIndexer &&density_mat, RealIndexer &&sie_mat, LambdaIndexer &&lambdas,
      Real *scratch
  ) {

  }

  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    // Call PTE solver when density and energy are known
  }
  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    // Call PTE solver when density and energy are known
  }
  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_FUNCTION Real
  MinInternalEnergyFromDensity(const Real rho, Indexer_t &&lambda = nullptr) const {
    // Return mass fraction weighted minimum energy
  }
  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    // Call PTE solver when density and energy are known
  }
  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    // Call PTE solver when density and energy are known
  }
  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    // Call PTE solver when density and energy are known
  }
  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, Indexer_t &&lambda = nullptr) const {
    // Call PTE solver when density and energy are known
  }
  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    // Call PTE solver when density and temperature are known (energy is unknown)
  }
  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    // Call PTE solver when density and temperature are known (energy is unknown)
  }
  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_FUNCTION Real EntropyFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    // Call PTE solver when density and temperature are known (energy is unknown)
  }
  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    // Call PTE solver when density and temperature are known (energy is unknown)
  }
  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    // Call PTE solver when density and temperature are known (energy is unknown)
  }
  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, Indexer_t &&lambda = nullptr) const {
    // Call PTE solver when density and temperature are known (energy is unknown)
  }
  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                 Real &cv, Real &bmod, const unsigned long output,
                                 Indexer_t &&lambda = nullptr) const {
    // Depending on the independent variables, call a different PTE solver and
    // then fill then properly average the rest of the derivative quantities
  }
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept {
    // Add lambdas from EOS and add mass fraction lambda
  }

  static constexpr unsigned long PreferredInput() {
    // Allow EOS to vote on preferred input
  }

  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    // Add up scratch requirements from all EOS for this method as well as the
    // appropriate PTE solver requirements
    unsigned long scratch_size = 0;

    // Find the maximum scratch required for the PTE solvers. Note the "fall
    // through for cases that have the same independent variables"
    switch (method) {
    case "PressureFromDensityTemperature":
      [[fallthrough]];
    case "InternalEnergyFromDensityTemperature":
      [[fallthrough]];
    case "SpecificHeatFromDensityTemperature":
      [[fallthrough]];
    case "BulkModulusFromDensityTemperature":
      [[fallthrough]];
    case "GruneisenParamFromDensityTemperature":
      scratch_size += ClosureRT.RequiredScratchInBytes(nmat_);
      break;
    case "TemperatureFromDensityInternalEnergy":
      [[fallthrough]];
    case "PressureFromDensityInternalEnergy":
      [[fallthrough]];
    case "SpecificHeatFromDensityInternalEnergy":
      [[fallthrough]];
    case "BulkModulusFromDensityInternalEnergy":
      [[fallthrough]];
    case "GruneisenParamFromDensityInternalEnergy":
      scratch_size += ClosureRE.RequiredScratchInBytes(nmat_);
      break;
    case "MinInternalEnergyFromDensity":
      // No scratch needed
      break;
    case "FillEos":
      // Without knowing the exact inputs, we have to assume a worst-case
      // scenario
      scratch_size += std::max(ClosureRT.RequiredScratchInBytes(nmat_),
                               std::max(ClosureRE.RequiredScratchInBytes(nmat_),
                                        ClsoureRP.RequiredScratchInBytes(nmat_))
      );
      break;
    default:
      PORTABLE_ALWAYS_ABORT("No lookup method matches the requested method");
    }

    // Since execution could be asynchronous, require that all scratch for PTE
    // be available for all states simultaneously. This mirrors the EOS scratch
    // size behavior for e.g. EOSPAC
    scratch_size *= nelements;

    // Add in the requirements of the individual EOS. Note the `nelements`
    // argument
    std::apply([](const auto&... eos_models) {
      (scratch_size + ... + eos_models.scratch_size(method, nelements));
    }, models_);

    // Contribution from lambdas... does this change the meaning of
    // MAX_NUM_LAMBDAS?
    scratch_size += nmat_ * MAX_NUM_LAMBDAS * nelements;

    return scratch_size;
  }

  static inline unsigned long max_scratch_size(unsigned int nelements) {
    // Worst case scenario for PTE solve
    unsigned long scratch_size =
        std::max(ClosureRT.RequiredScratchInBytes(nmat_),
                 std::max(ClosureRE.RequiredScratchInBytes(nmat_),
                          ClsoureRP.RequiredScratchInBytes(nmat_))
      );

    // Since execution could be asynchronous, require that all scratch for PTE
    // be available for all states simultaneously. This mirrors the EOS scratch
    // size behavior for e.g. EOSPAC
    scratch_size *= nelements;

    // Add in the maximum requirements of the individual EOS. Note the
    // `nelements` argument
    std::apply([](const auto&... eos_models) {
      (scratch_size + ... + eos_models.max_scratch_size(nelements));
    }, models_);

    // Contribution from lambdas... does this change the meaning of
    // MAX_NUM_LAMBDAS?
    scratch_size += nmat_ * MAX_NUM_LAMBDAS * nelements;

    return scratch_size;

  }

  PORTABLE_FUNCTION void PrintParams() const {
    // print the parameters for each EOS via fold expression
    std::apply([](const auto&... eos_models) {
      (eos_models.PrintParams(), ...);
    }, models_);
  }

  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const {
    // Average the contributions from each EOS
  }

  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_FUNCTION void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie,
                                                Real &press, Real &cv, Real &bmod,
                                                Real &dpde, Real &dvdt,
                                                Indexer_t &&lambda = nullptr) const {
    // Average the contributions from each EOS
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumDensity() const {
    // Max over the contributions from each EOS
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumTemperature() const {
    // Max over the contributions from each EOS
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MaximumDensity() const {
    // Min over the contributions from each EOS
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumPressure() const {
    // Max over the contributions from each EOS
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MaximumPressureAtTemperature(const Real temp) const {
    // Min over the contributions from each EOS
  }

  // 3T member functions
  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicMassFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    // Average over the contributions from each EOS
  }
  template <typename Indexer_t = Real *, LAMBDA_HAS_MASS_FRACTION_INDICES>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicNumberFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    // Average over the contributions from each EOS
  }
  PORTABLE_INLINE_FUNCTION
  Real MeanAtomicMass() const {
    // Average over the contributions from each EOS
  }
  PORTABLE_INLINE_FUNCTION
  Real MeanAtomicNumber() const {
    // Average over the contributions from each EOS
  }

  // Modifier member functions
  static inline constexpr bool IsModified() {
    return true;
  }
  inline constexpr T UnmodifyOnce() {
    // Return the tuple of EOS objects maybe?
  }

  // Dynamic and shared memory member functions (normally included via macro)
  std::size_t DynamicMemorySizeInBytes() const {
    // Add the contributions from each EOS
  }
  std::size_t SharedMemorySizeInBytes() const {
    // Add the contributions from each EOS
  }
  std::size_t DumpDynamicMemory(char *dst) {
    // It seems like we'll need to dump from each EOS and then advance the
    // pointer to the next stretch of available memory
  }
  std::size_t SetDynamicMemory(char *src,
                               const SharedMemSettings &stngs = DEFAULT_SHMEM_STNGS) {
    // I think I should loop through the source memory and provide each EOS with
    // a pointer to where its memory begins
  }
  constexpr bool AllDynamicMemoryIsShareable() const {
    // If any is false, then all are false
  }


 private:
  std::tuple<EOSModelsT...> models_;
  static constexpr std::size_t nmat_ = sizeof...(EOSModelsT);  // for convenience

};

// TODO: Template deduction guide to allow providing the EOS models without
// providing the corresponding template parameters



#endif // #ifndef _SINGULARITY_EOS_EOS_MULTI_EOS_
