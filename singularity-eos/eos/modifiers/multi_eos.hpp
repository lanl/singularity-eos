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
#include <utility>
#include <algorithm>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/base/tuple_utils.hpp>
#include <singularity-eos/base/indexable_types.hpp>
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>
#include <singularity-eos/eos/eos_base.hpp>
#include <singularity-eos/eos/eos_variant.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>

// SFINAE helper macro to make sure the lambda has mass fraction information for
// all EOS in the MultiEOS
// TODO: Does this trigger an issue where all EOS in the Variant must be passed
// mass fraction information? Do we want this behavior? Or do we want a runtime
// error when this EOS is in the variant but _isn't_ passed mass fraction
// information?
#define LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)                                  \
  typename std::enable_if<                                                               \
        variadic_utils::is_indexable_v<                                                  \
            LambdaIndexer, IndexableTypes::MassFraction<nmat_>                           \
        >                                                                                \
      >::type

namespace singularity {

// This functor essentially wraps a fold expression so that each term in the
// fold can be injected via the functor f. The functor should take the EOS
// followed by the density, internal energy, and temperature so that it can
// perform either density-energy or density-temperature lookups with the EOS.
// The functor should also return the _inverse_ term so that a harmonic average
// is properly formed. Note that similar averaging functors would need to obey
// this basic API.
struct VolumeFracHarmonicAverageFunctor {
  template<typename Func, typename EOSmodelsT, typename RealIndexer,
           typename LambdaIndexer, std::size_t... Is>
  constexpr Real operator()(Func&& f, EOSmodelsT models, RealIndexer density_arr,
                            RealIndexer sie_arr, Real density, Real temperature,
                            LambdaIndexer const& lambda,
                            std::index_sequence<Is...>) const {
    // Note that the function should appropriately invert any quantity that
    return robust_ratio(
      1.0,
      (
        ( /* volume fraction */
          lambda[IndexableTypes::MassFraction<Is>{}] * density / density_arr[Is]
          /* Inverse value to be volume-avreaged */
           * robust_ratio(
              1.0,
              std::invoke(
                std::forward<Func>(f),
                std::get<Is>(models),
                density_arr[Is],
                sie_arr[Is],
                temperature,
                lambda
              )
           )
        )
      + ... )
    );
  }
};

struct MassFracAverageFunctor {
  template<typename Func, typename EOSmodelsT, typename RealIndexer,
           typename LambdaIndexer, std::size_t... Is>
  constexpr Real operator()(Func&& f, EOSmodelsT models, RealIndexer density_arr,
                            RealIndexer sie_arr, [[maybe_unused]] Real density,
                            Real temperature, LambdaIndexer const& lambda,
                            std::index_sequence<Is...>) const {
    return
      (
        (
          lambda[IndexableTypes::MassFraction<Is>{}]
          /* Value to be volume-avreaged */
          * std::invoke(
             std::forward<Func>(f),
             std::get<Is>(models),
             density_arr[Is],
             sie_arr[Is],
             temperature,
             lambda
          )
        )
      + ... );
  }
};

// Index into an arbitrary array with an arbitrary map
template<typename arrT, typename mapT>
struct GenericIndexer {
  arrT arr_;
  mapT map_;

  template<typename arrT_, typename mapT_>
  constexpr
  GenericIndexer(arrT_&& arr_in, mapT_&& map_in)
      : arr_{std::forward<arrT_>(arr_in)}
      , map_{std::forward<mapT_>(map_in)}
  {}

  template<typename idxT>
  constexpr
  auto &operator[](const idxT i) { return arr_[map_[i]]; }

  template<typename idxT>
  constexpr
  const auto &operator[](const idxT i) const { return arr_[map_[i]]; }
};
// CTAD
template<typename arrT_, typename mapT_>
GenericIndexer(arrT_&& arr_in, mapT_&& map_in) -> GenericIndexer<arrT_, mapT_>;

using namespace eos_base;

/*
This is a modifer that assumes three closure models with an arbitrary number of
EOS and combines them so that they can essentially be grouped as a single
material. Mass fractions are then provided through the lambda parameter and a
closure rule determines the combined material response.

Note that three closure models are needed because the the different known
mixture quantities, i.e. pressure, density, temperature, and internal energy,
can be combined in different ways to specify the system. The PTE solver for
each case _must_ be different unless an iterative wrapper is employed.

From a design perspective, this code involves a LOT of fold expressions to
expand "loops" over the EOS at compile time. Note that anything static can use
a fold expression by itself, but non-static member functions need std::apply to
unpack the tuple of models.

TODO: It would be nice to be able to inject the PTE solvers rather than hard
coding them. Right now this is difficult because the entire class is templated
on the PTE state arrays, so specifying the type is extremely difficult to do
when _this_ class would be instantiated. Ideally, the types passed to the PTE
solver would only matter when the PTE solver is called rather than when it's
instantiated. This would allow injection of the PTE solver methods, allowing
analytic PTE solvers where possible.
*/

template<
      typename BulkModAvgT = VolumeFracHarmonicAverageFunctor
    , typename GrunesienAvgT = VolumeFracHarmonicAverageFunctor
    , typename... EOSModelsT
>
class MultiEOS : public EosBase<MultiEOS<BulkModAvgT, GruneisenAvgT, EOSModelsT...>> {
 private:
  Real mass_frac_cutoff_;
  std::tuple<EOSModelsT...> models_;
  using eosT_ = Variant<EOSModelsT*...>;

  static constexpr std::size_t nmat_ = sizeof...(EOSModelsT);  // for convenience

  // NOTE: These functions are private since they're details of the
  // implementation

  // Implementation function to unpack the lambda mass fraction indices into an
  // array. I believe this has to be separate from the public implementation to
  // get the 'Is' template parameter to specify the index sequence
  template<typename RealIndexer, typename LambdaIndexer, size_t... Is,
           LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)
  >
  constexpr void assign_mass_fractions(RealIndexer mass_fracs, const LambdaIndexer lambda,
                                       std::index_sequence<Is...>) {
    ( (mass_fracs[Is] = lambda[IndexableTypes::MassFraction<Is>{}]), ...);
  }

  // Implementation function to call DensityEnergyFromPressureTemperature() on
  // each EOS for a common pressure and temperature and populate material arrays
  template<typename RealIndexer, typename LambdaIndexer, size_t... Is,
           LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)
  >
  constexpr void callDensityEnergyFromPressureTemperatureAll(
      const Real pressure, const Real temperature, LambdaIndexer&& lambdas, RealIndexer& density_mat,
      RealIndexer& energy_mat, std::index_sequence<Is...>) {
    ( (std::get<Is>(models_).DensityEnergyFromPressureTemperature(
          pressure, temperature, std::forward<LambdaIndexer>(lambdas), density_mat[Is],
          energy_mat[Is])), ...
    );
  }

  // Implementation of mass fraction average over all EOS but at a single state
  // given by input1 and input2
  template<typename Func, typename LambdaIndexer,
           std::size_t... Is>
  Real massFracAvgQuantityAtOneState(Func&& f, Real input1, Real input2,
                                  LambdaIndexer const& lambda,
                                  std::index_sequence<Is...>) const {
    // Note that the function signature should take the EOS as its first argument
    // and then the two independent variables and the lambda
    return (  ( lambda[IndexableTypes::MassFraction<Is>{}]
                 * std::invoke(std::forward<Func>(f),
                               std::get<Is>(models_),
                               input1,
                               input2,
                               lambda)
              )
            + ... );
  }

  // TODO: Is this wrapper necessary? Just use bare functor...
  // Implementation of mass fraction average over all EOS, each at its own state
  // Both internal energy and temperature are provided along with density so that
  // the function wrapper can include logic to decide which lookup to use
  template<typename Func, typename RealIndexer, typename LambdaIndexer,
           std::size_t... Is>
  Real massFracAverageQuantityAtManyStates(Func&& f, RealIndexer density_arr, RealIndexer sie_arr,
                                    Real temperature,
                                    LambdaIndexer const& lambda) const {
    MassFracAverageFunctor avg_funct{};
    // Note that the f function signature should take the EOS as its first argument
    // followed by the lookup state and lambda
    return avg_funct(f, models_, density_arr, sie_arr, temperature, /* density */ 1.0,
                    density_arr, sie_arr, lambda, make_index_sequence<nmat_>)
  }

  // Implementation function to call BulkModulusFromDensityTemperature()
  // on each EOS with a common temperature and array of densities
  template<typename RealIndexer, typename LambdaIndexer, size_t... Is,
           LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)
  >
  constexpr void getMaterialBulkModulusFromDensityTemperature(
      const Real temperature, LambdaIndexer&& lambdas, RealIndexer& density_mat,
      RealIndexer& bmod_mat, std::index_sequence<Is...>) {
    ( (bmod_mat[Is] = std::get<Is>(models_).BulkModulusFromDensityTemperature(
          density_mat[Is], temperature, std::forward<LambdaIndexer>(lambdas))), ...
    );
  }

  // Implementation function to call BulkModulusFromDensityInternalEnergy()
  // on each EOS with arrays of densities and internal energies
  template<typename RealIndexer, typename LambdaIndexer, size_t... Is,
           LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)
  >
  constexpr void getMaterialBulkModulusFromDensityInternalEnergy(
      LambdaIndexer&& lambdas, RealIndexer& density_mat, RealIndexer& sie_mat,
      RealIndexer& bmod_mat, std::index_sequence<Is...>) {
    ( (bmod_mat[Is] = std::get<Is>(models_).BulkModulusFromDensityInternalEnergy(
          density_mat[Is], sie_mat[Is], std::forward<LambdaIndexer>(lambdas))), ...
    );
  }

 public:
  SG_ADD_BASE_CLASS_USINGS(MultiEOS)

  static std::string EosType() {
    std::string all_models;
    // Use a fold to evaluate the expression in parentheses with all elements of
    // the pack (C++17)
    ((all_models += (all_models.empty() ? "" : ", ") + EOSModelsT::EosType()), ...);
    return "MultiEOS<" + all_models + ">";
  }

  static std::string EosPyType() {
    // Use a left fold expression to sum the model names (C++17)
    std::string all_models = (std::string() + ... + EOSModelsT::EosType());
    return "MultiEOS" + all_models;
  }

  MultiEOS() = default;

  template<typename... EOSModelsT_>
  MultiEOS(Real mass_frac_cutoff_in, EOSModelsT_&&... eos_models)
      : models_(std::forward<EOSModelsT_>(eos_models)...)
      , mass_frac_cutoff_{std::forward<Real>(mass_frac_cutoff_in)}
  {
    CheckParams();
  }

  template<typename... EOSModelsT_>
  MultiEOS(EOSModelsT_&&... eos_models)
    : MultiEOS(
        /*mass_frac_cutoff=*/ 1e-8,
        std::forward<EOSModelsT_>(eos_models)...
      )
  {}

  // Overload that makes copying the models easier (e.g. GetOnDevice()). Make
  // explicit to avoid accidentally storing a tuple of a tuple
  explicit MultiEOS(Real mass_frac_cutoff_in, std::tuple<EOSModelsT...> model_tuple)
      : models_(std::move(model_tuple))
      , mass_frac_cutoff_{mass_frac_cutoff_in}
  {
    CheckParams();
  }

  PORTABLE_INLINE_FUNCTION void CheckParams() const {
    // Again we're using a fold expression now with std::apply to call
    // CheckParams() on each model (C++17)
    PORTABLE_ALWAYS_REQUIRE(mass_frac_cutoff_ > 0., "Mass fracton cutoff must be non-negative");
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
    return MultiEOS<BulkModAvgT, EOSModelsT...>(mass_frac_cutoff_, models_on_device);
  }

  inline void Finalize() {
    // All the fold expressions!
    std::apply([](const auto&... eos_models) {
      (eos_models.Finalize(), ...);
    }, models_);
  }

  // Compute volume fractions and possibily material densities from mass
  // fractions and total density
  template<typename MassFracT, typename VolFracT, typename DensityT, typename SieT,
           typename PTEMatT>
  PORTABLE_FORCEINLINE_FUNCTION
  void SetUpPTE(MassFracT mass_fracs, VolFracT vol_fracs,
                DensityT density_mat, SieT sie_mat, PTEMatT pte_mats, size_t &npte,
                Real &vfrac_tot, const Real density_tot, const Real sie_tot) const {
    npte = 0;
    vfrac_tot = 0;
    for (size_t m = 0; m < nmat_; m++) {

      // Check whether a material is present or not
      if (mass_fracs[m] <= mass_frac_cutoff_) {
        // Make sure sane values are returned for material density/energy even
        // if material isn't participating
        vol_fracs[m] = 0.;
        density_mat[m] = density_tot;
        sie_mat[m] = sie_tot;
        continue;
      }

      // Material will participate in PTE
      pte_mats[npte] = m;
      npte += 1;

      // Calculate volume fraction from input density to decide whether it was
      // sane or not. The volume fraction lower bound is essentially the volume
      // fraction of a material at the mass fraction limit whose density is
      // `density_factor_max` more dense than the average material density.
      const auto vol_frac_test = robust::ratio(density_tot, density_mat[m])
          / mass_fracs[m];
      constexpr auto density_factor_max = 1.0e6;
      const auto vol_frac_cutoff = mass_frac_cutoff_ / density_factor_max;
      if (!std::isnormal(density_mat[m]) || vol_frac_test > 1.0
          || vol_frac_test < vol_frac_cutoff) {
        // Bad density guess... fall back to average value
        density_mat[m] = density_tot;
        sie_mat[m] = sie_tot;
        vol_fracs[m] = mass_fracs[m]; // Consistent with density_mat[m] = density_tot
      } else {
        vol_fracs[m] = vol_frac_test;
      }

      // This ensures volume fractions are consistent with material densities
      vfrac_tot += vol_fracs[m];
    }

    // JHP: I think we should always test this unless it is proven to be a performance
    // bottleneck
    PORTABLE_ALWAYS_REQUIRE(npte > 0, "No material mass fraction is greater than cutoff");
  }

  // Helper function to volume-average an array
  template<typename arrT, typename vfracT>
  PORTABLE_FORCEINLINE_FUNCTION
  Real VolumeAverageMatArray(arrT arr, vfracT vfracs, int num = nmat_) {
    Real avg = 0;
    for (size_t i=0; i < num; i++) {
      avg += vfracs[i] * arr[i];
    }
    return avg;
  }

  // Wrappers for the various types of independent variables for the mixture
  template<typename RealIndexer, typename LambdaIndexer,
           LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_INLINE_FUNCTION
  SolverStatus GetStatesFromDensityEnergy(
      const Real density_tot, const Real sie_tot, Real &pressure, Real &temperature,
      RealIndexer &&density_mat, RealIndexer &&sie_mat, LambdaIndexer &&lambdas
  ) const {
    // Create temporary arrays
    std::array<Real, nmat_> mass_fracs{};
    std::array<Real, nmat_> vol_fracs{};
    std::array<size_t, nmat_> pte_mats{};
    std::array<Real, nmat_> pres_mat{};
    std::array<Real, nmat_> temp_mat{};

    // Extract the mass fractions from the lambdas
    PopulateMassFracArray(mass_fracs, lambdas); // Guarantees sum = 1

    // Compute volume fractions from mass fractions and densities. Mark
    // materials to participate in PTE
    size_t npte = 0;
    Real vfrac_tot = 0;
    SetUpPTE(mass_fracs, vol_fracs, density_mat, sie_mat, pte_mats, npte, vfrac_tot,
             density_tot, sie_tot);

    // Initialize solver scratch memory (including cache)
    constexpr int pte_solver_scratch_size = PTESolverPTRequiredScratch(nmat_) + nmat_
        * MAX_NUM_LAMBDAS;
    std::array<Real, pte_solver_scratch_size> solver_scratch{};

    // Get cache from offsets into scratch
    const int neq = npte + 1;
    singularity::mix_impl::CacheAccessor cache(solver_scratch.data() +
                                               neq * (neq + 4) + 2 * npte);

    // Create indexers for the EOS
    const auto eos_idxr = GenericIndexer(CreateEOSArray(), pte_mats);

    // Solve for PTE state
    SolverStatus status;
    if (npte > 1) {
      // Create indexers for the various arrays. Note that we have to use
      // references since the PTE solver expects them all to be the same type
      auto density_idxr = GenericIndexer(&density_mat[0], pte_mats);
      auto vfrac_idxr = GenericIndexer(&vol_fracs[0], pte_mats);
      auto sie_idxr = GenericIndexer(&sie_mat[0], pte_mats);
      auto temp_idxr = GenericIndexer(&temp_mat[0], pte_mats);
      auto pres_idxr = GenericIndexer(&pres_mat[0], pte_mats);
      PTESolverPT<decltype(eos_idxr), decltype(density_idxr), decltype(cache)> method(
              npte, eos_idxr, vfrac_tot, sie_tot, density_idxr, vfrac_idxr, sie_idxr,
              temp_idxr, pres_idxr, cache, solver_scratch.data());
      status = PTESolver(method);
      temperature = VolumeAverageMatArray(temp_idxr, vfrac_idxr, npte);
      pressure = VolumeAverageMatArray(pres_idxr, vfrac_idxr, npte);
    } else {
      // Single material
      temperature = eos_idxr[0].TemperatureFromDensityInternalEnergy(density_tot, sie_tot,
                                                                    cache[0]);
      // density-temperature lookup if preferred
      if (eos_idxr[0].PreferredInput() == thermalqs::density | thermalqs::temperature) {
        pressure = eos_idxr[0].PressureFromDensityTemperature(density_tot, temperature,
                                                                 cache[0]);
      } else {
        pressure = eos_idxr[0].PressureFromDensityInternalEnergy(density_tot, sie_tot,
                                                                 cache[0]);
      }
      density_mat[pte_mats[0]] = density_tot;
      sie_mat[pte_mats[0]] = sie_tot;
      status.converged = true;
      status.residual = 0.;
    }

    return status;
  }

  template<typename RealIndexer, typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_INLINE_FUNCTION
  SolverStatus GetStatesFromDensityTemperature(
      const Real density_tot, Real &sie_tot, Real &pressure, const Real temperature,
      RealIndexer &&density_mat, RealIndexer &&sie_mat, LambdaIndexer &&lambdas
  ) const {
    // Create temporary arrays
    std::array<Real, nmat_> mass_fracs{};
    std::array<Real, nmat_> vol_fracs{};
    std::array<size_t, nmat_> pte_mats{};
    std::array<Real, nmat_> pres_mat{};
    std::array<Real, nmat_> temp_mat{};

    // Extract the mass fractions from the lambdas
    PopulateMassFracArray(mass_fracs, lambdas); // Guarantees sum = 1

    // Compute volume fractions from mass fractions and densities. Mark
    // materials to participate in PTE
    size_t npte = 0;
    Real vfrac_tot = 0;
    sie_tot = 0;  // Initialize to some value for bad guesses or non-participating
                  // materials
    SetUpPTE(mass_fracs, vol_fracs, density_mat, pte_mats, npte, vfrac_tot, density_tot,
             sie_tot);

    // Initialize solver scratch memory (including cache)
    constexpr int pte_solver_scratch_size = PTESolverFixedTRequiredScratch(nmat_) + nmat_
        * MAX_NUM_LAMBDAS;
    std::array<Real, pte_solver_scratch_size> solver_scratch{};

    // Get cache from offsets into scratch
    const int neq = npte + 1;
    singularity::mix_impl::CacheAccessor cache(solver_scratch.data() +
                                               neq * (neq + 4) + 2 * npte);

    // Create indexers for the EOS
    const auto eos_idxr = GenericIndexer(CreateEOSArray(), pte_mats);


    // Initialize tempeature array to the fixed temperature. This likely isn't
    // necessary since the scalar temperature guess _should_ be used in favor
    // of any initial values from these arrays
    for (size_t i = 0; i < nmat_; i++) {
      temp_mat[i] = temperature;
    }

    // Solve for PTE state
    SolverStatus status;
    if (npte > 1) {
      // Create indexers for the various arrays. Note that we have to use
      // references since the PTE solver expects them all to be the same type
      auto density_idxr = GenericIndexer(&density_mat[0], pte_mats);
      auto vfrac_idxr = GenericIndexer(&vol_fracs[0], pte_mats);
      auto sie_idxr = GenericIndexer(&sie_mat[0], pte_mats);
      const auto temp_idxr = GenericIndexer(&temp_mat[0], pte_mats);
      auto pres_idxr = GenericIndexer(&pres_mat[0], pte_mats);
      PTESolverFixedT<decltype(eos_idxr), decltype(density_idxr), decltype(cache)>
      method(npte, eos_idxr, vfrac_tot, temperature, density_idxr, vfrac_idxr, sie_idxr,
             temp_idxr, pres_idxr, cache, solver_scratch.data());
      status = PTESolver(method);
      pressure = VolumeAverageMatArray(pres_idxr, vfrac_idxr, npte);
    } else {
      // Single material
      sie_tot = eos_idxr[0].InternalEnergyFromDensityTemperature(density_tot, temperature,
                                                                 cache[0]);
      // density-temperature lookup if preferred
      if (eos_idxr[0].PreferredInput() == thermalqs::density | thermalqs::temperature) {
        pressure = eos_idxr[0].PressureFromDensityTemperature(density_tot, temperature,
                                                                 cache[0]);
      } else {
        pressure = eos_idxr[0].PressureFromDensityInternalEnergy(density_tot, sie_tot,
                                                                 cache[0]);
      }
      density_mat[pte_mats[0]] = density_tot;
      sie_mat[pte_mats[0]] = sie_tot;
      status.converged = true;
      status.residual = 0.;
    }

    // Recalculate sie_tot as the sum of all material energies
    sie_tot = 0;
    for (size_t m = 0; m < nmat_; m++) {
      sie_tot += sie_mat[m] * mass_fracs[m];
    }

    return status;
  }

  template<typename RealIndexer, typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_INLINE_FUNCTION
  SolverStatus GetStatesFromDensityPressure(
      const Real density_tot, Real &sie_tot, const Real pressure, Real &temperature,
      RealIndexer &&density_mat, RealIndexer &&sie_mat, LambdaIndexer &&lambdas
  ) const {
    // Create temporary arrays
    std::array<Real, nmat_> mass_fracs{};
    std::array<Real, nmat_> vol_fracs{};
    std::array<size_t, nmat_> pte_mats{};
    std::array<Real, nmat_> pres_mat{};
    std::array<Real, nmat_> temp_mat{};

    // Extract the mass fractions from the lambdas
    PopulateMassFracArray(mass_fracs, lambdas); // Guarantees sum = 1

    // Compute volume fractions from mass fractions and densities. Mark
    // materials to participate in PTE
    size_t npte = 0;
    Real vfrac_tot = 0;
    sie_tot = 0;  // Initialize to some value for bad guesses or non-participating
                  // materials
    SetUpPTE(mass_fracs, vol_fracs, density_mat, pte_mats, npte, vfrac_tot, density_tot,
             sie_tot);

    // Initialize solver scratch memory (including cache)
    constexpr int pte_solver_scratch_size = PTESolverFixedPRequiredScratch(nmat_) + nmat_
        * MAX_NUM_LAMBDAS;
    std::array<Real, pte_solver_scratch_size> solver_scratch{};

    // Get cache from offsets into scratch
    const int neq = npte + 1;
    singularity::mix_impl::CacheAccessor cache(solver_scratch.data() +
                                               neq * (neq + 4) + 2 * npte);

    // Create indexers for the EOS
    const auto eos_idxr = GenericIndexer(CreateEOSArray(), pte_mats);


    // Initialize pressure array to the fixed pressure. This likely isn't
    // necessary since the pressure will be overwritten in the residual
    // calculation
    for (size_t i = 0; i < nmat_; i++) {
      pres_mat[i] = pressure;
    }

    // Solve for PTE state
    SolverStatus status;
    if (npte > 1) {
      // Create indexers for the various arrays. Note that we have to use
      // references since the PTE solver expects them all to be the same type.
      // Note that the pressure array can't be const since the PTE solver uses
      // it in the residual calculation
      auto density_idxr = GenericIndexer(&density_mat[0], pte_mats);
      auto vfrac_idxr = GenericIndexer(&vol_fracs[0], pte_mats);
      auto sie_idxr = GenericIndexer(&sie_mat[0], pte_mats);
      auto temp_idxr = GenericIndexer(&temp_mat[0], pte_mats);
      auto pres_idxr = GenericIndexer(&pres_mat[0], pte_mats);
      PTESolverFixedP<decltype(eos_idxr), decltype(density_idxr), decltype(cache)>
      method(npte, eos_idxr, vfrac_tot, pressure, density_idxr, vfrac_idxr, sie_idxr,
             temp_idxr, pres_idxr, cache, solver_scratch.data());
      status = PTESolver(method);
      temperature = VolumeAverageMatArray(temp_idxr, vfrac_idxr, npte);
    } else {
      // Single material, so use root find to find the temperature that achieves
      // the desired pressure at this density
      auto p_from_t = [&](const Real &t_i) {
        return eos_idxr[0].PressureFromDensityTemperature(density_tot, t_i, cache[0]);
      };

      // Use reference state to bound root
      Real r_rho{}, r_temp{}, r_sie{}, r_press{}, r_cv{}, r_bmod{}, r_dpde{}, r_dvdt{};
      eos_idxr[0].ValuesAtReferenceState(r_rho, r_temp, r_sie, r_press, r_cv, r_bmod,
                                         r_dpde, r_dvdt);
      Real temperature_root;
      RootFinding1D::Status root_status = RootFinding1D::findRoot(
          p_from_t, pressure, r_temp, 1.e-10 * r_temp, 1.e10 * r_temp, 1.e-12, 1.e-12,
          temperature_root);

      // Populate output variables
      temperature = temperature_root;
      sie_tot = eos_idxr[0].InternalEnergyFromDensityTemperature(density_tot, temperature,
                                                                 cache[0]);
      density_mat[pte_mats[0]] = density_tot;
      sie_mat[pte_mats[0]] = sie_tot;

      // Transfer root finder status to PTE solver status
      if (root_status == RootFinding1D::Status::SUCCESS) {
        status.converged = true;
        status.residual = 0.;
      } else {
        status.converged = false;
        status.residual = robust::ratio(p_from_t(temperature) - pressure, pressure);
      }
    }

    // Recalculate sie_tot as the sum of all material energies
    sie_tot = 0;
    for (size_t m = 0; m < nmat_; m++) {
      sie_tot += sie_mat[m] * mass_fracs[m];
    }

    return status;
  }

  template<typename RealIndexer, typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_INLINE_FUNCTION
  void GetStatesFromPressureTemperature(
      Real &density_tot, Real &sie_tot, const Real pressure, const Real temperature,
      RealIndexer &density_mat, RealIndexer &sie_mat, LambdaIndexer &lambdas
  ) const {
    // Extract the mass fractions from the lambdas
    std::array<Real, nmat_> mass_fracs{};
    PopulateMassFracArray(mass_fracs, lambdas);

    // Call DensityEnergyFromPressureTemperature() on all EOS
    callDensityEnergyFromPressureTemperatureAll(
        pressure, temperature, lambdas, density_mat, sie_mat,
        std::make_index_sequence<nmat_>{}
    );

    // Mass average the total density and internal energy
    sie_tot = 0;
    density_tot = 0;
    Real spvol_tot = 0;
    for (size_t m = 0; m < nmat_; m++) {
      spvol_tot += robust::ratio(mass_fracs, density_mat[m]);
      sie_tot += robust::ratio(mass_fracs, sie_mat[m]);
    }
    PORTABLE_REQUIRE(spvol_tot > 0., "Material volumes sum to zero");
    density_tot = robust::ratio(1.0, spvol_tot);
  }

  PORTABLE_FORCEINLINE_FUNCTION
  auto CreateEOSArray() {
    return std::apply(
      [this](auto&&... eos_models){
        return std::array<eosT_, nmat_>{ &eos_models...};
      },
      models_
    );
  }

  template<typename RealIndexer, typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_INLINE_FUNCTION
  void PopulateMassFracArray(RealIndexer &mass_fracs, LambdaIndexer &lambdas) const {

    // Wrap the implementation function so we don't have index sequences elsewhere
    assign_mass_fractions(mass_fracs, lambdas, std::make_index_sequence<nmat_>{});

    // Normalize the mass fractions. JHP: Is this necessary?
    Real mass_frac_sum = 0.;
    for (size_t m = 0; m < nmat_ ; m++) {
      mass_frac_sum += mass_fracs[m];
    }
    PORTABLE_REQUIRE(mass_frac_sum > 0., "Mass fractions sum to zero");
    for (size_t m = 0; m < nmat_; m++) {
      mass_fracs[m] /= mass_frac_sum;
    }
  }

  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, LambdaIndexer &&lambda) const {

    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, temperature;
    SolverStatus status = GetStatesFromDensityEnergy(rho, sie, pressure, temperature,
                                                     density_mat, sie_mat, lambda);

    return temperature;
  }
  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, temperature;
    SolverStatus status = GetStatesFromDensityEnergy(rho, sie, pressure, temperature,
                                                     density_mat, sie_mat, lambda);

    return pressure;
  }
  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_FUNCTION Real
  MinInternalEnergyFromDensity(const Real rho, LambdaIndexer &&lambda) const {
    // Technically this question is ill-defined since we need another thermodynamic
    // variable to constrain the volume fractions of the individual components
    // that combine to give the desired density. Instead, we return the mass-fraction
    // weighted average of the values returned by each EOS
    return massFracAverageQuantityAtOneState(
      [](auto const& eos, Real rho, Real dummy_value, LambdaIndexer lambda) {
          return eos.MinInternalEnergyFromDensity(rho, lambda);
      },
      rho,
      1.0, /* Dummy additional input */
      lambda
    );
  }
  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, temperature;
    SolverStatus status = GetStatesFromDensityEnergy(rho, sie, pressure, temperature,
                                                     density_mat, sie_mat, lambda);
    return massFracAverageQuantityAtManyStates(
      [](auto const& eos, Real const rho, Real const sie, Real const temperature,
         LambdaIndexer lambda) {
        if constexpr (eos.PreferredInput() == thermalqs::density | thermalqs::temperature) {
          return eos.EntropyFromDensityTemperature(rho, temperature, lambda);
        } else {
          return eos.EntropyFromDensityInternalEnergy(rho, sie, lambda);
        }
      },
      density_mat,
      sie_mat,
      temperature,
      lambda
    );
  }
  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, temperature;
    SolverStatus status = GetStatesFromDensityEnergy(rho, sie, pressure, temperature,
                                                     density_mat, sie_mat);

    return massFracAverageQuantityAtManyStates(
      [](auto const& eos, Real const rho, Real const sie, Real const temperature, LambdaIndexer &&lambda) {
        if constexpr (eos.PreferredInput() == thermalqs::density | thermalqs::temperature) {
          return eos.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
        } else {
          return eos.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
        }
      },
      density_mat,
      sie_mat,
      temperature,
      lambda
    );
  }
  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, temperature;
    SolverStatus status = GetStatesFromDensityEnergy(rho, sie, pressure, temperature,
                                                     density_mat, sie_mat, lambda);

    BulkModAvgT avg_funct{};
    return avg_funct(
      [](auto const& eos, Real const rho, Real const sie, Real const temperature,
         LambdaIndexer &&lambda) {
        if constexpr (eos.PreferredInput() == thermalqs::density | thermalqs::temperature) {
          return eos.BulkModulusFromDensityTemperature(rho, temperature, lambda);
        } else {
          return eos.BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
        }
      },
      density_mat,
      sie_mat,
      temperature,
      lambda,
      std::make_index_sequence<nmat_>{}
    );
  }
  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, temperature;
    SolverStatus status = GetStatesFromDensityEnergy(rho, sie, pressure, temperature,
                                                     density_mat, sie_mat, lambda);

    GrunesienAvgT avg_funct{};
    return avg_funct(
      [](auto const& eos, Real const rho, Real const sie, Real const temperature,
         LambdaIndexer &&lambda) {
        if constexpr (eos.PreferredInput() == thermalqs::density | thermalqs::temperature) {
          return eos.GruneisenParamFromDensityTemperature(rho, temperature, lambda);
        } else {
          return eos.GruneisenParamFromDensityInternalEnergy(rho, sie, lambda);
        }
      },
      density_mat,
      sie_mat,
      temperature,
      lambda,
      std::make_index_sequence<nmat_>{}
    );
  }
  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, sie;
    SolverStatus status = GetStatesFromDensityTemperature(rho, sie, pressure, temperature,
                                                          density_mat, sie_mat, lambda);

    return sie;
  }
  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, sie;
    SolverStatus status = GetStatesFromDensityTemperature(rho, sie, pressure, temperature,
                                                          density_mat, sie_mat, lambda);

    return pressure;
  }
  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_FUNCTION Real EntropyFromDensityTemperature(
      const Real rho, const Real temperature, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, sie;
    SolverStatus status = GetStatesFromDensityTemperature(rho, sie, pressure, temperature,
                                                          density_mat, sie_mat, lambda);
    return massFracAverageQuantityAtManyStates(
      [](auto const& eos, Real const rho, Real const sie, Real const temperature, LambdaIndexer lambda) {
        if constexpr (eos.PreferredInput() == thermalqs::density | thermalqs::temperature) {
          return eos.EntropyFromDensityTemperature(rho, temperature, lambda);
        } else {
          return eos.EntropyFromDensityInternalEnergy(rho, sie, lambda);
        }
      },
      density_mat,
      sie_mat,
      temperature,
      lambda
    );
  }
  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, sie;
    SolverStatus status = GetStatesFromDensityTemperature(rho, sie, pressure, temperature,
                                                          density_mat, sie_mat, lambda);
    return massFracAverageQuantityAtManyStates(
      [](auto const& eos, Real const rho, Real const sie, Real const temperature, LambdaIndexer lambda) {
        if constexpr (eos.PreferredInput() == thermalqs::density | thermalqs::temperature) {
          return eos.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
        } else {
          return eos.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
        }
      },
      density_mat,
      sie_mat,
      temperature,
      lambda
    );
  }
  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, sie;
    SolverStatus status = GetStatesFromDensityTemperature(rho, sie, pressure, temperature,
                                                          density_mat, sie_mat, lambda);
    BulkModAvgT avg_funct{};
    return avg_funct(
      [](auto const& eos, Real const rho, Real const sie, Real const temperature, LambdaIndexer lambda) {
        if constexpr (eos.PreferredInput() == thermalqs::density | thermalqs::temperature) {
          return eos.BulkModulusFromDensityTemperature(rho, temperature, lambda);
        } else {
          return eos.BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
        }
      },
      density_mat,
      sie_mat,
      temperature,
      lambda,
      std::make_index_sequence<nmat_>{}
    );
  }
  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, sie;
    SolverStatus status = GetStatesFromDensityTemperature(rho, sie, pressure, temperature,
                                                          density_mat, sie_mat, lambda);
    GrunesienAvgT avg_funct{};
    return avg_funct(
      [](auto const& eos, Real const rho, Real const sie, Real const temperature, LambdaIndexer lambda) {
        if constexpr (eos.PreferredInput() == thermalqs::density | thermalqs::temperature) {
          return eos.GruneisenParamFromDensityTemperature(rho, temperature, lambda);
        } else {
          return eos.GruneisenParamFromDensityInternalEnergy(rho, sie, lambda);
        }
      },
      density_mat,
      sie_mat,
      temperature,
      lambda,
      std::make_index_sequence<nmat_>{}
    );
  }
  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_FUNCTION void FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                 Real &cv, Real &bmod, const unsigned long output,
                                 LambdaIndexer &&lambda) const {
    const unsigned long input = ~output;
    PORTABLE_ALWAYS_REQUIRE(input != thermalqs::none, "No input provided");

    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};

    if (input & thermalqs::pressure && input & thermalqs::temperature) {
      PORTABLE_REQUIRE(
        output & thermalqs::density && output & thermalqs::specific_internal_energy,
        "Internal energy and density must be defined as outputs when pressure and "
        "temperature are inputs");
      GetStatesFromPressureTemperature(energy, rho, press, temp, density_mat,
                                       sie_mat, lambda);
    } else if (input & thermalqs::density && input & thermalqs::temperature) {
      PORTABLE_REQUIRE(
        output & thermalqs::pressure && output & thermalqs::specific_internal_energy,
        "Internal energy and pressure must be defined as outputs when density and "
        "temperature are inputs");
      GetStatesFromDensityTemperature(energy, rho, press, temp, density_mat,
                                      sie_mat, lambda);
    } else if (input & thermalqs::density && input & thermalqs::specific_internal_energy) {
      PORTABLE_REQUIRE(
        output & thermalqs::pressure && output & thermalqs::specific_internal_energy,
        "Temperature and pressure must be defined as outputs when density and "
        "internal energy are inputs");
      GetStatesFromDensityEnergy(energy, rho, press, temp, density_mat,
                                 sie_mat, lambda);
    } else if (input & thermalqs::density ** input * thermalqs::pressure) {
      PORTABLE_REQUIRE(
        output & thermalqs::pressure && output & thermalqs::specific_internal_energy,
        "Temperature and pressure must be defined as outputs when density and "
        "internal energy are inputs");
      GetStatesFromDensityPressure(energy, rho, press, temp, density_mat,
                                   sie_mat, lambda);
    }

    if (output & thermalqs::specific_heat) {
      cv = massFracAverageQuantityAtManyStates(
        [](auto const& eos, Real const rho, Real const sie, Real const temperature,
           LambdaIndexer lambda) {
          if constexpr (eos.PreferredInput() == thermalqs::density | thermalqs::temperature) {
            return eos.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
          } else {
            return eos.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
          }
        },
        density_mat,
        sie_mat,
        temperature,
        lambda,
      )
    }

    if (output & thermalqs::bulk_modulus) {
      BulkModAvgT avg_funct{};
      bmod = avg_funct(
        [](auto const& eos, Real const rho, Real const sie, Real const temperature,
           LambdaIndexer lambda) {
          if constexpr (eos.PreferredInput() == thermalqs::density | thermalqs::temperature) {
            return eos.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
          } else {
            return eos.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
          }
        },
        density_mat,
        sie_mat,
        temperature,
        lambda,
        std::make_index_sequence<nmat_>{}
      )
    }
  }
  constexpr static inline int nlambda() noexcept {
    // Sum the lambda requirements from each EOS with the number of mass
    // fractions required.
    // TODO: Is this really what we want to do? Really some EOS will share
    // lambda variables and we probably want the superset of all possible lambda
    // variables
    return (EOSModelsT::nlambda() + ...) + nmat_;
  }

  static constexpr unsigned long PreferredInput() {
    // Store the "vote" for each EOS on what the preferred input should be
    std::array<unsigned long, nmat_> desired_inputs{{ EOSModelsT::PreferredInput()... }};

    // Assume the best input is the mode of all the EOS preferred inputs
    unsigned long chosen_inputs = thermalqs::none;
    size_t chosen_count = 0;
    for (size_t i = 0; i < nmat_; i++) {
      size_t this_choice_count = 0;
      // Loop over other elements to count how many are equal to this one's
      // desired inputs
      for (size_t j = 0; j < nmat_; ++j) {
        if (desired_inputs[i] == desired_inputs[j]) {
          ++this_choice_count;
        }
      }
      if (this_choice_count > chosen_count) {
        // Replace the current leader if needed
        chosen_count = this_choice_count;
        chosen_inputs = desired_inputs[i];
      }
    }

    // Assertion: there should always be an appropriate choice of inputs
    PORTABLE_REQUIRE(chosen_inputs != thermalqs::none,
                     "PreferredInput: No valid preferred input found");

    return chosen_inputs;
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
      scratch_size += PTESolverFixedTRequiredScratchInBytes(nmat_);
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
      scratch_size += PTESolverPTRequiredScratchInBytes(nmat_);
      break;
    case "MinInternalEnergyFromDensity":
      // No scratch needed since not a PTE solve
      break;
    case "FillEos":
      // Without knowing the exact inputs, we have to assume a worst-case
      // scenario
      scratch_size += std::max(PTESolverFixedTRequiredScratchInBytes(nmat_),
                               std::max(PTESolverPTRequiredScratchInBytes(nmat_),
                                        PTESolverFixedPRequiredScratchInBytes(nmat_))
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
    (scratch_size + ... + EOSModelsT::scratch_size(method, nelements));

    // Contribution from lambdas
    scratch_size += nmat_ * nlambda() * nelements;

    return scratch_size;
  }

  static inline unsigned long max_scratch_size(unsigned int nelements) {
    // Worst case scenario for PTE solve
    unsigned long scratch_size =
        std::max(PTESolverFixedTRequiredScratchInBytes(nmat_),
                 std::max(PTESolverPTRequiredScratchInBytes(nmat_),
                          PTESolverFixedPRequiredScratchInBytes(nmat_))
      );

    // Since execution could be asynchronous, require that all scratch for PTE
    // be available for all states simultaneously. This mirrors the EOS scratch
    // size behavior for e.g. EOSPAC
    scratch_size *= nelements;

    // Add in the maximum requirements of the individual EOS. Note the
    // `nelements` argument
    (scratch_size + ... + EOSModelsT::max_scratch_size(nelements));

    // Contribution from lambdas... does this change the meaning of
    // MAX_NUM_LAMBDAS?
    scratch_size += nmat_ * nlambda() * nelements;

    return scratch_size;

  }

  PORTABLE_FUNCTION void PrintParams() const {
    // print the parameters for each EOS via fold expression
    std::apply([](const auto&... eos_models) {
      (eos_models.PrintParams(), ...);
    }, models_);
  }

  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       LambdaIndexer &&lambda, Real &rho, Real &sie) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    GetStatesFromPressureTemperature(sie, rho, press, temp, density_mat,
                                     sie_mat, lambda);

  }

  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_FUNCTION void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie,
                                                Real &press, Real &cv, Real &bmod,
                                                Real &dpde, Real &dvdt,
                                                LambdaIndexer &&lambda) const {
    // Assume a reference state of 298K and 1 bar
    temp = 298;
    press = 1e5;

    // Get the mixture properties
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    GetStatesFromPressureTemperature(sie, rho, press, temp, density_mat,
                                     sie_mat, lambda);

    // Mass-fraction average the gruneisen and heat capacity
    cv = massFracAverageQuantityAtManyStates(
      [](auto const& eos, Real const rho, Real const sie, Real const temperature,
         LambdaIndexer lambda) {
        if constexpr (eos.PreferredInput() == thermalqs::density | thermalqs::temperature) {
          return eos.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
        } else {
          return eos.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
        }
      },
      density_mat,
      sie_mat,
      temp,
      lambda,
    )

    dpde = rho * massFracAverageQuantityAtManyStates(
      [](auto const& eos, Real const rho, Real const sie, Real const temperature,
         LambdaIndexer lambda) {
        if constexpr (eos.PreferredInput() == thermalqs::density | thermalqs::temperature) {
          return eos.GruneisenParamFromDensityTemperature(rho, temperature, lambda);
        } else {
          return eos.GruneisenParamFromDensityInternalEnergy(rho, sie, lambda);
        }
      },
      density_mat,
      sie_mat,
      temp,
      lambda,
    )

    // Extract the mass fractions from the lambdas
    std::array<Real, nmat_> mass_fracs{};
    PopulateMassFracArray(mass_fracs, lambda);

    // Populate material derivative information
    std::array<Real, nmat_> bmod_mat{};
    std::array<Real, nmat_> gruneisen_mat{};
    std::array<Real, nmat_> cv_mat{};
    // TODO: Use fold expressions instead

    // Mass fraction average dpde, dvdt, and cv. Volume fraction average the
    // bulk modulus
    cv = 0.;
    dpde = 0.;
    dvdt = 0.;
    for (size_t m = 0; m < nmat_; m++) {
      Real const this_dpde = density_mat[m] * gruneisen_mat[m];
      Real const this_dvdt = gruneisen_mat[m] * cv_mat[m] / bmod_mat[m];
      cv += mass_fracs[m] * cv_mat[m];
      dpde += mass_fracs[m] * this_dpde;
      dvdt += mass_fracs[m] * this_dvdt;
    }

    // Bulk modulus averaging is injected because there are multiple ways to do
    // it. Harmonic average is appropriate, but may not always be appropriate.
    constexpr auto bmod_avg_func = BulkModAvgT();
    bmod = bmod_avg_func(mass_fracs, density_mat, bmod_mat, rho, nmat_);
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumDensity() const {
    // TODO
    // Max over the contributions from each EOS
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumTemperature() const {
    // TODO
    // Max over the contributions from each EOS
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MaximumDensity() const {
    // TODO
    // Min over the contributions from each EOS
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumPressure() const {
    // TODO
    // Max over the contributions from each EOS
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MaximumPressureAtTemperature(const Real temp) const {
    // TODO
    // Min over the contributions from each EOS
  }

  // 3T member functions
  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicMassFromDensityTemperature(
      const Real rho, const Real temperature,
      LambdaIndexer &&lambda = static_cast<Real *>(nullptr)) const {
    // TODO
    // Average over the contributions from each EOS
  }
  template <typename LambdaIndexer, LAMBDA_HAS_MASS_FRACTION_INDICES(LambdaIndexer)>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicNumberFromDensityTemperature(
      const Real rho, const Real temperature,
      LambdaIndexer &&lambda = static_cast<Real *>(nullptr)) const {
    // TODO
    // Average over the contributions from each EOS
  }
  PORTABLE_INLINE_FUNCTION
  Real MeanAtomicMass() const {
    // TODO
    // Average over the contributions from each EOS
  }
  PORTABLE_INLINE_FUNCTION
  Real MeanAtomicNumber() const {
    // TODO
    // Average over the contributions from each EOS
  }

  // Modifier member functions
  static inline constexpr bool IsModified() {
    return true;
  }
  inline constexpr const std::tuple<EOSModelsT...> UnmodifyOnce() {
    return models_;
  }

  // Dynamic and shared memory member functions (normally included via macro)
  std::size_t DynamicMemorySizeInBytes() const {
    // TODO
    // Add the contributions from each EOS
  }
  std::size_t SharedMemorySizeInBytes() const {
    // TODO
    // Add the contributions from each EOS
  }
  std::size_t DumpDynamicMemory(char *dst) {
    // TODO
    // It seems like we'll need to dump from each EOS and then advance the
    // pointer to the next stretch of available memory
  }
  std::size_t SetDynamicMemory(char *src,
                               const SharedMemSettings &stngs = DEFAULT_SHMEM_STNGS) {
    // TODO
    // I think I should loop through the source memory and provide each EOS with
    // a pointer to where its memory begins
  }
  constexpr bool AllDynamicMemoryIsShareable() const {
    return std::apply(
    [](auto const&... eos){
      // AND operation across all individual models. Returns true IFF all models
      // have sharable dynamic memory
      return (eos.AllDynamicMemoryIsShareable() && ...);
    },
    models_
  );
  }
};

// Factory function helper to group EOS together and provide default averaging
// behavior
template<
      typename BulkModAvgT = VolumeFracHarmonicAverageFunctor
    , typename GrunesienAvgT = VolumeFracHarmonicAverageFunctor
    , typename... EOSModelsT
>
inline
auto make_MultiEOS(EOSModelsT&&... eos_models) {
  return MultiEOS<
    std::remove_cv_t<std::remove_reference_t<BulkModAvgT>>,
    std::remove_cv_T<std::remove_reference_t<GrunesienAvgT>>,
    std::remove_cv_t<std::remove_reference_t<EOSModelsT>>...
  >(
    std::forward<EOSModelsT>(eos_models)...
  );
}

// Factory function helper to group EOS together and provide default averaging
// behavior
template<
      typename BulkModAvgT = VolumeFracHarmonicAverageFunctor
    , typename GrunesienAvgT = VolumeFracHarmonicAverageFunctor
    , typename... EOSModelsT
>
inline
auto make_MultiEOS(Real mass_frac_cutoff_in, EOSModelsT&&... eos_models) {
  return MultiEOS<
    std::remove_cv_t<std::remove_reference_t<BulkModAvgT>>,
    std::remove_cv_T<std::remove_reference_t<GrunesienAvgT>>,
    std::remove_cv_t<std::remove_reference_t<EOSModelsT>>...
  >(
    mass_frac_cutoff_in,
    std::forward<EOSModelsT>(eos_models)...
  );
}

} // namespace singularity

#endif // #ifndef _SINGULARITY_EOS_EOS_MULTI_EOS_
