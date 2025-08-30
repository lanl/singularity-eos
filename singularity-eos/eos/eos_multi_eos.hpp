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
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <tuple>
#include <utility>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/generic_indexer.hpp>
#include <singularity-eos/base/indexable_types.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>
#include <singularity-eos/base/tuple_utils.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>
#include <singularity-eos/eos/eos_base.hpp>
#include <singularity-eos/eos/eos_variant.hpp>

namespace singularity {

// This functor essentially wraps a fold expression so that each term in the
// fold can be injected via the functor f. The functor should take the EOS
// followed by the density, internal energy, and temperature so that it can
// perform either density-energy or density-temperature lookups with the EOS.
// The functor should also return the _inverse_ term so that a harmonic average
// is properly formed. Note that similar averaging functors would need to obey
// this basic API.
struct VolumeFracHarmonicAverageFunctor {
  template <typename Func, typename EOSmodelsT, typename RealIndexer,
            typename LambdaIndexer, std::size_t... Is>
  constexpr Real operator()(Func &&f, EOSmodelsT models, RealIndexer density_arr,
                            RealIndexer sie_arr, Real density, Real temperature,
                            LambdaIndexer const &lambda,
                            std::index_sequence<Is...>) const {
    // Note that the function should appropriately invert any quantity that
    return robust::ratio(1.0,
                         ((/* volume fraction */
                           lambda[IndexableTypes::MassFraction<Is>{}] * density /
                           density_arr[Is]
                           /* Inverse value to be volume-averaged */
                           * robust::ratio(1.0, f(std::get<Is>(models), density_arr[Is],
                                                  sie_arr[Is], temperature, lambda))) +
                          ...));
  }
};

struct MassFracAverageFunctor {
  template <typename Func, typename EOSmodelsT, typename RealIndexer,
            typename LambdaIndexer, std::size_t... Is>
  constexpr Real operator()(Func &&f, EOSmodelsT models, RealIndexer density_arr,
                            RealIndexer sie_arr, [[maybe_unused]] Real density,
                            Real temperature, LambdaIndexer const &lambda,
                            std::index_sequence<Is...>) const {
    return ((lambda[IndexableTypes::MassFraction<Is>{}]
             /* Value to be volume-avreaged */
             *
             f(std::get<Is>(models), density_arr[Is], sie_arr[Is], temperature, lambda)) +
            ...);
  }
};

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

template <typename BulkModAvgT = VolumeFracHarmonicAverageFunctor,
          typename GruneisenAvgT = VolumeFracHarmonicAverageFunctor,
          typename... EOSModelsT>
class MultiEOS : public EosBase<MultiEOS<BulkModAvgT, GruneisenAvgT, EOSModelsT...>> {
 private:
  Real mass_frac_cutoff_;
  std::tuple<EOSModelsT...> models_;
  using eosT_ = Variant<EOSModelsT *...>;

  static constexpr std::size_t nmat_ = sizeof...(EOSModelsT); // for convenience

  // NOTE: These functions are private since they're details of the
  // implementation

  // Implementation function to unpack the lambda mass fraction indices into an
  // array. I believe this has to be separate from the public implementation to
  // get the 'Is' template parameter to specify the index sequence
  template <typename RealIndexer, typename LambdaIndexer, size_t... Is,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  constexpr void assign_mass_fractions(RealIndexer mass_fracs, const LambdaIndexer lambda,
                                       std::index_sequence<Is...>) {
    ((mass_fracs[Is] = lambda[IndexableTypes::MassFraction<Is>{}]), ...);
  }

  // Implementation function to call DensityEnergyFromPressureTemperature() on
  // each EOS for a common pressure and temperature and populate material arrays
  template <typename RealIndexer, typename LambdaIndexer, size_t... Is,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  constexpr void callDensityEnergyFromPressureTemperatureAll(
      const Real pressure, const Real temperature, LambdaIndexer &&lambdas,
      RealIndexer &density_mat, RealIndexer &energy_mat, std::index_sequence<Is...>) {
    ((std::get<Is>(models_).DensityEnergyFromPressureTemperature(
         pressure, temperature, std::forward<LambdaIndexer>(lambdas), density_mat[Is],
         energy_mat[Is])),
     ...);
  }

  // Implementation of mass fraction average over all EOS but at a single state
  // given by input1 and input2
  template <typename Func, typename LambdaIndexer, std::size_t... Is>
  Real massFracAverageQuantityAtOneState(Func &&f, Real input1, Real input2,
                                         LambdaIndexer const &lambda,
                                         std::index_sequence<Is...>) const {
    // Note that the function signature should take the EOS as its first argument
    // and then the two independent variables and the lambda
    return ((lambda[IndexableTypes::MassFraction<Is>{}] *
             f(std::get<Is>(models_), input1, input2, lambda)) +
            ...);
  }

  // TODO: Is this wrapper necessary? Just use bare functor...
  // Implementation of mass fraction average over all EOS, each at its own state
  // Both internal energy and temperature are provided along with density so that
  // the function wrapper can include logic to decide which lookup to use
  template <typename Func, typename RealIndexer, typename LambdaIndexer,
            std::size_t... Is>
  Real massFracAverageQuantityAtManyStates(Func &&f, RealIndexer density_arr,
                                           RealIndexer sie_arr, Real temperature,
                                           LambdaIndexer const &lambda) const {
    MassFracAverageFunctor avg_funct{};
    // Note that the f function signature should take the EOS as its first argument
    // followed by the lookup state and lambda
    return avg_funct(f, models_, density_arr, sie_arr, temperature, /* density */ 1.0,
                     density_arr, sie_arr, lambda, std::make_index_sequence<nmat_>{});
  }

  // Implementation function to call BulkModulusFromDensityTemperature()
  // on each EOS with a common temperature and array of densities
  template <typename RealIndexer, typename LambdaIndexer, size_t... Is,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  constexpr void getMaterialBulkModulusFromDensityTemperature(
      const Real temperature, LambdaIndexer &&lambdas, RealIndexer &density_mat,
      RealIndexer &bmod_mat, std::index_sequence<Is...>) {
    ((bmod_mat[Is] = std::get<Is>(models_).BulkModulusFromDensityTemperature(
          density_mat[Is], temperature, std::forward<LambdaIndexer>(lambdas))),
     ...);
  }

  // Implementation function to call BulkModulusFromDensityInternalEnergy()
  // on each EOS with arrays of densities and internal energies
  template <typename RealIndexer, typename LambdaIndexer, size_t... Is,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  constexpr void getMaterialBulkModulusFromDensityInternalEnergy(
      LambdaIndexer &&lambdas, RealIndexer &density_mat, RealIndexer &sie_mat,
      RealIndexer &bmod_mat, std::index_sequence<Is...>) {
    ((bmod_mat[Is] = std::get<Is>(models_).BulkModulusFromDensityInternalEnergy(
          density_mat[Is], sie_mat[Is], std::forward<LambdaIndexer>(lambdas))),
     ...);
  }

  // Generic "max" across models for arbitrary callable f(model, args...)
  template <typename Func, typename... Args>
  Real max_over_models(Func &&f, Args &&...args) const {
    static_assert((std::is_invocable_v<Func &, const EOSModelsT &, Args...> && ...),
                  "Callable must be invocable with (const EOS&, Args...) for every EOS");

    Real max_val = std::numeric_limits<Real>::lowest();

    // Fold expression over all EOS models to return maximum
    std::apply(
        [&](auto const &...eos) {
          ((max_val = (std::max)(max_val,
                                 static_cast<Real>(f(eos, std::forward<Args>(args)...)))),
           ...);
        },
        models_);
    return max_val;
  }

  // Generic "min" across models for arbitrary callable f(model, args...)
  template <typename Func, typename... Args>
  Real min_over_models(Func &&f, Args &&...args) const {
    static_assert((std::is_invocable_v<Func &, const EOSModelsT &, Args...> && ...),
                  "Callable must be invocable with (const EOS&, Args...) for every EOS");

    Real min_val = std::numeric_limits<Real>::max();

    // Fold expression over all EOS models to return maximum
    std::apply(
        [&](auto const &...eos) {
          ((min_val = (std::min)(min_val,
                                 static_cast<Real>(f(eos, std::forward<Args>(args)...)))),
           ...);
        },
        models_);
    return min_val;
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
    std::string all_models = (std::string() + ... + EOSModelsT::EosPyType());
    return "MultiEOS" + all_models;
  }

  MultiEOS() = default;

  template <typename... EOSModelsT_>
  MultiEOS(Real mass_frac_cutoff_in, EOSModelsT_ &&...eos_models)
      : mass_frac_cutoff_{std::forward<Real>(mass_frac_cutoff_in)},
        models_(std::forward<EOSModelsT_>(eos_models)...) {
    CheckParams();
  }

  template <typename... EOSModelsT_>
  MultiEOS(EOSModelsT_ &&...eos_models)
      : MultiEOS(
            /*mass_frac_cutoff=*/1e-8, std::forward<EOSModelsT_>(eos_models)...) {}

  // Overload that makes copying the models easier (e.g. GetOnDevice()). Make
  // explicit to avoid accidentally storing a tuple of a tuple
  explicit MultiEOS(Real mass_frac_cutoff_in, std::tuple<EOSModelsT...> model_tuple)
      : models_(std::move(model_tuple)), mass_frac_cutoff_{mass_frac_cutoff_in} {
    CheckParams();
  }

  PORTABLE_INLINE_FUNCTION void CheckParams() const {
    // Again we're using a fold expression now with std::apply to call
    // CheckParams() on each model (C++17)
    PORTABLE_ALWAYS_REQUIRE(mass_frac_cutoff_ > 0.,
                            "Mass fracton cutoff must be non-negative");
    PORTABLE_ALWAYS_REQUIRE(mass_frac_cutoff_ < 1., "Mass fractons must be less than 1");
    std::apply([](const auto &...eos_models) { (eos_models.CheckParams(), ...); },
               models_);
  }

  auto GetOnDevice() {
    // Because GetOnDevice() returns a copy of the EOS, we can't really just
    // use std::apply, and instead we need to also populate a new copy of the
    // model tuple
    auto models_on_device =
        tuple_transform(models_, [](const auto &eos) { return eos.GetOnDevice(); });
    return MultiEOS<BulkModAvgT, EOSModelsT...>(mass_frac_cutoff_, models_on_device);
  }

  inline void Finalize() {
    // All the fold expressions!
    std::apply([](const auto &...eos_models) { (eos_models.Finalize(), ...); }, models_);
  }

  // Compute volume fractions and possibily material densities from mass
  // fractions and total density
  template <typename MassFracT, typename VolFracT, typename DensityT, typename SieT,
            typename PTEMatT, typename EosArrT>
  PORTABLE_FORCEINLINE_FUNCTION void
  SetUpPTE(MassFracT mass_fracs, VolFracT vol_fracs, DensityT density_mat, SieT sie_mat,
           PTEMatT pte_mats, size_t &npte, Real &vfrac_tot, const Real density_tot,
           const Real sie_tot, EosArrT eos_arr) const {
    npte = 0;
    for (size_t m = 0; m < nmat_; m++) {

      // Check whether a material is present or not
      if (mass_fracs[m] <= mass_frac_cutoff_) {
        // Make sure sane values are returned for material density/energy even
        // if material isn't participating. These can be overwritten with more
        // accurate values alter, but importantly we ensure the volume fraction
        // is zero for this material
        vol_fracs[m] = 0.;
        density_mat[m] = density_tot;
        sie_mat[m] = sie_tot;
        continue;
      }

      // Material will participate in PTE
      pte_mats[npte] = m;
      npte += 1;

      // Calculate volume fraction from input density to decide whether it was
      // sane or not. The volume fraction lower bound is given by the maximum
      // density from the EOS
      const auto vol_frac_test =
          mass_fracs[m] * robust::ratio(density_tot, density_mat[m]);
      const auto vol_frac_cutoff =
          mass_frac_cutoff_ * robust::ratio(density_tot, eos_arr[m].MaximumDensity());
      if (!std::isnormal(density_mat[m]) || vol_frac_test > 1.0 ||
          vol_frac_test < vol_frac_cutoff) {
        // Bad density guess... fall back to average value
        density_mat[m] = density_tot;
        sie_mat[m] = sie_tot;
        vol_fracs[m] = mass_fracs[m]; // Consistent with density_mat[m] = density_tot
      } else {
        vol_fracs[m] = vol_frac_test;
      }
    }

    vfrac_tot = math_utils::sum_neumaier(vol_fracs, nmat_);

    // JHP: I think we should always test this unless it is proven to be a performance
    // bottleneck
    PORTABLE_ALWAYS_REQUIRE(npte > 0, "No material mass fraction is greater than cutoff");
  }

  // Helper function to volume-average an array
  template <typename arrT, typename vfracT>
  PORTABLE_FORCEINLINE_FUNCTION Real VolumeAverageMatArray(arrT arr, vfracT vfracs,
                                                           int num = nmat_) {
    Real avg = 0;
    for (size_t i = 0; i < num; i++) {
      avg += vfracs[i] * arr[i];
    }
    return avg;
  }

  // Wrappers for the various types of independent variables for the mixture
  template <typename RealIndexer, typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_INLINE_FUNCTION SolverStatus GetStatesFromDensityEnergy(
      const Real density_tot, const Real sie_tot, Real &pressure, Real &temperature,
      RealIndexer &&density_mat, RealIndexer &&sie_mat, LambdaIndexer &&lambdas,
      bool small_mass_mat_consistency = false) const {
    // Create temporary arrays
    std::array<Real, nmat_> mass_fracs{};
    std::array<Real, nmat_> vol_fracs{};
    std::array<size_t, nmat_> pte_mats{};
    std::array<Real, nmat_> pres_mat{};
    std::array<Real, nmat_> temp_mat{};

    // Extract the mass fractions from the lambdas
    PopulateMassFracArray(mass_fracs, lambdas); // Guarantees sum = 1

    // Create array of pointers to EOS models
    const auto eos_arr = CreateEOSArray();

    // Compute volume fractions from mass fractions and densities. Mark
    // materials to participate in PTE
    size_t npte = 0;
    Real vfrac_tot = 0;
    SetUpPTE(mass_fracs, vol_fracs, density_mat, sie_mat, pte_mats, npte, vfrac_tot,
             density_tot, sie_tot, eos_arr);

    // Initialize solver scratch memory (including cache)
    constexpr int pte_solver_scratch_size =
        PTESolverPTRequiredScratch(nmat_) + nmat_ * MAX_NUM_LAMBDAS;
    std::array<Real, pte_solver_scratch_size> solver_scratch{};

    // Get cache from offsets into scratch
    const int neq = npte + 1;
    singularity::mix_impl::CacheAccessor cache(solver_scratch.data() + neq * (neq + 4) +
                                               2 * npte);

    // Create indexers for the EOS
    const auto eos_idxr = GenericIndexer(eos_arr, pte_mats);

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
      if (eos_idxr[0].PreferredInput() == (thermalqs::density | thermalqs::temperature)) {
        pressure = eos_idxr[0].PressureFromDensityTemperature(density_tot, temperature,
                                                              cache[0]);
      } else {
        pressure =
            eos_idxr[0].PressureFromDensityInternalEnergy(density_tot, sie_tot, cache[0]);
      }
      density_mat[pte_mats[0]] = density_tot;
      sie_mat[pte_mats[0]] = sie_tot;
      status.converged = true;
      status.residual = 0.;
    }

    // Since this is going to be rarely used and could be unneccesarily
    // expensive, it's wrapped in an input flag.
    // Loop through non-participating materials and ensure their energy and
    // density are consistent with P-T state. We also skip materials for which
    // the PTE pressure is below the minimum pressure (i.e. tension)
    if (small_mass_mat_consistency and npte != nmat_) {
      for (size_t m = 0; m < nmat_; m++) {
        // Zero volume fraction materials didn't participate in PTE
        if (vol_fracs[m] > 0 || eos_arr[m].MinimumPressure() > pressure) {
          continue;
        }
        eos_arr[m].DensityEnergyFromPressureTemperature(
            pressure, temperature, std::forward<LambdaIndexer>(lambdas), density_mat[m],
            sie_mat[m]);
      }
    }

    return status;
  }

  template <typename RealIndexer, typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_INLINE_FUNCTION SolverStatus GetStatesFromDensityTemperature(
      const Real density_tot, Real &sie_tot, Real &pressure, const Real temperature,
      RealIndexer &&density_mat, RealIndexer &&sie_mat, LambdaIndexer &&lambdas) const {
    // Create temporary arrays
    std::array<Real, nmat_> mass_fracs{};
    std::array<Real, nmat_> vol_fracs{};
    std::array<size_t, nmat_> pte_mats{};
    std::array<Real, nmat_> pres_mat{};
    std::array<Real, nmat_> temp_mat{};

    // Extract the mass fractions from the lambdas
    PopulateMassFracArray(mass_fracs, lambdas); // Guarantees sum = 1

    // Create array of pointers to EOS models
    const auto eos_arr = CreateEOSArray();

    // Compute volume fractions from mass fractions and densities. Mark
    // materials to participate in PTE
    size_t npte = 0;
    Real vfrac_tot = 0;
    sie_tot = 0; // Initialize to some value for bad guesses or non-participating
                 // materials
    SetUpPTE(mass_fracs, vol_fracs, density_mat, pte_mats, npte, vfrac_tot, density_tot,
             sie_tot, eos_arr);

    // Initialize solver scratch memory (including cache)
    constexpr int pte_solver_scratch_size =
        PTESolverFixedTRequiredScratch(nmat_) + nmat_ * MAX_NUM_LAMBDAS;
    std::array<Real, pte_solver_scratch_size> solver_scratch{};

    // Get cache from offsets into scratch
    const int neq = npte + 1;
    singularity::mix_impl::CacheAccessor cache(solver_scratch.data() + neq * (neq + 4) +
                                               2 * npte);

    // Create indexers for the EOS
    const auto eos_idxr = GenericIndexer(eos_arr, pte_mats);

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
      PTESolverFixedT<decltype(eos_idxr), decltype(density_idxr), decltype(cache)> method(
          npte, eos_idxr, vfrac_tot, temperature, density_idxr, vfrac_idxr, sie_idxr,
          temp_idxr, pres_idxr, cache, solver_scratch.data());
      status = PTESolver(method);
      pressure = VolumeAverageMatArray(pres_idxr, vfrac_idxr, npte);
    } else {
      // Single material
      sie_tot = eos_idxr[0].InternalEnergyFromDensityTemperature(density_tot, temperature,
                                                                 cache[0]);
      // density-temperature lookup if preferred
      if (eos_idxr[0].PreferredInput() == (thermalqs::density | thermalqs::temperature)) {
        pressure = eos_idxr[0].PressureFromDensityTemperature(density_tot, temperature,
                                                              cache[0]);
      } else {
        pressure =
            eos_idxr[0].PressureFromDensityInternalEnergy(density_tot, sie_tot, cache[0]);
      }
      density_mat[pte_mats[0]] = density_tot;
      sie_mat[pte_mats[0]] = sie_tot;
      status.converged = true;
      status.residual = 0.;
    }

    // Loop through non-participating materials and ensure their energy and
    // density are consistent with P-T state. Otherwise, the material energy for
    // small-mass-fraction materials would be zero and this would change the
    // energy sum below. If a material's minimum pressure is above the PTE
    // pressure, then we keep the zero energy and cell density
    if (npte != nmat_) {
      for (size_t m = 0; m < nmat_; m++) {
        // Zero volume fraction materials didn't participate in PTE
        if (vol_fracs[m] > 0 || eos_arr[m].MinimumPressure() > pressure) {
          continue;
        }
        eos_arr[m].DensityEnergyFromPressureTemperature(
            pressure, temperature, std::forward<LambdaIndexer>(lambdas), density_mat[m],
            sie_mat[m]);
      }
    }

    // Recalculate sie_tot as the sum of all material energies
    sie_tot = 0;
    for (size_t m = 0; m < nmat_; m++) {
      sie_tot += sie_mat[m] * mass_fracs[m];
    }

    return status;
  }

  template <typename RealIndexer, typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_INLINE_FUNCTION SolverStatus GetStatesFromDensityPressure(
      const Real density_tot, Real &sie_tot, const Real pressure, Real &temperature,
      RealIndexer &&density_mat, RealIndexer &&sie_mat, LambdaIndexer &&lambdas) const {
    // Create temporary arrays
    std::array<Real, nmat_> mass_fracs{};
    std::array<Real, nmat_> vol_fracs{};
    std::array<size_t, nmat_> pte_mats{};
    std::array<Real, nmat_> pres_mat{};
    std::array<Real, nmat_> temp_mat{};

    // Extract the mass fractions from the lambdas
    PopulateMassFracArray(mass_fracs, lambdas); // Guarantees sum = 1

    // Create array of pointers to EOS models
    const auto eos_arr = CreateEOSArray();

    // Compute volume fractions from mass fractions and densities. Mark
    // materials to participate in PTE
    size_t npte = 0;
    Real vfrac_tot = 0;
    sie_tot = 0; // Initialize to some value for bad guesses or non-participating
                 // materials
    SetUpPTE(mass_fracs, vol_fracs, density_mat, pte_mats, npte, vfrac_tot, density_tot,
             sie_tot, eos_arr);

    // Initialize solver scratch memory (including cache)
    constexpr int pte_solver_scratch_size =
        PTESolverFixedPRequiredScratch(nmat_) + nmat_ * MAX_NUM_LAMBDAS;
    std::array<Real, pte_solver_scratch_size> solver_scratch{};

    // Get cache from offsets into scratch
    const int neq = npte + 1;
    singularity::mix_impl::CacheAccessor cache(solver_scratch.data() + neq * (neq + 4) +
                                               2 * npte);

    // Create indexers for the EOS
    const auto eos_idxr = GenericIndexer(eos_arr, pte_mats);

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
      PTESolverFixedP<decltype(eos_idxr), decltype(density_idxr), decltype(cache)> method(
          npte, eos_idxr, vfrac_tot, pressure, density_idxr, vfrac_idxr, sie_idxr,
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
      RootFinding1D::Status root_status =
          RootFinding1D::findRoot(p_from_t, pressure, r_temp, 1.e-10 * r_temp,
                                  1.e10 * r_temp, 1.e-12, 1.e-12, temperature_root);

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

    // Loop through non-participating materials and ensure their energy and
    // density are consistent with P-T state. Otherwise, the material energy for
    // small-mass-fraction materials would be zero and this would change the
    // energy sum below. If a material's minimum pressure is above the PTE
    // pressure, then we keep the zero energy and cell density
    if (npte != nmat_) {
      for (size_t m = 0; m < nmat_; m++) {
        // Zero volume fraction materials didn't participate in PTE
        if (vol_fracs[m] > 0 || eos_arr[m].MinimumPressure() > pressure) {
          continue;
        }
        eos_arr[m].DensityEnergyFromPressureTemperature(
            pressure, temperature, std::forward<LambdaIndexer>(lambdas), density_mat[m],
            sie_mat[m]);
      }
    }

    // Recalculate sie_tot as the sum of all material energies
    sie_tot = 0;
    for (size_t m = 0; m < nmat_; m++) {
      sie_tot += sie_mat[m] * mass_fracs[m];
    }

    return status;
  }

  template <typename RealIndexer, typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_INLINE_FUNCTION void
  GetStatesFromPressureTemperature(Real &density_tot, Real &sie_tot, const Real pressure,
                                   const Real temperature, RealIndexer &density_mat,
                                   RealIndexer &sie_mat, LambdaIndexer &lambdas) const {
    // Extract the mass fractions from the lambdas
    std::array<Real, nmat_> mass_fracs{};
    PopulateMassFracArray(mass_fracs, lambdas);

    // Call DensityEnergyFromPressureTemperature() on all EOS
    callDensityEnergyFromPressureTemperatureAll(pressure, temperature, lambdas,
                                                density_mat, sie_mat,
                                                std::make_index_sequence<nmat_>{});

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
        [this](auto &&...eos_models) { return std::array<eosT_, nmat_>{&eos_models...}; },
        models_);
  }

  template <typename RealIndexer, typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_INLINE_FUNCTION void PopulateMassFracArray(RealIndexer &mass_fracs,
                                                      LambdaIndexer &lambdas) const {

    // Wrap the implementation function so we don't have index sequences elsewhere
    assign_mass_fractions(mass_fracs, lambdas, std::make_index_sequence<nmat_>{});

    // Normalize the mass fractions. JHP: Is this necessary?
    Real mass_frac_sum = math_utils::sum_neumaier(mass_fracs, nmat_);
    PORTABLE_REQUIRE(mass_frac_sum > 0., "Mass fractions sum to zero");
    for (size_t m = 0; m < nmat_; m++) {
      mass_fracs[m] /= mass_frac_sum;
    }
  }

  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie, LambdaIndexer &&lambda) const {

    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, temperature;
    SolverStatus status = GetStatesFromDensityEnergy(rho, sie, pressure, temperature,
                                                     density_mat, sie_mat, lambda);

    return temperature;
  }
  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_FUNCTION Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                                           LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, temperature;
    SolverStatus status = GetStatesFromDensityEnergy(rho, sie, pressure, temperature,
                                                     density_mat, sie_mat, lambda);

    return pressure;
  }
  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_FUNCTION Real MinInternalEnergyFromDensity(const Real rho,
                                                      LambdaIndexer &&lambda) const {
    // Technically this question is ill-defined since we need another thermodynamic
    // variable to constrain the volume fractions of the individual components
    // that combine to give the desired density. Instead we can interpret this
    // function as returning the energy along the lowest isotherm of an EOS.
    const auto min_isotherm = MinimumTemperature();
    return InternalEnergyFromDensityTemperature(rho, min_isotherm, lambda);
  }
  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_FUNCTION Real EntropyFromDensityInternalEnergy(const Real rho, const Real sie,
                                                          LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, temperature;
    SolverStatus status = GetStatesFromDensityEnergy(rho, sie, pressure, temperature,
                                                     density_mat, sie_mat, lambda);
    return massFracAverageQuantityAtManyStates(
        [](auto const &eos, Real const rho, Real const sie, Real const temperature,
           LambdaIndexer lambda) {
          if constexpr (eos.PreferredInput() ==
                        (thermalqs::density | thermalqs::temperature)) {
            return eos.EntropyFromDensityTemperature(rho, temperature, lambda);
          } else {
            return eos.EntropyFromDensityInternalEnergy(rho, sie, lambda);
          }
        },
        density_mat, sie_mat, temperature, lambda);
  }
  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, temperature;
    SolverStatus status =
        GetStatesFromDensityEnergy(rho, sie, pressure, temperature, density_mat, sie_mat);

    return massFracAverageQuantityAtManyStates(
        [](auto const &eos, Real const rho, Real const sie, Real const temperature,
           LambdaIndexer &&lambda) {
          if constexpr (eos.PreferredInput() ==
                        (thermalqs::density | thermalqs::temperature)) {
            return eos.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
          } else {
            return eos.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
          }
        },
        density_mat, sie_mat, temperature, lambda);
  }
  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, temperature;
    SolverStatus status = GetStatesFromDensityEnergy(rho, sie, pressure, temperature,
                                                     density_mat, sie_mat, lambda);

    BulkModAvgT avg_funct{};
    return avg_funct(
        [](auto const &eos, Real const rho, Real const sie, Real const temperature,
           LambdaIndexer &&lambda) {
          if constexpr (eos.PreferredInput() ==
                        (thermalqs::density | thermalqs::temperature)) {
            return eos.BulkModulusFromDensityTemperature(rho, temperature, lambda);
          } else {
            return eos.BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
          }
        },
        density_mat, sie_mat, temperature, lambda, std::make_index_sequence<nmat_>{});
  }
  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, temperature;
    SolverStatus status = GetStatesFromDensityEnergy(rho, sie, pressure, temperature,
                                                     density_mat, sie_mat, lambda);

    GruneisenAvgT avg_funct{};
    return avg_funct(
        [](auto const &eos, Real const rho, Real const sie, Real const temperature,
           LambdaIndexer &&lambda) {
          if constexpr (eos.PreferredInput() ==
                        (thermalqs::density | thermalqs::temperature)) {
            return eos.GruneisenParamFromDensityTemperature(rho, temperature, lambda);
          } else {
            return eos.GruneisenParamFromDensityInternalEnergy(rho, sie, lambda);
          }
        },
        density_mat, sie_mat, temperature, lambda, std::make_index_sequence<nmat_>{});
  }
  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, sie;
    SolverStatus status = GetStatesFromDensityTemperature(rho, sie, pressure, temperature,
                                                          density_mat, sie_mat, lambda);

    return sie;
  }
  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_FUNCTION Real PressureFromDensityTemperature(const Real rho,
                                                        const Real temperature,
                                                        LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, sie;
    SolverStatus status = GetStatesFromDensityTemperature(rho, sie, pressure, temperature,
                                                          density_mat, sie_mat, lambda);

    return pressure;
  }
  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_FUNCTION Real EntropyFromDensityTemperature(const Real rho,
                                                       const Real temperature,
                                                       LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, sie;
    SolverStatus status = GetStatesFromDensityTemperature(rho, sie, pressure, temperature,
                                                          density_mat, sie_mat, lambda);
    return massFracAverageQuantityAtManyStates(
        [](auto const &eos, Real const rho, Real const sie, Real const temperature,
           LambdaIndexer lambda) {
          if constexpr (eos.PreferredInput() ==
                        (thermalqs::density | thermalqs::temperature)) {
            return eos.EntropyFromDensityTemperature(rho, temperature, lambda);
          } else {
            return eos.EntropyFromDensityInternalEnergy(rho, sie, lambda);
          }
        },
        density_mat, sie_mat, temperature, lambda);
  }
  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, sie;
    SolverStatus status = GetStatesFromDensityTemperature(rho, sie, pressure, temperature,
                                                          density_mat, sie_mat, lambda);
    return massFracAverageQuantityAtManyStates(
        [](auto const &eos, Real const rho, Real const sie, Real const temperature,
           LambdaIndexer lambda) {
          if constexpr (eos.PreferredInput() ==
                        (thermalqs::density | thermalqs::temperature)) {
            return eos.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
          } else {
            return eos.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
          }
        },
        density_mat, sie_mat, temperature, lambda);
  }
  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_FUNCTION Real BulkModulusFromDensityTemperature(const Real rho,
                                                           const Real temperature,
                                                           LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, sie;
    SolverStatus status = GetStatesFromDensityTemperature(rho, sie, pressure, temperature,
                                                          density_mat, sie_mat, lambda);
    BulkModAvgT avg_funct{};
    return avg_funct(
        [](auto const &eos, Real const rho, Real const sie, Real const temperature,
           LambdaIndexer lambda) {
          if constexpr (eos.PreferredInput() ==
                        (thermalqs::density | thermalqs::temperature)) {
            return eos.BulkModulusFromDensityTemperature(rho, temperature, lambda);
          } else {
            return eos.BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
          }
        },
        density_mat, sie_mat, temperature, lambda, std::make_index_sequence<nmat_>{});
  }
  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature, LambdaIndexer &&lambda) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    Real pressure, sie;
    SolverStatus status = GetStatesFromDensityTemperature(rho, sie, pressure, temperature,
                                                          density_mat, sie_mat, lambda);
    GruneisenAvgT avg_funct{};
    return avg_funct(
        [](auto const &eos, Real const rho, Real const sie, Real const temperature,
           LambdaIndexer lambda) {
          if constexpr (eos.PreferredInput() ==
                        (thermalqs::density | thermalqs::temperature)) {
            return eos.GruneisenParamFromDensityTemperature(rho, temperature, lambda);
          } else {
            return eos.GruneisenParamFromDensityInternalEnergy(rho, sie, lambda);
          }
        },
        density_mat, sie_mat, temperature, lambda, std::make_index_sequence<nmat_>{});
  }
  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
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
      GetStatesFromPressureTemperature(energy, rho, press, temp, density_mat, sie_mat,
                                       lambda);
    } else if (input & thermalqs::density && input & thermalqs::temperature) {
      PORTABLE_REQUIRE(
          output & thermalqs::pressure && output & thermalqs::specific_internal_energy,
          "Internal energy and pressure must be defined as outputs when density and "
          "temperature are inputs");
      GetStatesFromDensityTemperature(energy, rho, press, temp, density_mat, sie_mat,
                                      lambda);
    } else if (input & thermalqs::density &&
               input & thermalqs::specific_internal_energy) {
      PORTABLE_REQUIRE(
          output & thermalqs::pressure && output & thermalqs::specific_internal_energy,
          "Temperature and pressure must be defined as outputs when density and "
          "internal energy are inputs");
      GetStatesFromDensityEnergy(energy, rho, press, temp, density_mat, sie_mat, lambda);
    } else if (input & thermalqs::density && input & thermalqs::pressure) {
      PORTABLE_REQUIRE(
          output & thermalqs::pressure && output & thermalqs::specific_internal_energy,
          "Temperature and pressure must be defined as outputs when density and "
          "internal energy are inputs");
      GetStatesFromDensityPressure(energy, rho, press, temp, density_mat, sie_mat,
                                   lambda);
    }

    if (output & thermalqs::specific_heat) {
      cv = massFracAverageQuantityAtManyStates(
          [](auto const &eos, Real const rho, Real const sie, Real const temperature,
             LambdaIndexer lambda) {
            if constexpr (eos.PreferredInput() ==
                          (thermalqs::density | thermalqs::temperature)) {
              return eos.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
            } else {
              return eos.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
            }
          },
          density_mat, sie_mat, temp, lambda);
    }

    if (output & thermalqs::bulk_modulus) {
      BulkModAvgT avg_funct{};
      bmod = avg_funct(
          [](auto const &eos, Real const rho, Real const sie, Real const temperature,
             LambdaIndexer lambda) {
            if constexpr (eos.PreferredInput() ==
                          (thermalqs::density | thermalqs::temperature)) {
              return eos.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
            } else {
              return eos.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
            }
          },
          density_mat, sie_mat, temp, lambda, std::make_index_sequence<nmat_>{});
    }
  }
  constexpr static inline int nlambda() noexcept {
    // Sum the lambda requirements from each EOS with the number of mass
    // fractions required.

    // TODO: Really some EOS will share lambda variables and we probably want
    // the superset of all possible lambda variables
    return (EOSModelsT::nlambda() + ...) + nmat_;
  }

  static constexpr unsigned long PreferredInput() {
    // A prevoius version of this function had each EOS "vote" for the
    // preferred input.

    // Really though, this is a complex question from both a performance and
    // accuracy perspective. Any lookup other than P-T will trigger a PTE
    // solve, and among the PTE solves, density-energy currently uses the 2x2
    // P-T formulation, which is cheaper than density-temperature (the other
    // lookups). There could also be concerns about the underlying EOS and
    // their ability to accurately (and efficiently) produce a state from
    // pressure and temperature, but this also applies to the 2x2 solver since
    // it uses P-T lookups.

    // Given these considerations, it seems prudent to just assume P-T lookups
    // are "preferred"
    return (thermalqs::pressure | thermalqs::temperature);
  }

  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    // NOTE: we are currently statically allocating the memory for the PTE
    // solvers so we don't need to include this in the scratch memory
    // consideration. Also, we aren't explicitly able to pass scratch memory to
    // the individual EOS so their scratch requirements are null anyway

    // In the future we may want to provide dynamic pools via the lambda
    // parameter though, in which case we'll need this function to basically
    // provide the PTE solver scratch requirements with additional memory for
    // cache memory

    return 0;
  }

  static inline unsigned long max_scratch_size(unsigned int nelements) {
    // See scratch_size() for an explanation of why this is zero
    return 0;
  }

  PORTABLE_FUNCTION void PrintParams() const {
    // print the parameters for each EOS via fold expression
    std::apply([](const auto &...eos_models) { (eos_models.PrintParams(), ...); },
               models_);
  }

  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       LambdaIndexer &&lambda, Real &rho,
                                       Real &sie) const {
    std::array<Real, nmat_> density_mat{};
    std::array<Real, nmat_> sie_mat{};
    GetStatesFromPressureTemperature(sie, rho, press, temp, density_mat, sie_mat, lambda);
  }

  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
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
    GetStatesFromPressureTemperature(sie, rho, press, temp, density_mat, sie_mat, lambda);

    // Mass-fraction average the heat capacity
    cv = massFracAverageQuantityAtManyStates(
        [](auto const &eos, Real const rho, Real const sie, Real const temperature,
           LambdaIndexer lambda) {
          if constexpr (eos.PreferredInput() ==
                        (thermalqs::density | thermalqs::temperature)) {
            return eos.SpecificHeatFromDensityTemperature(rho, temperature, lambda);
          } else {
            return eos.SpecificHeatFromDensityInternalEnergy(rho, sie, lambda);
          }
        },
        density_mat, sie_mat, temp, lambda);

    // Use injected average functions for gruneisen and bulk modulus averaging
    GruneisenAvgT g_avg_funct{};
    auto avg_gruneisen = g_avg_funct(
        [](auto const &eos, Real const rho, Real const sie, Real const temperature,
           LambdaIndexer lambda) {
          if constexpr (eos.PreferredInput() ==
                        (thermalqs::density | thermalqs::temperature)) {
            return eos.GruneisenParamFromDensityTemperature(rho, temperature, lambda);
          } else {
            return eos.GruneisenParamFromDensityInternalEnergy(rho, sie, lambda);
          }
        },
        density_mat, sie_mat, temp, lambda, std::make_index_sequence<nmat_>{});
    dpde = rho * avg_gruneisen;

    BulkModAvgT b_avg_funct{};
    bmod = b_avg_funct(
        [](auto const &eos, Real const rho, Real const sie, Real const temperature,
           LambdaIndexer lambda) {
          if constexpr (eos.PreferredInput() ==
                        (thermalqs::density | thermalqs::temperature)) {
            return eos.BulkModulusFromDensityTemperature(rho, temperature, lambda);
          } else {
            return eos.BulkModulusFromDensityInternalEnergy(rho, sie, lambda);
          }
        },
        density_mat, sie_mat, temp, lambda, std::make_index_sequence<nmat_>{});
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumDensity() const {
    return max_over_models([](auto const &eos) { return eos.MinimumDensity(); });
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumTemperature() const {
    return max_over_models([](auto const &eos) { return eos.MinimumTemperature(); });
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MaximumDensity() const {
    return min_over_models([](auto const &eos) { return eos.MaximumDensity(); });
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumPressure() const {
    return max_over_models([](auto const &eos) { return eos.MinimumPressure(); });
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MaximumPressureAtTemperature(const Real temp) const {
    return min_over_models([](auto const &eos, Real const temp) {
      return eos.MaximumPressureAtTemperature(temp);
    });
  }

  // 3T member functions
  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicMassFromDensityTemperature(
      const Real density, const Real temperature,
      LambdaIndexer &&lambda = static_cast<Real *>(nullptr)) const {
    return massFracAverageQuantityAtOneState(
        [&](auto const &eos, Real R, Real T, LambdaIndexer lambda) {
          return eos.MeanAtomicMassFromDensityTemperature(R, T, lambda);
        },
        density, temperature, lambda, std::make_index_sequence<nmat_>());
  }
  template <typename LambdaIndexer,
            SINGULARITY_INDEXER_HAS_MASS_FRAC(LambdaIndexer, nmat_)>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicNumberFromDensityTemperature(
      const Real density, const Real temperature,
      LambdaIndexer &&lambda = static_cast<Real *>(nullptr)) const {
    return massFracAverageQuantityAtOneState(
        [&](auto const &eos, Real R, Real T, LambdaIndexer lambda) {
          return eos.MeanAtomicNumberFromDensityTemperature(R, T, lambda);
        },
        density, temperature, lambda, std::make_index_sequence<nmat_>());
  }
  PORTABLE_INLINE_FUNCTION
  Real MeanAtomicMass() const {
    // Since we have no mass fraction information, assume equal weights in this
    // instance
    return (0. + ... + EOSModelsT::MeanAtomicMass()) / nmat_;
  }
  PORTABLE_INLINE_FUNCTION
  Real MeanAtomicNumber() const {
    // Since we have no mass fraction information, assume equal weights in this
    // instance
    return (0. + ... + EOSModelsT::MeanAtomicNumber()) / nmat_;
  }

  // Modifier member functions
  static inline constexpr bool IsModified() {
    // Even though this class could be considered a modifier, there is no clean
    // way to return an "unmodified object" since this class takes a collection
    // of EOS objects. So we'll claim this isn't a modifier for this purpose
    return false;
  }
  inline constexpr const auto UnmodifyOnce() {
    // Even though this class could be considered a modifier, there is no clean
    // way to return an "unmodified object" since this class takes a collection
    // of EOS objects. So we'll claim this isn't a modifier for this purpose
    return *this;
  }

  // Dynamic and shared memory member functions (normally included via macro)
  std::size_t DynamicMemorySizeInBytes() const {
    return (0 + ... + EOSModelsT::DynamicMemorySizeInBytes());
  }
  std::size_t SharedMemorySizeInBytes() const {
    return (0 + ... + EOSModelsT::SharedMemorySizeInBytes());
  }

  std::size_t DumpDynamicMemory(char *dst) {
    // Call each model's DumpDynamicMemory(dst_i) in order and return total
    // bytes written
    std::size_t total = 0;

    std::apply(
        [&](auto &...eos) {
          // Each term updates 'total' before the next term sees it.
          ((total += eos.DumpDynamicMemory(dst + total)), ...);
        },
        models_);

    return total;
  }

  std::size_t SetDynamicMemory(char *src, SharedMemSettings stngs = DEFAULT_SHMEM_STNGS) {

    // Source and shared memory pointers
    char *p_src = src;
    char *p_shared = stngs.data;

    // Work on a local, mutable copy in order to increment the shared memory pointer
    SharedMemSettings current_shared_stngs = stngs;

    // We keep track of how far we increment in data space
    std::size_t total = 0;

    std::apply(
        [&p_src, &p_shared, &current_shared_stngs, &total](auto &...eos) {
          std::size_t shared_increment = 0;
          std::size_t src_increment;
          // One iteration per EOS, left-to-right.
          (((
                // 1) Increment the per-EOS shared-memory pointer
                current_shared_stngs.data =
                    p_shared == nullptr ? (p_shared + shared_increment) : nullptr,

                // 2) Set the dynamic memory and record how much was used
                src_increment = eos.SetDynamicMemory(p_src, current_shared_stngs),

                // 3) Bump the shared-memory offset by the model's declared size.
                //    This _may_ be different from the  src_increment
                shared_increment = eos.SharedMemorySizeInBytes(),

                // 4) We need to figure out how much to incremennt the src pointer
                //    depending on the output of `AllDynamicMemoryIsShareable()`:
                //    - Spiner either points to shared memory or the src memory, but
                //      not both, so its shared memory size is the same as its
                //      dynamic memory size. So if shared memory is used, the src
                //      memory shouldn't be incremented. For spiner,
                //      `AllDynamicMemoryIsShareable()` returns true.
                //    - EOSPAC reads from the src memory and then either writes some
                //      of this data to the shared memory pointer on the root node
                //      or just points to the shared memory. Either way, both the
                //      shared memory and src memory pointers need to be
                //      incremented. For EOSPAC, `AllDynamicMemoryIsShareable()`
                //      returns false.
                p_src = p_shared == nullptr || eos.AllDynamicMemoryIsShareable()
                            ? p_src
                            : p_src + src_increment,

                // 5) Our output should reflect how much we've read
                total += src_increment),
            void() // make the subexpression a valid comma-fold operand
            ),
           ...);
        },
        models_);

    return total;
  }

  constexpr bool AllDynamicMemoryIsShareable() const {
    return std::apply(
        [](auto const &...eos) {
          // AND operation across all individual models. Returns true IFF all models
          // have sharable dynamic memory. If this isn't true, then the user needs
          // to keep track of SharedMemorySizeInBytes and DynamicMemorySizeInBytes
          // separately
          return (eos.AllDynamicMemoryIsShareable() && ...);
        },
        models_);
  }
};

// Factory function helper to group EOS together and provide default averaging
// behavior
template <typename BulkModAvgT = VolumeFracHarmonicAverageFunctor,
          typename GruneisenAvgT = VolumeFracHarmonicAverageFunctor,
          typename... EOSModelsT>
inline auto make_MultiEOS(EOSModelsT &&...eos_models) {
  return MultiEOS<std::remove_cv_t<std::remove_reference_t<BulkModAvgT>>,
                  std::remove_cv_t<std::remove_reference_t<GruneisenAvgT>>,
                  std::remove_cv_t<std::remove_reference_t<EOSModelsT>>...>(
      std::forward<EOSModelsT>(eos_models)...);
}

// Factory function helper to group EOS together and provide default averaging
// behavior but including the mass fraction cutoff as a constructor argument
template <typename BulkModAvgT = VolumeFracHarmonicAverageFunctor,
          typename GruneisenAvgT = VolumeFracHarmonicAverageFunctor,
          typename... EOSModelsT>
inline auto make_MultiEOS(Real mass_frac_cutoff_in, EOSModelsT &&...eos_models) {
  return MultiEOS<std::remove_cv_t<std::remove_reference_t<BulkModAvgT>>,
                  std::remove_cv_t<std::remove_reference_t<GruneisenAvgT>>,
                  std::remove_cv_t<std::remove_reference_t<EOSModelsT>>...>(
      mass_frac_cutoff_in, std::forward<EOSModelsT>(eos_models)...);
}

} // namespace singularity

#endif // #ifndef _SINGULARITY_EOS_EOS_MULTI_EOS_
