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

#ifndef _SINGULARITY_EOS_EOS_EOS_SPINER_HPP_
#define _SINGULARITY_EOS_EOS_EOS_SPINER_HPP_
#include <type_traits>

#ifdef SINGULARITY_USE_SPINER_WITH_HDF5
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
// #include <iostream> // debug
// #include <stdio.h> // debug

#include <hdf5.h>
#include <hdf5_hl.h>

// ports-of-call
#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>

// base
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/fast-math/logs.hpp>
#include <singularity-eos/base/indexable_types.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>
#include <singularity-eos/base/sp5/singularity_eos_sp5.hpp>
#include <singularity-eos/base/spiner_table_utils.hpp>
#include <singularity-eos/base/variadic_utils.hpp>
#include <singularity-eos/eos/eos_base.hpp>

// spiner
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/sp5.hpp>
#include <spiner/spiner_types.hpp>

#define SPINER_EOS_VERBOSE (0)
#define ROOT_FINDER (RootFinding1D::regula_falsi)

namespace singularity {

using namespace eos_base;

/*
  Tables all have indep. variables log10(rho), log10(T)

  Extrapolation strategy:
  ------------------------------------------------------------
  We use cold curve for low temperatures/energies
  and assume ideal gas for large temperatures/energies.

  For low densities, we floor the density. For high densities, we
  we use log-linear extrapolation.
*/
class SpinerEOSDependsRhoT : public EosBase<SpinerEOSDependsRhoT> {
  friend class table_utils::SpinerTricks<SpinerEOSDependsRhoT>;
  using SpinerTricks = table_utils::SpinerTricks<SpinerEOSDependsRhoT>;

 public:
  static constexpr int NGRIDS = 3;
  using Grid_t = Spiner::PiecewiseGrid1D<Real, NGRIDS>;
  using DataBox = Spiner::DataBox<Real, Grid_t>;

  // A weakly typed index map for lambdas
  struct Lambda {
    enum Index { lRho = 0, lT = 1 };
  };
  SG_ADD_DEFAULT_MEAN_ATOMIC_FUNCTIONS(AZbar_)
  SG_ADD_BASE_CLASS_USINGS(SpinerEOSDependsRhoT);
  inline SpinerEOSDependsRhoT(const std::string &filename, int matid, TableSplit split,
                              bool reproducibility_mode = false);
  inline SpinerEOSDependsRhoT(const std::string &filename, int matid,
                              bool reproducibility_mode = false)
      : SpinerEOSDependsRhoT(filename, matid, TableSplit::Total, reproducibility_mode) {}
  inline SpinerEOSDependsRhoT(const std::string &filename,
                              const std::string &materialName, TableSplit split,
                              bool reproducibility_mode = false);
  inline SpinerEOSDependsRhoT(const std::string &filename,
                              const std::string &materialName,
                              bool reproducibility_mode = false)
      : SpinerEOSDependsRhoT(filename, materialName, TableSplit::Total,
                             reproducibility_mode) {}
  PORTABLE_INLINE_FUNCTION
  SpinerEOSDependsRhoT()
      : memoryStatus_(DataStatus::Deallocated), split_(TableSplit::Total) {}

  inline SpinerEOSDependsRhoT GetOnDevice();

  PORTABLE_INLINE_FUNCTION void CheckParams() const {
    PORTABLE_ALWAYS_REQUIRE(numRho_ > 0, "Finite number of density points");
    PORTABLE_ALWAYS_REQUIRE(numT_ > 0, "Finite number of temperature points");
    PORTABLE_ALWAYS_REQUIRE(!(std::isnan(lRhoMin_) || std::isnan(lRhoMax_)),
                            "Density bounds well defined");
    PORTABLE_ALWAYS_REQUIRE(lRhoMax_ > lRhoMin_, "Density bounds ordered");
    PORTABLE_ALWAYS_REQUIRE(rhoMax_ > 0, "Max density must be positive");
    PORTABLE_ALWAYS_REQUIRE(!(std::isnan(lTMin_) || std::isnan(lTMax_)),
                            "Temperature bounds well defined");
    PORTABLE_ALWAYS_REQUIRE(lTMax_ > lTMin_, "Temperature bounds ordered");
    PORTABLE_ALWAYS_REQUIRE(TMax_ > 0, "Max temperature must be positive");
  }

  std::size_t DynamicMemorySizeInBytes() const;
  std::size_t DumpDynamicMemory(char *dst);
  std::size_t SetDynamicMemory(char *src,
                               const SharedMemSettings &stngs = DEFAULT_SHMEM_STNGS);

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MinInternalEnergyFromDensity(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  PressureFromDensityTemperature(const Real rho, const Real temperature,
                                 Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
          const unsigned long output,
          Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                         Real &bmod, Real &dpde, Real &dvdt,
                         Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;

  PORTABLE_INLINE_FUNCTION
  Real RhoPmin(const Real temp) const;

  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return 0;
  }
  static inline unsigned long max_scratch_size(unsigned int nelements) { return 0; }
  int matid() const { return matid_; }
  PORTABLE_FORCEINLINE_FUNCTION Real lRhoOffset() const { return lRhoOffset_; }
  PORTABLE_FORCEINLINE_FUNCTION Real lTOffset() const { return lTOffset_; }
  PORTABLE_FORCEINLINE_FUNCTION Real rhoMin() const { return rho_(lRhoMin_); }
  PORTABLE_FORCEINLINE_FUNCTION Real rhoMax() const { return rhoMax_; }
  PORTABLE_FORCEINLINE_FUNCTION Real TMin() const { return T_(lTMin_); }
  PORTABLE_FORCEINLINE_FUNCTION Real TMax() const { return TMax_; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    static constexpr char s1[]{"SpinerEOS Parameters:"};
    static constexpr char s2[]{"depends on log_10(rho) and log_10(temp)"};
    static constexpr char s3[]{"EOS file   = "};
    static constexpr char s4[]{"EOS mat ID = "};
    static constexpr char s5[]{"EOS name   = "};
    printf("%s\n\t%s\n\t%s\n\t%s%i\n\t%s\n", s1, s2, s3, s4, matid_, s5);
    return;
  }
  PORTABLE_FORCEINLINE_FUNCTION Real MinimumDensity() const { return rhoMin(); }
  PORTABLE_FORCEINLINE_FUNCTION Real MinimumTemperature() const { return T_(lTMin_); }
  PORTABLE_FORCEINLINE_FUNCTION Real MaximumDensity() const { return rhoMax(); }
  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumPressure() const { return PMin_; }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return _n_lambda; }
  template <typename T>
  static inline constexpr bool NeedsLambda() {
    using namespace IndexableTypes;
    return std::is_same<T, LogDensity>::value || std::is_same<T, LogTemperature>::value;
  }
  RootFinding1D::RootCounts counts;
  inline void Finalize();
  static std::string EosType() { return std::string("SpinerEOSDependsRhoT"); }
  static std::string EosPyType() { return EosType(); }

 private:
  herr_t loadDataboxes_(const std::string &matid_str, hid_t file, hid_t lTGroup,
                        hid_t coldGroup);
  inline void fixBulkModulus_();
  inline void setlTColdCrit_();

  static PORTABLE_FORCEINLINE_FUNCTION Real toLog_(const Real x, const Real offset) {
    // return std::log10(x + offset + robust::EPS());
    // return std::log10(std::abs(std::max(x,-offset) + offset)+robust::EPS());
    return FastMath::log10(std::abs(std::max(x, -offset) + offset) + robust::EPS());
  }
  static PORTABLE_FORCEINLINE_FUNCTION Real fromLog_(const Real lx, const Real offset) {
    return FastMath::pow10(lx) - offset;
  }
  PORTABLE_FORCEINLINE_FUNCTION
  Real lRho_(const Real rho) const noexcept {
    Real out = toLog_(rho, lRhoOffset_);
    return out;
    // return out < lRhoMin_ ? lRhoMin_ : out;
  }
  PORTABLE_FORCEINLINE_FUNCTION
  Real lT_(const Real T) const noexcept { return toLog_(T, lTOffset_); }
  PORTABLE_FORCEINLINE_FUNCTION
  Real rho_(const Real lRho) const noexcept {
    Real rho = fromLog_(lRho, lRhoOffset_);
    return rho < 0 ? 0 : rho;
  }
  PORTABLE_FORCEINLINE_FUNCTION
  Real T_(const Real lT) const noexcept { return fromLog_(lT, lTOffset_); }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  lTFromlRhoSie_(const Real lRho, const Real sie, TableStatus &whereAmI,
                 Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  lTFromlRhoP_(const Real lRho, const Real press, TableStatus &whereAmI,
               Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  lRhoFromPlT_(const Real P, const Real lT, TableStatus &whereAmI,
               Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  getLogsRhoT_(const Real rho, const Real temperature, Real &lRho, Real &lT,
               Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  PORTABLE_INLINE_FUNCTION
  Real sieFromlRhoTlT_(const Real lRho, const Real T, const Real lT,
                       const TableStatus &whereAmI) const;
  PORTABLE_INLINE_FUNCTION
  Real PFromRholRhoTlT_(const Real rho, const Real lRho, const Real T, const Real lT,
                        const TableStatus &whereAmI) const;
  PORTABLE_INLINE_FUNCTION
  Real CvFromlRholT_(const Real lRho, const Real lT, const TableStatus &whereAmI) const;
  PORTABLE_INLINE_FUNCTION
  Real bModFromRholRhoTlT_(const Real rho, const Real lRho, const Real T, const Real lT,
                           const TableStatus &whereAmI) const;
  PORTABLE_INLINE_FUNCTION
  TableStatus getLocDependsRhoSie_(const Real lRho, const Real sie) const;
  PORTABLE_INLINE_FUNCTION
  TableStatus getLocDependsRhoT_(const Real lRho, const Real lT) const;
  // PORTABLE_INLINE_FUNCTION
  // Real sieToColdInterval_(const Real lRho, const Real sie) const;

  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::temperature;
  // static constexpr const char _eos_type[] {"SpinerEOSDependsRhoT"};

  // TODO(JMM): Could unify declarations and macro below by using
  // reference_wrapper instead of pointers... worth it?
  DataBox P_, sie_, bMod_, dPdRho_, dPdE_, dTdRho_, dTdE_, dEdRho_, dEdT_;
  DataBox PMax_, sielTMax_, dEdTMax_, gm1Max_;
  DataBox lTColdCrit_;
  DataBox PCold_, sieCold_, bModCold_;
  DataBox dPdRhoCold_, dPdECold_, dTdRhoCold_, dTdECold_, dEdTCold_;
  DataBox rho_at_pmin_;

  // TODO(JMM): Pointers here? or reference_wrapper? IMO the pointers are more clear
#define DBLIST                                                                           \
  &P_, &sie_, &bMod_, &dPdRho_, &dPdE_, &dTdRho_, &dTdE_, &dEdRho_, &dEdT_, &PMax_,      \
      &sielTMax_, &dEdTMax_, &gm1Max_, &lTColdCrit_, &PCold_, &sieCold_, &bModCold_,     \
      &dPdRhoCold_, &dPdECold_, &dTdRhoCold_, &dTdECold_, &dEdTCold_, &rho_at_pmin_
  auto GetDataBoxPointers_() const { return std::vector<const DataBox *>{DBLIST}; }
  auto GetDataBoxPointers_() { return std::vector<DataBox *>{DBLIST}; }
#undef DBLIST

  int numRho_, numT_;
  Real lRhoMin_, lRhoMax_, rhoMax_;
  Real lRhoMinSearch_;
  Real lTMin_, lTMax_, TMax_;
  Real PMin_;
  Real rhoNormal_, TNormal_, sieNormal_, PNormal_;
  Real CvNormal_, bModNormal_, dPdENormal_, dVdTNormal_;
  Real lRhoOffset_, lTOffset_; // offsets must be non-negative
  MeanAtomicProperties AZbar_;
  int matid_;
  TableSplit split_;
  bool reproducible_;
  static constexpr const Real ROOT_THRESH = 1e-14; // TODO: experiment
  static constexpr const Real SOFT_THRESH = 1e-8;
  DataStatus memoryStatus_ = DataStatus::Deallocated;
  static constexpr const int _n_lambda = 2;
  static constexpr const char *_lambda_names[2] = {"log(rho)", "log(T)"};
};

/*
  TODO(JMM): Extrapolation Strategy
  ----------------------------------
  Currently the bottom of the table is the bound.
  We do constant extrapolation off the bottom of the table
  and linear extrapolation off the top.

  Since the table is in log-log space, this means extrapolation off
  the top is a power law. Extrapolation off the bottom is constant.
  This worked for nubhlight. But it may or may not work here. We will
  potentially need to revisit this.

  The best solution might be a three-part EOS matched to the bottom of
  the table, containing:
  - A polytropic term
  - A photon pressure term
  - An ideal gas term
  mitigated by Ye and (1-Ye) to control how important each term is.
 */
class SpinerEOSDependsRhoSie : public EosBase<SpinerEOSDependsRhoSie> {
  friend class table_utils::SpinerTricks<SpinerEOSDependsRhoSie>;

 public:
  struct Lambda {
    enum Index { lRho = 0 };
  };
  using Grid_t = SpinerEOSDependsRhoT::Grid_t;
  using DataBox = SpinerEOSDependsRhoT::DataBox;
  struct SP5Tables {
    DataBox P, bMod, dPdRho, dPdE, dTdRho, dTdE, dEdRho;
  };
  using STricks = table_utils::SpinerTricks<SpinerEOSDependsRhoSie>;

  SG_ADD_DEFAULT_MEAN_ATOMIC_FUNCTIONS(AZbar_)
  SG_ADD_BASE_CLASS_USINGS(SpinerEOSDependsRhoSie);
  PORTABLE_INLINE_FUNCTION SpinerEOSDependsRhoSie()
      : memoryStatus_(DataStatus::Deallocated) {}
  inline SpinerEOSDependsRhoSie(const std::string &filename, int matid, TableSplit split,
                                bool reproducibility_mode = false);
  inline SpinerEOSDependsRhoSie(const std::string &filename, int matid,
                                bool reproducibility_mode = false)
      : SpinerEOSDependsRhoSie(filename, matid, TableSplit::Total, reproducibility_mode) {
  }
  inline SpinerEOSDependsRhoSie(const std::string &filename,
                                const std::string &materialName, TableSplit split,
                                bool reproducibility_mode = false);
  inline SpinerEOSDependsRhoSie(const std::string &filename,
                                const std::string &materialName,
                                bool reproducibility_mode = false)
      : SpinerEOSDependsRhoSie(filename, materialName, TableSplit::Total,
                               reproducibility_mode) {}
  inline SpinerEOSDependsRhoSie GetOnDevice();

  PORTABLE_INLINE_FUNCTION void CheckParams() const {
    PORTABLE_ALWAYS_REQUIRE(numRho_ > 0, "Finite number of density points");
    PORTABLE_ALWAYS_REQUIRE(!(std::isnan(lRhoMin_) || std::isnan(lRhoMax_)),
                            "Density bounds well defined");
    PORTABLE_ALWAYS_REQUIRE(lRhoMax_ > lRhoMin_, "Density bounds ordered");
    PORTABLE_ALWAYS_REQUIRE(rhoMax_ > 0, "Max density must be positive");
  }

  std::size_t DynamicMemorySizeInBytes() const;
  std::size_t DumpDynamicMemory(char *dst);
  std::size_t SetDynamicMemory(char *src,
                               const SharedMemSettings &stngs = DEFAULT_SHMEM_STNGS);

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real T,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  PressureFromDensityTemperature(const Real rho, const Real T,
                                 Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MinInternalEnergyFromDensity(
      const Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real T,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real T,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real T,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
          const unsigned long output,
          Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                         Real &bmod, Real &dpde, Real &dvdt,
                         Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;

  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return 0;
  }
  static inline unsigned long max_scratch_size(unsigned int nelements) { return 0; }
  int matid() const { return matid_; }
  PORTABLE_FORCEINLINE_FUNCTION Real lRhoOffset() const { return lRhoOffset_; }
  PORTABLE_FORCEINLINE_FUNCTION Real lTOffset() const { return lTOffset_; }
  PORTABLE_FORCEINLINE_FUNCTION Real lEOffset() const { return lEOffset_; }
  PORTABLE_FORCEINLINE_FUNCTION Real rhoMin() const {
    return fromLog_(lRhoMin_, lRhoOffset_);
  }
  PORTABLE_FORCEINLINE_FUNCTION Real rhoMax() const { return rhoMax_; }
  PORTABLE_FORCEINLINE_FUNCTION Real TMin() const {
    return fromLog_(T_.range(0).min(), lTOffset_);
  }
  PORTABLE_FORCEINLINE_FUNCTION Real TMax() const {
    return fromLog_(T_.range(0).max(), lTOffset_);
  }
  PORTABLE_FORCEINLINE_FUNCTION Real sieMin() const {
    return fromLog_(sie_.range(0).min(), lEOffset_);
  }
  PORTABLE_FORCEINLINE_FUNCTION Real sieMax() const {
    return fromLog_(sie_.range(0).max(), lEOffset_);
  }

  PORTABLE_FORCEINLINE_FUNCTION Real MinimumDensity() const { return rhoMin(); }
  PORTABLE_FORCEINLINE_FUNCTION Real MinimumTemperature() const { return TMin(); }
  PORTABLE_FORCEINLINE_FUNCTION Real MaximumDensity() const { return rhoMax(); }
  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumPressure() const { return PMin_; }
  PORTABLE_INLINE_FUNCTION
  Real RhoPmin(const Real temp) const {
    return rho_at_pmin_.interpToReal(toLog_(temp, lTOffset_));
  }

  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return _n_lambda; }
  template <typename T>
  static inline constexpr bool NeedsLambda() {
    return std::is_same<T, IndexableTypes::LogDensity>::value;
  }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    static constexpr char s1[]{"SpinerEOS Parameters:"};
    static constexpr char s2[]{"depends on log_10(rho) and log_10(sie)"};
    static constexpr char s3[]{"EOS mat ID = "};
    printf("%s\n\t%s\n\t%s%i\n", s1, s2, s3, matid_);
    return;
  }
  RootFinding1D::RootCounts counts;
  static std::string EosType() { return std::string("SpinerEOSDependsRhoSie"); }
  static std::string EosPyType() { return EosType(); }
  inline void Finalize();

 private:
  inline herr_t loadDataboxes_(const std::string &matid_str, hid_t file, hid_t lTGroup,
                               hid_t lEGroup);
  inline void calcBMod_(SP5Tables &tables);

  static PORTABLE_FORCEINLINE_FUNCTION Real toLog_(const Real x, const Real offset) {
    // return std::log10(std::abs(std::max(x,-offset) + offset)+robust::EPS());
    return FastMath::log10(std::abs(std::max(x, -offset) + offset) + robust::EPS());
  }
  static PORTABLE_FORCEINLINE_FUNCTION Real fromLog_(const Real lx, const Real offset) {
    return FastMath::pow10(lx) - offset;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  interpRhoT_(const Real rho, const Real T, const DataBox &db,
              Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  interpRhoSie_(const Real rho, const Real sie, const DataBox &db,
                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real lRhoFromPlT_(const Real P, const Real lT,
                                             Indexer_t &&lambda) const;

  DataBox sie_; // depends on (rho,T)
  DataBox T_;   // depends on (rho, sie)
  DataBox rho_at_pmin_;
  SP5Tables dependsRhoT_;
  SP5Tables dependsRhoSie_;
  int numRho_, numT_;
  Real rhoNormal_, TNormal_, sieNormal_, PNormal_;
  Real CvNormal_, bModNormal_, dPdENormal_, dVdTNormal_;
  Real lRhoMin_, lRhoMax_, rhoMax_;
  DataBox PlRhoMax_, dPdRhoMax_;
  Real PMin_;
  Real lRhoOffset_, lTOffset_, lEOffset_; // offsets must be non-negative

#define DBLIST                                                                           \
  &sie_, &T_, &rho_at_pmin_, &(dependsRhoT_.P), &(dependsRhoT_.bMod),                    \
      &(dependsRhoT_.dPdRho), &(dependsRhoT_.dPdE), &(dependsRhoT_.dTdRho),              \
      &(dependsRhoT_.dTdE), &(dependsRhoT_.dEdRho), &(dependsRhoSie_.P),                 \
      &(dependsRhoSie_.bMod), &(dependsRhoSie_.dPdRho), &(dependsRhoSie_.dPdE),          \
      &(dependsRhoSie_.dTdRho), &(dependsRhoSie_.dTdE), &(dependsRhoSie_.dEdRho),        \
      &PlRhoMax_, &dPdRhoMax_
  std::vector<const DataBox *> GetDataBoxPointers_() const {
    return std::vector<const DataBox *>{DBLIST};
  }
  std::vector<DataBox *> GetDataBoxPointers_() { return std::vector<DataBox *>{DBLIST}; }
#undef DBLIST

  static constexpr unsigned long _preferred_input =
      thermalqs::density | thermalqs::temperature;
  // static constexpr const char _eos_type[] = "SpinerEOSDependsRhoSie";
  int matid_;
  TableSplit split_;
  MeanAtomicProperties AZbar_;
  bool reproducible_;
  static constexpr const int _n_lambda = 1;
  static constexpr const char *_lambda_names[1] = {"log(rho)"};
  DataStatus memoryStatus_ = DataStatus::Deallocated;
};

// implementation details below
// ======================================================================

// TODO: we're using log-linear interpolation, not log-log
// this may be suboptimal. We may want a way to do some variables
// in log-log. In particular, it might be good to do:
// pressure, energy, and bulk modulus in log-log.
// ~JMM

// replace lambdas with callable
namespace callable_interp {
using DataBox = SpinerEOSDependsRhoT::DataBox;
class l_interp {
 private:
  const DataBox &field;
  const Real fixed;

 public:
  PORTABLE_INLINE_FUNCTION
  l_interp(const DataBox &field_, const Real fixed_) : field{field_}, fixed{fixed_} {}

  PORTABLE_INLINE_FUNCTION Real operator()(const Real x) const {
    return field.interpToReal(x, fixed);
  }
};

class r_interp {
 private:
  const DataBox &field;
  const Real fixed;

 public:
  PORTABLE_INLINE_FUNCTION
  r_interp(const DataBox &field_, const Real fixed_) : field{field_}, fixed{fixed_} {}

  PORTABLE_INLINE_FUNCTION Real operator()(const Real x) const {
    return field.interpToReal(fixed, x);
  }
};

class prod_interp_1d {
 private:
  const DataBox &field1, field2;
  const Real r;

 public:
  PORTABLE_INLINE_FUNCTION
  prod_interp_1d(const DataBox &field1_, const DataBox &field2_, const Real r_)
      : field1{field1_}, field2{field2_}, r{r_} {}
  PORTABLE_INLINE_FUNCTION Real operator()(const Real x) const {
    return field1.interpToReal(x) * field2.interpToReal(x) * r * x;
  }
};

class interp {
 private:
  const DataBox &field;

 public:
  PORTABLE_INLINE_FUNCTION
  interp(const DataBox &field_) : field(field_) {}
  PORTABLE_INLINE_FUNCTION Real operator()(const Real x) const {
    return field.interpToReal(x);
  }
};
} // namespace callable_interp

inline SpinerEOSDependsRhoT::SpinerEOSDependsRhoT(const std::string &filename, int matid,
                                                  TableSplit split,
                                                  bool reproducibility_mode)
    : matid_(matid), split_(split), reproducible_(reproducibility_mode),
      memoryStatus_(DataStatus::OnHost) {

  std::string matid_str = std::to_string(matid);
  hid_t file, matGroup, lTGroup, coldGroup;
  herr_t status = H5_SUCCESS;

  file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  int log_type = FastMath::LogType::NQT1;
  if (H5LTfind_attribute(file, SP5::logType)) {
    H5LTget_attribute_int(file, "/", SP5::logType, &log_type);
  }
  PORTABLE_ALWAYS_REQUIRE(
      log_type == FastMath::Settings::log_type,
      "Log mode used at runtime must be identical to the one used to generate the file!");

  matGroup = H5Gopen(file, matid_str.c_str(), H5P_DEFAULT);

  std::string lTGroupName = SP5::Depends::logRhoLogT;
  if (split == TableSplit::ElectronOnly) {
    lTGroupName += (std::string("/") + SP5::SubTable::electronOnly);
  } else if (split == TableSplit::IonCold) {
    lTGroupName += (std::string("/") + SP5::SubTable::ionCold);
  }
  lTGroup = H5Gopen(matGroup, lTGroupName.c_str(), H5P_DEFAULT);
  coldGroup = H5Gopen(matGroup, SP5::Depends::coldCurve, H5P_DEFAULT);

  status += loadDataboxes_(matid_str, file, lTGroup, coldGroup);

  status += H5Gclose(lTGroup);
  status += H5Gclose(coldGroup);
  status += H5Gclose(matGroup);
  status += H5Fclose(file);

  if (status != H5_SUCCESS) {
    EOS_ERROR("SpinerDependsRHoT: HDF5 error\n"); // TODO: make this better
  }

  CheckParams();
}

inline SpinerEOSDependsRhoT::SpinerEOSDependsRhoT(const std::string &filename,
                                                  const std::string &materialName,
                                                  TableSplit split,
                                                  bool reproducibility_mode)
    : split_(split), reproducible_(reproducibility_mode),
      memoryStatus_(DataStatus::OnHost) {

  std::string matid_str;
  hid_t file, matGroup, lTGroup, subGroup, coldGroup;
  herr_t status = H5_SUCCESS;

  file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  matGroup = H5Gopen(file, materialName.c_str(), H5P_DEFAULT);

  std::string lTGroupName = SP5::Depends::logRhoLogT;
  if (split == TableSplit::ElectronOnly) {
    lTGroupName += (std::string("/") + SP5::SubTable::electronOnly);
  } else if (split == TableSplit::IonCold) {
    lTGroupName += (std::string("/") + SP5::SubTable::ionCold);
  }
  lTGroup = H5Gopen(matGroup, lTGroupName.c_str(), H5P_DEFAULT);
  coldGroup = H5Gopen(matGroup, SP5::Depends::coldCurve, H5P_DEFAULT);

  status +=
      H5LTget_attribute_int(file, materialName.c_str(), SP5::Material::matid, &matid_);
  matid_str = std::to_string(matid_);

  status += loadDataboxes_(matid_str, file, lTGroup, coldGroup);

  status += H5Gclose(lTGroup);
  status += H5Gclose(coldGroup);
  status += H5Gclose(matGroup);
  status += H5Fclose(file);

  if (status != H5_SUCCESS) {
    EOS_ERROR("SpinerDependsRhoT: HDF5 error\n");
  }

  CheckParams();
}

inline SpinerEOSDependsRhoT SpinerEOSDependsRhoT::GetOnDevice() {
  return SpinerTricks::GetOnDevice(this);
}

void SpinerEOSDependsRhoT::Finalize() { SpinerTricks::Finalize(this); }

inline std::size_t SpinerEOSDependsRhoT::DynamicMemorySizeInBytes() const {
  return SpinerTricks::DynamicMemorySizeInBytes(this);
}

inline std::size_t SpinerEOSDependsRhoT::DumpDynamicMemory(char *dst) {
  return SpinerTricks::DumpDynamicMemory(dst, this);
}

inline std::size_t
SpinerEOSDependsRhoT::SetDynamicMemory(char *src, const SharedMemSettings &stngs) {
  return SpinerTricks::SetDynamicMemory((stngs.data == nullptr) ? src : stngs.data, this);
}

inline herr_t SpinerEOSDependsRhoT::loadDataboxes_(const std::string &matid_str,
                                                   hid_t file, hid_t lTGroup,
                                                   hid_t coldGroup) {
  herr_t status = H5_SUCCESS;

  // offsets
  status +=
      H5LTget_attribute_double(file, matid_str.c_str(), SP5::Offsets::rho, &lRhoOffset_);
  status +=
      H5LTget_attribute_double(file, matid_str.c_str(), SP5::Offsets::T, &lTOffset_);
  lRhoOffset_ = std::abs(lRhoOffset_);
  lTOffset_ = std::abs(lTOffset_);
  // normal density
  status += H5LTget_attribute_double(file, matid_str.c_str(),
                                     SP5::Material::normalDensity, &rhoNormal_);
  rhoNormal_ = std::abs(rhoNormal_);
  // Mean atomic mass and mean atomic number
  status += H5LTget_attribute_double(file, matid_str.c_str(),
                                     SP5::Material::meanAtomicMass, &(AZbar_.Abar));
  status += H5LTget_attribute_double(file, matid_str.c_str(),
                                     SP5::Material::meanAtomicNumber, &(AZbar_.Zbar));

  // tables
  status += P_.loadHDF(lTGroup, SP5::Fields::P);
  status += sie_.loadHDF(lTGroup, SP5::Fields::sie);
  status += bMod_.loadHDF(lTGroup, SP5::Fields::bMod);
  status += dPdRho_.loadHDF(lTGroup, SP5::Fields::dPdRho);
  status += dPdE_.loadHDF(lTGroup, SP5::Fields::dPdE);
  status += dTdRho_.loadHDF(lTGroup, SP5::Fields::dTdRho);
  status += dTdE_.loadHDF(lTGroup, SP5::Fields::dTdE);
  status += dEdRho_.loadHDF(lTGroup, SP5::Fields::dEdRho);
  status += dEdT_.loadHDF(lTGroup, SP5::Fields::dEdT);

  // cold curves
  status += PCold_.loadHDF(coldGroup, SP5::Fields::P);
  status += sieCold_.loadHDF(coldGroup, SP5::Fields::sie);
  status += bModCold_.loadHDF(coldGroup, SP5::Fields::bMod);
  status += dPdRhoCold_.loadHDF(coldGroup, SP5::Fields::dPdRho);

  numRho_ = bMod_.dim(2);
  numT_ = bMod_.dim(1);

  // set bounds
  lRhoMin_ = P_.range(1).min();
  lRhoMax_ = P_.range(1).max();
  rhoMax_ = fromLog_(lRhoMax_, lRhoOffset_);
  lTMin_ = P_.range(0).min();
  lTMax_ = P_.range(0).max();
  TMax_ = fromLog_(lTMax_, lTOffset_);

  Real rhoMin = fromLog_(lRhoMin_, lRhoOffset_);
  Real rhoMinSearch = std::max(
      rhoMin, std::max(std::abs(robust::EPS()) * 10, std::abs(robust::EPS() * rhoMin)));
  lRhoMinSearch_ = toLog_(rhoMinSearch, lRhoOffset_);

  // bulk modulus can be wrong in the tables. Use FLAG's approach to
  // fix the table.
  fixBulkModulus_();

  // find critical temperature Tcrit(rho)
  // where sie(rho,Tcrit(rho)) = sieCold(rho)
  setlTColdCrit_();

  // fill in minimum pressure as a function of temperature
  rho_at_pmin_.resize(numT_);
  rho_at_pmin_.setRange(0, P_.range(0));
  for (int i = 0; i < numT_; i++) {
    PMin_ = std::numeric_limits<Real>::max();
    int jmax = -1;
    for (int j = 0; j < numRho_; j++) {
      if (P_(j, i) < PMin_) {
        PMin_ = P_(j, i);
        jmax = j;
      }
    }
    if (jmax < 0) printf("Failed to find minimum pressure.\n");
    rho_at_pmin_(i) = rho_(P_.range(1).x(jmax));
  }

  // fill in Gruneisen parameter and bulk modulus on cold curves
  // unfortunately, EOSPAC's output for these parameters appears
  // unreliable we fix it by using constant extrapolation of our
  // values from lTColdCrit
  // dT/dRho and dT/dSie are obviously not physically defined on the
  // cold curve, or they're physically zero. But we have to return
  // something sensible. We extrapolate use constant extrapolation
  // from lTColdCrit.
  // TODO(JMM): the right thing to do here depends on the PTE
  // solver. Maybe think about PTE-solver dependent settings.
  dPdECold_.copyMetadata(bModCold_);
  dTdRhoCold_.copyMetadata(bModCold_);
  dTdECold_.copyMetadata(bModCold_);
  dEdTCold_.copyMetadata(bModCold_);
  for (int j = 0; j < numRho_; j++) {
    Real lRho = bModCold_.range(0).x(j);
    Real lT = lTColdCrit_(j);
    Real rho = rho_(lRho);
    bModCold_(j) = bMod_.interpToReal(lRho, lT);
    dPdECold_(j) = dPdE_.interpToReal(lRho, lT);
    dTdRhoCold_(j) = dTdRho_.interpToReal(lRho, lT);
    dTdECold_(j) = dTdE_.interpToReal(lRho, lT);
    dEdTCold_(j) = dEdT_.interpToReal(lRho, lT);
  }

  // major vs. minor axes change, so this must be done by hand
  PMax_.resize(numRho_);
  PMax_.setRange(0, bMod_.range(1));
  dEdTMax_.copyMetadata(PMax_);
  gm1Max_.copyMetadata(PMax_);
  sielTMax_.copyMetadata(PMax_);
  for (int j = 0; j < numRho_; j++) {
    Real lRho = PMax_.range(0).x(j);
    Real rho = rho_(lRho);
    PMax_(j) = P_(j, numT_ - 1);
    dEdTMax_(j) = dEdT_(j, numT_ - 1);
    gm1Max_(j) = robust::ratio(dPdE_(j, numT_ - 1), rho); // max gruneisen
    sielTMax_(j) = sie_(j, numT_ - 1);
  }

  // reference state
  Real lRhoNormal = lRho_(rhoNormal_);
  // if rho normal not on the table, set it to the middle
  if (!(lRhoMin_ < lRhoNormal && lRhoNormal < lRhoMax_)) {
    lRhoNormal = 0.5 * (lRhoMin_ + lRhoMax_);
    rhoNormal_ = rho_(lRhoNormal);
  }
  // Same for temperature. Use room temperature if it's available
  TNormal_ = ROOM_TEMPERATURE;
  Real lTNormal = lT_(TNormal_);
  if (!(lTMin_ < lTNormal && lTNormal < lTMax_)) {
    lTNormal = 0.5 * (lTMin_ + lTMax_);
    TNormal_ = T_(lTNormal);
  }
  sieNormal_ = sie_.interpToReal(lRhoNormal, lTNormal);
  PNormal_ = P_.interpToReal(lRhoNormal, lTNormal);
  CvNormal_ = dEdT_.interpToReal(lRhoNormal, lTNormal);
  bModNormal_ = bMod_.interpToReal(lRhoNormal, lTNormal);
  dPdENormal_ = dPdE_.interpToReal(lRhoNormal, lTNormal);
  Real dPdR = dPdRho_.interpToReal(lRhoNormal, lTNormal);
  dVdTNormal_ = dPdENormal_ * CvNormal_ / (rhoNormal_ * rhoNormal_ * dPdR);

  return status;
}

inline void SpinerEOSDependsRhoT::fixBulkModulus_() {
  // assumes all databoxes are the same size
  // TODO: do we need to smooth this data with a median filter
  // or something like that?
  for (int j = 0; j < numRho_; j++) {
    Real lRho = bMod_.range(1).x(j);
    Real rho = rho_(lRho);
    for (int i = 0; i < numT_; i++) {
      Real lT = bMod_.range(0).x(i);
      Real press = P_.interpToReal(lRho, lT);
      Real DPDR_E = dPdRho_.interpToReal(lRho, lT);
      Real DPDE_R = dPdE_.interpToReal(lRho, lT);
      Real DEDR_T = dEdRho_.interpToReal(lRho, lT);
      Real DPDR_T = DPDR_E + DPDE_R * DEDR_T;
      Real bMod;
      if (DPDE_R > 0.0 && rho > 0.0) {
        bMod = rho * DPDR_E + DPDE_R * (press / rho);
      } else if (rho > 0.0) {
        bMod = std::max(rho * DPDR_T, 0.0);
      } else {
        bMod = 0.0;
      }
      bMod_(j, i) = std::max(bMod, std::abs(robust::EPS()));
    }
  }
}

inline void SpinerEOSDependsRhoT::setlTColdCrit_() {
  lTColdCrit_.copyMetadata(bModCold_);
  for (int j = 0; j < numRho_; j++) {
    Real lRho = bModCold_.range(0).x(j);
    Real sieCold = sieCold_.interpToReal(lRho);
    // First find all "zero-crossings"
    // where sie - sieCold = 0
    // there could be:
    // 1. None
    // 2. Exactly 1
    // 3. Many
    // depending on how well-behaved the EOS is
    // None happens if sieCold is off the grid or at the min of the grid.
    // exactly 1 happens if sieCold is on the grid
    // many happens if the EOS has non-monotone sie and/or crossing isotherms
    // TODO(JMM): I do this in linear time because it's done once at startup.
    // if this turns out to be too costly, use modified binary search.
    std::vector<int> crossings;
    std::vector<int> directions;
    int last_pos_crossing = -1;
    for (int i = 0; i < numT_ - 1; i++) {
      Real sign_curr = sie_(j, i) - sieCold;
      Real sign_next = sie_(j, i + 1) - sieCold;
      if (sign_curr < 0 && sign_next > 0) {
        crossings.push_back(i);
        directions.push_back(1);
        last_pos_crossing = crossings.size() - 1;
      } else if (sign_curr > 0 && sign_next < 0) {
        crossings.push_back(i);
        directions.push_back(-1);
      }
    }
    // react accordingly
    // TODO(JMM): If we have multiple zero crossings,
    // there's no good uniquely defined solution.
    // We choose the highest-temperature crossing,
    // which results in the least code. But this may be wrong.
    // JMM: <= 0 required here so that a cold curve on the min T
    // isotherm doesn't trigger a root find.
    if (last_pos_crossing <= 0) { // off the grid
      lTColdCrit_(j) = lTMin_;
    } else { // at least one pos crossing. Use last one.
      const callable_interp::r_interp sieFunc(sie_, lRho);
      Real lT;
      int ilast = crossings[last_pos_crossing];
      // expand bounds by +/- 1.e-14 to help with round-off
      Real lTlower = bMod_.range(0).x(ilast) - 1.0e-14;
      Real lTupper = bMod_.range(0).x(ilast + 1) + 1.0e-14;
      Real lTGuess = 0.5 * (lTlower + lTupper);
      auto status = ROOT_FINDER(sieFunc, sieCold, lTGuess, lTlower, lTupper, ROOT_THRESH,
                                ROOT_THRESH, lT);
      if (status != RootFinding1D::Status::SUCCESS) {
        lT = lTGuess;
      }
      lTColdCrit_(j) = lT;
      // lTColdCrit_(j) = lTMin_;
    }
    lTColdCrit_(j) = std::max(lTMin_, lTColdCrit_(j));
  }
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoT::TemperatureFromDensityInternalEnergy(
    const Real rho, const Real sie, Indexer_t &&lambda) const {
  TableStatus whereAmI;
  const Real lRho = lRho_(rho);
  const Real lT = lTFromlRhoSie_(lRho, sie, whereAmI, lambda);
  return T_(lT);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoT::InternalEnergyFromDensityTemperature(
    const Real rho, const Real temperature, Indexer_t &&lambda) const {
  Real lRho, lT;
  getLogsRhoT_(rho, temperature, lRho, lT, lambda);
  TableStatus whereAmI = getLocDependsRhoT_(lRho, lT);
  return sieFromlRhoTlT_(lRho, temperature, lT, whereAmI);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoT::PressureFromDensityTemperature(
    const Real rho, const Real temperature, Indexer_t &&lambda) const {
  Real lRho, lT;
  getLogsRhoT_(rho, temperature, lRho, lT, lambda);
  TableStatus whereAmI = getLocDependsRhoT_(lRho, lT);
  return PFromRholRhoTlT_(rho, lRho, temperature, lT, whereAmI);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoT::PressureFromDensityInternalEnergy(
    const Real rho, const Real sie, Indexer_t &&lambda) const {
  TableStatus whereAmI;
  Real lRho = lRho_(rho);
  Real lT = lTFromlRhoSie_(lRho, sie, whereAmI, lambda);
  Real P;
  if (whereAmI == TableStatus::OffBottom) { // cold curve
    P = PCold_.interpToReal(lRho);
  } else if (whereAmI == TableStatus::OffTop) { // ideal gas
    const Real gm1 = gm1Max_.interpToReal(lRho);
    P = gm1 * rho * sie;
  } else { // on table
    P = P_.interpToReal(lRho, lT);
  }
  return P;
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoT::MinInternalEnergyFromDensity(
    const Real rho, Indexer_t &&lambda) const {
  const Real lRho = lRho_(rho);
  return sieCold_.interpToReal(lRho);
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoT::EntropyFromDensityTemperature(
    const Real rho, const Real temperature, Indexer_t &&lambda) const {
  EntropyIsNotEnabled("SpinerEOSDependsRhoT");
  return 1.0;
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoT::EntropyFromDensityInternalEnergy(
    const Real rho, const Real sie, Indexer_t &&lambda) const {
  EntropyIsNotEnabled("SpinerEOSDependsRhoT");
  return 1.0;
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoT::SpecificHeatFromDensityTemperature(
    const Real rho, const Real temperature, Indexer_t &&lambda) const {
  Real lRho, lT;
  getLogsRhoT_(rho, temperature, lRho, lT, lambda);
  TableStatus whereAmI = getLocDependsRhoT_(lRho, lT);
  return CvFromlRholT_(lRho, lT, whereAmI);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoT::SpecificHeatFromDensityInternalEnergy(
    const Real rho, const Real sie, Indexer_t &&lambda) const {
  TableStatus whereAmI;
  Real Cv;
  const Real lRho = lRho_(rho);
  const Real lT = lTFromlRhoSie_(lRho, sie, whereAmI, lambda);
  if (whereAmI == TableStatus::OffBottom) { // cold curve
    // on cold curve. Currently, we assume constant extrapolation.
    // TODO(JMM): Do something more sophisticated
    Cv = dEdTCold_.interpToReal(lRho);
  } else if (whereAmI == TableStatus::OffTop) { // ideal gas
    Cv = dEdTMax_.interpToReal(lRho);           // Cv assumed constant in T
  } else {                                      // on table
    Cv = dEdT_.interpToReal(lRho, lT);
  }
  return Cv > robust::EPS() ? Cv : robust::EPS();
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoT::BulkModulusFromDensityTemperature(
    const Real rho, const Real temperature, Indexer_t &&lambda) const {
  Real lRho, lT;
  getLogsRhoT_(rho, temperature, lRho, lT, lambda);
  TableStatus whereAmI = getLocDependsRhoT_(lRho, lT);
  return bModFromRholRhoTlT_(rho, lRho, temperature, lT, whereAmI);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoT::GruneisenParamFromDensityTemperature(
    const Real rho, const Real temp, Indexer_t &&lambda) const {
  Real lRho, lT, gm1;
  getLogsRhoT_(rho, temp, lRho, lT, lambda);

  TableStatus whereAmI = getLocDependsRhoT_(lRho, lT);
  if (whereAmI == TableStatus::OffBottom) {
    // use cold curves
    Real dpde = dPdECold_.interpToReal(lRho);
    gm1 = robust::ratio(std::abs(dpde), std::abs(rho));
  } else if (whereAmI == TableStatus::OffTop) {
    gm1 = gm1Max_.interpToReal(lRho);
  } else { // on table
    const Real dpde = dPdE_.interpToReal(lRho, lT);
    gm1 = robust::ratio(std::abs(dpde), std::abs(rho));
  }
  return gm1;
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoT::BulkModulusFromDensityInternalEnergy(
    const Real rho, const Real sie, Indexer_t &&lambda) const {
  TableStatus whereAmI;
  Real bMod;
  const Real lRho = lRho_(rho);
  const Real lT = lTFromlRhoSie_(lRho, sie, whereAmI, lambda);
  if (whereAmI == TableStatus::OffBottom) {
    bMod = bModCold_.interpToReal(lRho);
  } else if (whereAmI == TableStatus::OffTop) {
    const Real gm1 = gm1Max_.interpToReal(lRho);
    bMod = (gm1 + 1) * gm1 * rho * sie;
  } else { // on table
    bMod = bMod_.interpToReal(rho, sie);
  }
  return bMod;
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real
SpinerEOSDependsRhoT::GruneisenParamFromDensityInternalEnergy(const Real rho,
                                                              const Real sie,
                                                              Indexer_t &&lambda) const {
  TableStatus whereAmI;
  Real gm1;
  const Real lRho = lRho_(rho);
  const Real lT = lTFromlRhoSie_(lRho, sie, whereAmI, lambda);
  if (whereAmI == TableStatus::OffBottom) {
    Real dpde = dPdECold_.interpToReal(lRho);
    gm1 = robust::ratio(std::abs(dpde), std::abs(rho));
  } else if (whereAmI == TableStatus::OffTop) {
    gm1 = gm1Max_.interpToReal(lRho);
  } else {
    const Real dpde = dPdE_.interpToReal(lRho, lT);
    gm1 = robust::ratio(std::abs(dpde), std::abs(rho));
  }
  return gm1;
}

// TODO(JMM): This would be faster with hand-tuned code
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void SpinerEOSDependsRhoT::DensityEnergyFromPressureTemperature(
    const Real press, const Real temp, Indexer_t &&lambda, Real &rho, Real &sie) const {
  TableStatus whereAmI;
  Real lT = lT_(temp);
  Real lRho = lRhoFromPlT_(press, lT, whereAmI, lambda);
  rho = rho_(lRho);
  sie = InternalEnergyFromDensityTemperature(rho, temp, lambda);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void
SpinerEOSDependsRhoT::FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv,
                              Real &bmod, const unsigned long output,
                              Indexer_t &&lambda) const {
  Real lRho, lT;
  TableStatus whereAmI;
  const unsigned long input = ~output;
  if (output == thermalqs::none) {
    UNDEFINED_ERROR;
  }
  if (output & thermalqs::density) {
    if (input & thermalqs::pressure && input & thermalqs::temperature) {
      lT = lT_(temp);
      lRho = lRhoFromPlT_(press, lT, whereAmI, lambda);
      rho = rho_(lRho);
    } else {
      UNDEFINED_ERROR;
    }
  } else {
    lRho = lRho_(rho);
  }
  if (output & thermalqs::temperature) {
    if (input & thermalqs::density && input & thermalqs::specific_internal_energy) {
      lT = lTFromlRhoSie_(lRho, energy, whereAmI, lambda);
    } else if (input & thermalqs::density && input & thermalqs::pressure) {
      lT = lTFromlRhoP_(lRho, press, whereAmI, lambda);
    } else {
      UNDEFINED_ERROR;
    }
    temp = T_(lT);
  } else {
    lT = lT_(temp);
  }
  whereAmI = getLocDependsRhoT_(lRho, lT);
  if (output & thermalqs::specific_internal_energy) {
    energy = sieFromlRhoTlT_(lRho, temp, lT, whereAmI);
  }
  if (output & thermalqs::pressure) {
    press = PFromRholRhoTlT_(rho, lRho, temp, lT, whereAmI);
  }
  if (output & thermalqs::specific_heat) {
    cv = CvFromlRholT_(lRho, lT, whereAmI);
  }
  if (output & thermalqs::bulk_modulus) {
    bmod = bModFromRholRhoTlT_(rho, lRho, temp, lT, whereAmI);
  }
  if (!variadic_utils::is_nullptr(lambda)) {
    IndexerUtils::Get<IndexableTypes::LogDensity>(lambda, Lambda::lRho) = lRho;
    IndexerUtils::Get<IndexableTypes::LogTemperature>(lambda, Lambda::lT) = lT;
  }
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void SpinerEOSDependsRhoT::ValuesAtReferenceState(
    Real &rho, Real &temp, Real &sie, Real &press, Real &cv, Real &bmod, Real &dpde,
    Real &dvdt, Indexer_t &&lambda) const {
  rho = rhoNormal_;
  temp = TNormal_;
  sie = sieNormal_;
  press = PNormal_;
  cv = CvNormal_;
  bmod = bModNormal_;
  dpde = dPdENormal_;
  dvdt = dVdTNormal_;
}

PORTABLE_INLINE_FUNCTION
Real SpinerEOSDependsRhoT::RhoPmin(const Real temp) const {
  const Real lT = lT_(temp);
  if (lT <= lTMin_) return rho_at_pmin_(0);
  if (lT >= lTMax_) return 0.0;
  return rho_at_pmin_.interpToReal(lT);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void
SpinerEOSDependsRhoT::getLogsRhoT_(const Real rho, const Real temperature, Real &lRho,
                                   Real &lT, Indexer_t &&lambda) const {
  lRho = lRho_(rho);
  lT = lT_(temperature);
  if (!variadic_utils::is_nullptr(lambda)) {
    IndexerUtils::Get<IndexableTypes::LogDensity>(lambda, Lambda::lRho) = lRho;
    IndexerUtils::Get<IndexableTypes::LogTemperature>(lambda, Lambda::lT) = lT;
  }
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoT::lRhoFromPlT_(
    const Real P, const Real lT, TableStatus &whereAmI, Indexer_t &&lambda) const {
  RootFinding1D::Status status = RootFinding1D::Status::SUCCESS;

  Real lRho;
  Real lRhoGuess = reproducible_ ? lRhoMax_ : 0.5 * (lRhoMin_ + lRhoMax_);
  // Real lRhoGuess = lRhoMin_ + 0.9*(lRhoMax_ - lRhoMin_);
  const RootFinding1D::RootCounts *pcounts =
      (memoryStatus_ == DataStatus::OnDevice) ? nullptr : &counts;

  if (!variadic_utils::is_nullptr(lambda)) {
    Real lRho_cache = IndexerUtils::Get<IndexableTypes::LogDensity>(lambda, Lambda::lRho);
    if ((lRhoMin_ <= lRho_cache) && (lRho_cache <= lRhoMax_)) {
      lRhoGuess = lRho_cache;
    }
  }

  if (lT <= lTMin_) { // cold curve
    whereAmI = TableStatus::OffBottom;
    const callable_interp::interp PFunc(PCold_);
    status =
        ROOT_FINDER(PFunc, P, lRhoGuess,
                    // lRhoMin_, lRhoMax_,
                    lRhoMinSearch_, lRhoMax_, ROOT_THRESH, ROOT_THRESH, lRho, pcounts);
  } else if (lT >= lTMax_) { // ideal gas
    whereAmI = TableStatus::OffTop;
    const callable_interp::prod_interp_1d PFunc(gm1Max_, dEdTMax_, lT);
    status =
        ROOT_FINDER(PFunc, P, lRhoGuess,
                    // lRhoMin_, lRhoMax_,
                    lRhoMinSearch_, lRhoMax_, ROOT_THRESH, ROOT_THRESH, lRho, pcounts);
  } else { // on table
    whereAmI = TableStatus::OnTable;
    const callable_interp::l_interp PFunc(P_, lT);
    status =
        ROOT_FINDER(PFunc, P, lRhoGuess,
                    // lRhoMin_, lRhoMax_,
                    lRhoMinSearch_, lRhoMax_, ROOT_THRESH, ROOT_THRESH, lRho, pcounts);
  }
  if (status != RootFinding1D::Status::SUCCESS) {
#if SPINER_EOS_VERBOSE
    std::stringstream errorMessage;
    errorMessage << "inverting P table for logRho failed\n"
                 << "matid     = " << matid_ << "\n"
                 << "lT        = " << lT << "\n"
                 << "P         = " << P << "\n"
                 << "lRhoGuess = " << lRhoGuess << std::endl;
    EOS_ERROR(errorMessage.str().c_str());
#endif // SPINER_EOS_VERBOSE
    lRho = reproducible_ ? lRhoMax_ : lRhoGuess;
  }
  if (!variadic_utils::is_nullptr(lambda)) {
    IndexerUtils::Get<IndexableTypes::LogDensity>(lambda, Lambda::lRho) = lRho;
    IndexerUtils::Get<IndexableTypes::LogTemperature>(lambda, Lambda::lT) = lT;
    if constexpr (variadic_utils::is_indexable_v<Indexer_t, IndexableTypes::RootStatus>) {
      lambda[IndexableTypes::RootStatus()] = static_cast<Real>(status);
    }
    if constexpr (variadic_utils::is_indexable_v<Indexer_t,
                                                 IndexableTypes::TableStatus>) {
      lambda[IndexableTypes::TableStatus()] = static_cast<Real>(whereAmI);
    }
  }
  return lRho;
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoT::lTFromlRhoSie_(
    const Real lRho, const Real sie, TableStatus &whereAmI, Indexer_t &&lambda) const {

  const RootFinding1D::RootCounts *pcounts =
      (memoryStatus_ == DataStatus::OnDevice) ? nullptr : &counts;
  RootFinding1D::Status status = RootFinding1D::Status::SUCCESS;
  Real lT;

  whereAmI = getLocDependsRhoSie_(lRho, sie);
  if (whereAmI == TableStatus::OffBottom) {
    // On the cold curve. No unique temperature.
    // return the minimum temperature in the table.
    lT = lTMin_;
    if (pcounts != nullptr) {
      pcounts->increment(0);
    }
  } else if (whereAmI == TableStatus::OffTop) { // Assume ideal gas
    const Real Cv = dEdTMax_.interpToReal(lRho);
    const Real e0 = sielTMax_.interpToReal(lRho);
    const Real T = TMax_ + robust::ratio(sie - e0, Cv);
    lT = lT_(T);
    if (pcounts != nullptr) {
      pcounts->increment(0);
    }
  } else {
    Real lTGuess = reproducible_ ? lTMin_ : 0.5 * (lTMin_ + lTMax_);
    if (!variadic_utils::is_nullptr(lambda)) {
      Real lT_cache =
          IndexerUtils::Get<IndexableTypes::LogTemperature>(lambda, Lambda::lT);
      if ((lTMin_ <= lT_cache) && (lT_cache <= lTMax_)) {
        lTGuess = lT_cache;
      }
    }
    const callable_interp::r_interp sieFunc(sie_, lRho);
    status = ROOT_FINDER(sieFunc, sie, lTGuess, lTMin_, lTMax_, ROOT_THRESH, ROOT_THRESH,
                         lT, pcounts);

    if (status != RootFinding1D::Status::SUCCESS) {
#if SPINER_EOS_VERBOSE
      std::stringstream errorMessage;
      errorMessage << std::scientific << std::setprecision(14)
                   << "inverting sie table for logT failed\n"
                   << "matid   = " << matid_ << "\n"
                   << "lRho    = " << lRho << "\n"
                   << "sie     = " << sie << "\n"
                   << "lTGuess = " << lTGuess << "\n"
                   << "sielTMax = " << sielTMax_.interpToReal(lRho) << "\n"
                   << "sieCold = " << sieCold_.interpToReal(lRho) << std::endl;
      EOS_ERROR(errorMessage.str().c_str());
#endif // SPINER_EOS_VERBOSE
      lT = reproducible_ ? lTMin_ : lTGuess;
    }
  }
  if (!variadic_utils::is_nullptr(lambda)) {
    IndexerUtils::Get<IndexableTypes::LogDensity>(lambda, Lambda::lRho) = lRho;
    IndexerUtils::Get<IndexableTypes::LogTemperature>(lambda, Lambda::lT) = lT;
    if constexpr (variadic_utils::is_indexable_v<Indexer_t, IndexableTypes::RootStatus>) {
      lambda[IndexableTypes::RootStatus()] = static_cast<Real>(status);
    }
    if constexpr (variadic_utils::is_indexable_v<Indexer_t,
                                                 IndexableTypes::TableStatus>) {
      lambda[IndexableTypes::TableStatus()] = static_cast<Real>(whereAmI);
    }
  }
  return lT;
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoT::lTFromlRhoP_(
    const Real lRho, const Real press, TableStatus &whereAmI, Indexer_t &&lambda) const {
  const RootFinding1D::RootCounts *pcounts =
      (memoryStatus_ == DataStatus::OnDevice) ? nullptr : &counts;
  RootFinding1D::Status status = RootFinding1D::Status::SUCCESS;
  Real lT, lTGuess;

  // Assumes P is monotone in T
  const Real PCold = PCold_.interpToReal(lRho);
  const Real PMax = PMax_.interpToReal(lRho);
  if (press <= PCold) {
    whereAmI = TableStatus::OffBottom;
    lT = lTMin_;
    if (pcounts != nullptr) {
      pcounts->increment(0);
    }
  } else if (press >= PMax) {
    whereAmI = TableStatus::OffTop;
    lT = lTMax_;
    if (pcounts != nullptr) {
      pcounts->increment(0);
    }
  } else {
    whereAmI = TableStatus::OnTable;
    lTGuess = 0.5 * (lTMin_ + lTMax_);
    if (!variadic_utils::is_nullptr(lambda)) {
      Real lT_cache = IndexerUtils::Get<IndexableTypes::LogTemperature>(lambda, lT);
      if ((lTMin_ <= lT_cache) && (lT_cache <= lTMax_)) {
        lTGuess = lT_cache;
      }
    }
    const callable_interp::r_interp PFunc(P_, lRho);
    status = ROOT_FINDER(PFunc, press, lTGuess, lTMin_, lTMax_, ROOT_THRESH, ROOT_THRESH,
                         lT, pcounts);
    if (status != RootFinding1D::Status::SUCCESS) {
#if SPINER_EOS_VERBOSE
      std::stringstream errorMessage;
      errorMessage << "inverting P table for logT failed\n"
                   << "matid   = " << matid_ << "\n"
                   << "lRho    = " << lRho << "\n"
                   << "P       = " << press << "\n"
                   << "lTGuess = " << lTGuess << std::endl;
      EOS_ERROR(errorMessage.str().c_str());
#endif // SPINER_EOS_VERBOSE
      lT = reproducible_ ? lTMin_ : lTGuess;
    }
  }
  if (!variadic_utils::is_nullptr(lambda)) {
    IndexerUtils::Get<IndexableTypes::LogDensity>(lambda, Lambda::lRho) = lRho;
    IndexerUtils::Get<IndexableTypes::LogTemperature>(lambda, Lambda::lT) = lT;
    if constexpr (variadic_utils::is_indexable_v<Indexer_t, IndexableTypes::RootStatus>) {
      lambda[IndexableTypes::RootStatus()] = static_cast<Real>(status);
    }
    if constexpr (variadic_utils::is_indexable_v<Indexer_t,
                                                 IndexableTypes::TableStatus>) {
      lambda[IndexableTypes::TableStatus()] = static_cast<Real>(whereAmI);
    }
  }
  return lT;
}

PORTABLE_INLINE_FUNCTION
Real SpinerEOSDependsRhoT::sieFromlRhoTlT_(const Real lRho, const Real T, const Real lT,
                                           const TableStatus &whereAmI) const {
  Real sie;
  if (whereAmI == TableStatus::OffBottom) { // cold curve
    sie = sieCold_.interpToReal(lRho);
  } else if (whereAmI == TableStatus::OffTop) { // ideal gas
    const Real Cv = dEdTMax_.interpToReal(lRho);
    const Real e0 = sielTMax_.interpToReal(lRho);
    sie = e0 + Cv * (T - TMax_);
  } else { // on table
    sie = sie_.interpToReal(lRho, lT);
  }
  return sie;
}

PORTABLE_INLINE_FUNCTION
Real SpinerEOSDependsRhoT::PFromRholRhoTlT_(const Real rho, const Real lRho, const Real T,
                                            const Real lT,
                                            const TableStatus &whereAmI) const {
  Real P;
  if (whereAmI == TableStatus::OffBottom) { // cold curve
    P = PCold_.interpToReal(lRho);
  } else if (whereAmI == TableStatus::OffTop) { // ideal gas
    const Real Cv = dEdTMax_.interpToReal(lRho);
    const Real gm1 = gm1Max_.interpToReal(lRho);
    const Real e0 = sielTMax_.interpToReal(lRho);
    const Real e = e0 + Cv * (T - TMax_);
    P = gm1 * rho * e;
  } else { // if ( whereAmI == TableStatus::OnTable) {
    P = P_.interpToReal(lRho, lT);
  }
  return P;
}

PORTABLE_INLINE_FUNCTION
Real SpinerEOSDependsRhoT::CvFromlRholT_(const Real lRho, const Real lT,
                                         const TableStatus &whereAmI) const {
  Real Cv;
  if (whereAmI == TableStatus::OffBottom) { // cold curve
    // on cold curve. Currently, we assume constant extrapolation.
    // TODO(JMM): Do something more sophisticated
    Cv = dEdTCold_.interpToReal(lRho);
  } else if (whereAmI == TableStatus::OffTop) { // ideal gas
    Cv = dEdTMax_.interpToReal(lRho);           // Cv assumed constant in T
  } else {                                      // on table
    Cv = dEdT_.interpToReal(lRho, lT);
  }
  return Cv > robust::EPS() ? Cv : robust::EPS();
}

PORTABLE_INLINE_FUNCTION
Real SpinerEOSDependsRhoT::bModFromRholRhoTlT_(const Real rho, const Real lRho,
                                               const Real T, const Real lT,
                                               const TableStatus &whereAmI) const {
  Real bMod;
  if (whereAmI == TableStatus::OffBottom) {
    bMod = bModCold_.interpToReal(lRho);
  } else if (whereAmI == TableStatus::OffTop) {
    const Real Cv = dEdTMax_.interpToReal(lRho);
    const Real gm1 = gm1Max_.interpToReal(lRho);
    const Real e0 = sielTMax_.interpToReal(lRho);
    const Real e = e0 + Cv * (T - TMax_);
    bMod = (gm1 + 1) * gm1 * rho * e;
  } else { // on table
    bMod = bMod_.interpToReal(lRho, lT);
  }
  return bMod > robust::EPS() ? bMod : robust::EPS();
}

PORTABLE_INLINE_FUNCTION
TableStatus SpinerEOSDependsRhoT::getLocDependsRhoSie_(const Real lRho,
                                                       const Real sie) const {
  TableStatus whereAmI;
  const Real sielTMax = sielTMax_.interpToReal(lRho);
  const Real sieCold = sieCold_.interpToReal(lRho);
  // sie can be negative, so must make sign right
  if (sie >= sielTMax - SOFT_THRESH * std::abs(sielTMax)) {
    whereAmI = TableStatus::OffTop;
  } else if (sie <= sieCold + SOFT_THRESH * std::abs(sieCold)) {
    whereAmI = TableStatus::OffBottom;
  } else {
    whereAmI = TableStatus::OnTable;
  }
  return whereAmI;
}

PORTABLE_INLINE_FUNCTION TableStatus
SpinerEOSDependsRhoT::getLocDependsRhoT_(const Real lRho, const Real lT) const {
  TableStatus whereAmI;
  if (lT <= (1 + SOFT_THRESH) * lTMin_)
    whereAmI = TableStatus::OffBottom;
  else if (lT >= (1 - SOFT_THRESH) * lTMax_)
    whereAmI = TableStatus::OffTop;
  else
    whereAmI = TableStatus::OnTable;
  return whereAmI;
}

inline SpinerEOSDependsRhoSie::SpinerEOSDependsRhoSie(const std::string &filename,
                                                      int matid, TableSplit split,
                                                      bool reproducibility_mode)
    : matid_(matid), split_(split), reproducible_(reproducibility_mode),
      memoryStatus_(DataStatus::OnHost) {

  std::string matid_str = std::to_string(matid);
  hid_t file, matGroup, lTGroup, lEGroup;
  herr_t status = H5_SUCCESS;

  file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  matGroup = H5Gopen(file, matid_str.c_str(), H5P_DEFAULT);

  std::string lTGroupName = SP5::Depends::logRhoLogT;
  std::string lEGroupName = SP5::Depends::logRhoLogSie;
  if (split == TableSplit::ElectronOnly) {
    lTGroupName += (std::string("/") + SP5::SubTable::electronOnly);
    lEGroupName += (std::string("/") + SP5::SubTable::electronOnly);
  } else if (split == TableSplit::IonCold) {
    lTGroupName += (std::string("/") + SP5::SubTable::ionCold);
    lEGroupName += (std::string("/") + SP5::SubTable::ionCold);
  }

  lTGroup = H5Gopen(matGroup, lTGroupName.c_str(), H5P_DEFAULT);
  lEGroup = H5Gopen(matGroup, lEGroupName.c_str(), H5P_DEFAULT);

  status += loadDataboxes_(matid_str, file, lTGroup, lEGroup);

  status += H5Gclose(lTGroup);
  status += H5Gclose(lEGroup);
  status += H5Gclose(matGroup);
  status += H5Fclose(file);

  if (status != H5_SUCCESS) {
    EOS_ERROR("SpinerDependsTHoSIE: HDF5 error\n");
  }
}

inline SpinerEOSDependsRhoSie::SpinerEOSDependsRhoSie(const std::string &filename,
                                                      const std::string &materialName,
                                                      TableSplit split,
                                                      bool reproducibility_mode)
    : split_(split), reproducible_(reproducibility_mode),
      memoryStatus_(DataStatus::OnHost) {

  std::string matid_str;
  hid_t file, matGroup, lTGroup, lEGroup;
  herr_t status = H5_SUCCESS;

  file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  matGroup = H5Gopen(file, materialName.c_str(), H5P_DEFAULT);

  std::string lTGroupName = SP5::Depends::logRhoLogT;
  std::string lEGroupName = SP5::Depends::logRhoLogSie;
  if (split == TableSplit::ElectronOnly) {
    lTGroupName += (std::string("/") + SP5::SubTable::electronOnly);
    lEGroupName += (std::string("/") + SP5::SubTable::electronOnly);
  } else if (split == TableSplit::IonCold) {
    lTGroupName += (std::string("/") + SP5::SubTable::ionCold);
    lEGroupName += (std::string("/") + SP5::SubTable::ionCold);
  }

  lTGroup = H5Gopen(matGroup, lTGroupName.c_str(), H5P_DEFAULT);
  lEGroup = H5Gopen(matGroup, lEGroupName.c_str(), H5P_DEFAULT);

  status +=
      H5LTget_attribute_int(file, materialName.c_str(), SP5::Material::matid, &matid_);
  matid_str = std::to_string(matid_);

  status += loadDataboxes_(matid_str, file, lTGroup, lEGroup);

  status += H5Gclose(lTGroup);
  status += H5Gclose(lEGroup);
  status += H5Gclose(matGroup);
  status += H5Fclose(file);

  if (status != H5_SUCCESS) {
    EOS_ERROR("SpinerDependsRhoSie: HDF5 error\n");
  }
}

herr_t SpinerEOSDependsRhoSie::loadDataboxes_(const std::string &matid_str, hid_t file,
                                              hid_t lTGroup, hid_t lEGroup) {
  herr_t status = H5_SUCCESS;

  // offsets
  status +=
      H5LTget_attribute_double(file, matid_str.c_str(), SP5::Offsets::rho, &lRhoOffset_);
  status +=
      H5LTget_attribute_double(file, matid_str.c_str(), SP5::Offsets::T, &lTOffset_);
  status +=
      H5LTget_attribute_double(file, matid_str.c_str(), SP5::Offsets::sie, &lEOffset_);
  lRhoOffset_ = std::abs(lRhoOffset_);
  lTOffset_ = std::abs(lTOffset_);
  lEOffset_ = std::abs(lEOffset_);

  // normal density
  status += H5LTget_attribute_double(file, matid_str.c_str(),
                                     SP5::Material::normalDensity, &rhoNormal_);
  rhoNormal_ = std::abs(rhoNormal_);
  // Mean atomic mass and mean atomic number
  status += H5LTget_attribute_double(file, matid_str.c_str(),
                                     SP5::Material::meanAtomicMass, &(AZbar_.Abar));
  status += H5LTget_attribute_double(file, matid_str.c_str(),
                                     SP5::Material::meanAtomicNumber, &(AZbar_.Zbar));

  // sometimes independent variables
  status += sie_.loadHDF(lTGroup, SP5::Fields::sie);
  status += T_.loadHDF(lEGroup, SP5::Fields::T);

  // dependent variables
  // depends on rho and T
  status += dependsRhoT_.P.loadHDF(lTGroup, SP5::Fields::P);
  status += dependsRhoT_.bMod.loadHDF(lTGroup, SP5::Fields::bMod);
  status += dependsRhoT_.dPdRho.loadHDF(lTGroup, SP5::Fields::dPdRho);
  status += dependsRhoT_.dPdE.loadHDF(lTGroup, SP5::Fields::dPdE);
  status += dependsRhoT_.dTdRho.loadHDF(lTGroup, SP5::Fields::dTdRho);
  status += dependsRhoT_.dTdE.loadHDF(lTGroup, SP5::Fields::dTdE);
  status += dependsRhoT_.dEdRho.loadHDF(lTGroup, SP5::Fields::dEdRho);
  // depends on rho and e
  status += dependsRhoSie_.P.loadHDF(lEGroup, SP5::Fields::P);
  status += dependsRhoSie_.bMod.loadHDF(lEGroup, SP5::Fields::bMod);
  status += dependsRhoSie_.dPdRho.loadHDF(lEGroup, SP5::Fields::dPdRho);
  status += dependsRhoSie_.dPdE.loadHDF(lEGroup, SP5::Fields::dPdE);
  status += dependsRhoSie_.dTdRho.loadHDF(lEGroup, SP5::Fields::dTdRho);
  status += dependsRhoSie_.dTdE.loadHDF(lEGroup, SP5::Fields::dTdE);
  status += dependsRhoSie_.dEdRho.loadHDF(lEGroup, SP5::Fields::dEdRho);

  // Fix up bulk modulus
  calcBMod_(dependsRhoT_);
  calcBMod_(dependsRhoSie_);

  // Metadata for root finding extrapolation
  numRho_ = sie_.dim(2);
  numT_ = sie_.dim(1);
  lRhoMin_ = sie_.range(1).min();
  lRhoMax_ = sie_.range(1).max();
  rhoMax_ = fromLog_(lRhoMax_, lRhoOffset_);

  // slice to maximum of rho
  PlRhoMax_ = dependsRhoT_.P.slice(numRho_ - 1);
  dPdRhoMax_ = dependsRhoT_.dPdRho.slice(numRho_ - 1);

  // fill in minimum pressure as a function of temperature
  rho_at_pmin_.resize(numT_);
  rho_at_pmin_.setRange(0, sie_.range(0));
  for (int i = 0; i < numT_; i++) {
    PMin_ = std::numeric_limits<Real>::max();
    int jmax = -1;
    for (int j = 0; j < numRho_; j++) {
      if (dependsRhoT_.P(j, i) < PMin_) {
        PMin_ = dependsRhoT_.P(j, i);
        jmax = j;
      }
    }
    if (jmax < 0) printf("Failed to find minimum pressure.\n");
    rho_at_pmin_(i) = fromLog_(dependsRhoT_.P.range(1).x(jmax), lRhoOffset_);
  }

  // reference state
  Real lRhoNormal = toLog_(rhoNormal_, lRhoOffset_);
  // if rho normal not on the table, set it to the middle
  if (!(lRhoMin_ < lRhoNormal && lRhoNormal < lRhoMax_)) {
    lRhoNormal = 0.5 * (lRhoMin_ + lRhoMax_);
    rhoNormal_ = fromLog_(lRhoNormal, lRhoOffset_);
  }
  // Same for temperature. Use room temperature if it's available
  TNormal_ = ROOM_TEMPERATURE;
  Real lTNormal = toLog_(TNormal_, lTOffset_);
  Real lTMin = sie_.range(0).min();
  Real lTMax = sie_.range(0).max();
  if (!(lTMin < lTNormal && lTNormal < lTMax)) {
    lTNormal = 0.5 * (lTMin + lTMax);
    TNormal_ = fromLog_(lTNormal, lTOffset_);
  }
  sieNormal_ = sie_.interpToReal(lRhoNormal, lTNormal);
  PNormal_ = dependsRhoT_.P.interpToReal(lRhoNormal, lTNormal);
  CvNormal_ = 1. / dependsRhoT_.dTdE.interpToReal(lRhoNormal, lTNormal);
  bModNormal_ = dependsRhoT_.bMod.interpToReal(lRhoNormal, lTNormal);
  dPdENormal_ = dependsRhoT_.dPdE.interpToReal(lRhoNormal, lTNormal);
  Real dPdR = dependsRhoT_.dPdRho.interpToReal(lRhoNormal, lTNormal);
  dVdTNormal_ = dPdENormal_ * CvNormal_ / (rhoNormal_ * rhoNormal_ * dPdR);

  return status;
}

inline void SpinerEOSDependsRhoSie::calcBMod_(SP5Tables &tables) {
  for (int j = 0; j < tables.bMod.dim(2); j++) {
    Real lRho = tables.bMod.range(1).x(j);
    Real rho = fromLog_(lRho, lRhoOffset_);
    for (int i = 0; i < tables.bMod.dim(1); i++) {
      Real press = tables.P(j, i);
      Real DPDR_E = tables.dPdRho(j, i);
      Real DPDE_R = tables.dPdE(j, i);
      Real DEDR_T = tables.dEdRho(j, i);
      Real DPDR_T = DPDR_E + DPDE_R * DEDR_T;
      Real bMod;
      if (DPDE_R > 0.0 && rho > 0.0) {
        bMod = rho * DPDR_E + DPDE_R * (press / rho);
      } else if (rho > 0.0) {
        bMod = std::max(rho * DPDR_T, 0.0);
      } else {
        bMod = 0.0;
      }
      tables.bMod(j, i) = std::max(bMod, robust::EPS());
    }
  }
}

inline SpinerEOSDependsRhoSie SpinerEOSDependsRhoSie::GetOnDevice() {
  return STricks::GetOnDevice(this);
}

void SpinerEOSDependsRhoSie::Finalize() { STricks::Finalize(this); }

inline std::size_t SpinerEOSDependsRhoSie::DynamicMemorySizeInBytes() const {
  return STricks::DynamicMemorySizeInBytes(this);
}

inline std::size_t SpinerEOSDependsRhoSie::DumpDynamicMemory(char *dst) {
  return STricks::DumpDynamicMemory(dst, this);
}

inline std::size_t
SpinerEOSDependsRhoSie::SetDynamicMemory(char *src, const SharedMemSettings &stngs) {
  if (stngs.data != nullptr) src = stngs.data;
  return STricks::SetDynamicMemory(src, this);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real
SpinerEOSDependsRhoSie::TemperatureFromDensityInternalEnergy(const Real rho,
                                                             const Real sie,
                                                             Indexer_t &&lambda) const {
  return interpRhoSie_(rho, sie, T_, lambda);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real
SpinerEOSDependsRhoSie::InternalEnergyFromDensityTemperature(const Real rho, const Real T,
                                                             Indexer_t &&lambda) const {
  return interpRhoT_(rho, T, sie_, lambda);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoSie::PressureFromDensityTemperature(
    const Real rho, const Real T, Indexer_t &&lambda) const {
  return interpRhoT_(rho, T, dependsRhoT_.P, lambda);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoSie::PressureFromDensityInternalEnergy(
    const Real rho, const Real sie, Indexer_t &&lambda) const {
  return interpRhoSie_(rho, sie, dependsRhoSie_.P, lambda);
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoSie::MinInternalEnergyFromDensity(
    const Real rho, Indexer_t &&lambda) const {
  MinInternalEnergyIsNotEnabled("SpinerEOSDependsRhoSie");
  return 0.0;
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoSie::EntropyFromDensityTemperature(
    const Real rho, const Real temperature, Indexer_t &&lambda) const {
  EntropyIsNotEnabled("SpinerEOSDependsRhoSie");
  return 1.0;
}
template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoSie::EntropyFromDensityInternalEnergy(
    const Real rho, const Real sie, Indexer_t &&lambda) const {
  EntropyIsNotEnabled("SpinerEOSDependsRhoSie");
  return 1.0;
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoSie::SpecificHeatFromDensityTemperature(
    const Real rho, const Real T, Indexer_t &&lambda) const {
  return 1. / interpRhoT_(rho, T, dependsRhoT_.dTdE, lambda);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real
SpinerEOSDependsRhoSie::SpecificHeatFromDensityInternalEnergy(const Real rho,
                                                              const Real sie,
                                                              Indexer_t &&lambda) const {
  return 1. / interpRhoSie_(rho, sie, dependsRhoSie_.dTdE, lambda);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoSie::BulkModulusFromDensityTemperature(
    const Real rho, const Real T, Indexer_t &&lambda) const {
  return interpRhoT_(rho, T, dependsRhoT_.bMod, lambda);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real
SpinerEOSDependsRhoSie::BulkModulusFromDensityInternalEnergy(const Real rho,
                                                             const Real sie,
                                                             Indexer_t &&lambda) const {
  return interpRhoSie_(rho, sie, dependsRhoSie_.bMod, lambda);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real
SpinerEOSDependsRhoSie::GruneisenParamFromDensityTemperature(const Real rho, const Real T,
                                                             Indexer_t &&lambda) const {
  const Real dpde = interpRhoT_(rho, T, dependsRhoT_.dPdE, lambda);
  return dpde / rho;
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real
SpinerEOSDependsRhoSie::GruneisenParamFromDensityInternalEnergy(
    const Real rho, const Real sie, Indexer_t &&lambda) const {
  const Real lRho = toLog_(rho, lRhoOffset_);
  const Real lE = toLog_(sie, lEOffset_);
  const Real dpde = dependsRhoSie_.dPdE.interpToReal(lRho, lE);
  return dpde / rho;
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void
SpinerEOSDependsRhoSie::DensityEnergyFromPressureTemperature(const Real press,
                                                             const Real temp,
                                                             Indexer_t &&lambda,
                                                             Real &rho, Real &sie) const {
  Real lT = toLog_(temp, lTOffset_);
  Real lRho = lRhoFromPlT_(press, lT, lambda);
  rho = fromLog_(lRho, lRhoOffset_);
  sie = sie_.interpToReal(lRho, lT);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void
SpinerEOSDependsRhoSie::FillEos(Real &rho, Real &temp, Real &energy, Real &press,
                                Real &cv, Real &bmod, const unsigned long output,
                                Indexer_t &&lambda) const {
  Real lRho, lT, lE;
  if (output == thermalqs::none) {
    UNDEFINED_ERROR;
  }
  if (output & thermalqs::temperature && output & thermalqs::specific_internal_energy) {
    UNDEFINED_ERROR;
  }
  if (output & thermalqs::density) {
    if (!(output & thermalqs::pressure || output & thermalqs::temperature)) {
      lT = toLog_(temp, lTOffset_);
      lRho = lRhoFromPlT_(press, lT, lambda);
      rho = fromLog_(lRho, lRhoOffset_);
    } else {
      UNDEFINED_ERROR;
    }
  } else {
    lRho = toLog_(rho, lRhoOffset_);
    if (!variadic_utils::is_nullptr(lambda)) {
      IndexerUtils::Get<IndexableTypes::LogDensity>(lambda, Lambda::lRho) = lRho;
    }
  }
  if (output & thermalqs::temperature) {
    lE = toLog_(energy, lEOffset_);
    temp = T_.interpToReal(lRho, lE);
    if (output & thermalqs::pressure) {
      press = dependsRhoSie_.P.interpToReal(lRho, lE);
    }
    if (output & thermalqs::specific_heat) {
      cv = 1. / dependsRhoSie_.dTdE.interpToReal(lRho, lE);
    }
    if (output & thermalqs::bulk_modulus) {
      bmod = dependsRhoSie_.bMod.interpToReal(lRho, lE);
    }
  }
  if (output & thermalqs::specific_internal_energy) {
    lT = toLog_(temp, lTOffset_);
    energy = sie_.interpToReal(lRho, lT);
    if (output & thermalqs::pressure) {
      press = dependsRhoT_.P.interpToReal(lRho, lT);
    }
    if (output & thermalqs::specific_heat) {
      cv = 1. / dependsRhoT_.dTdE.interpToReal(lRho, lT);
    }
    if (output & thermalqs::bulk_modulus) {
      bmod = dependsRhoT_.bMod.interpToReal(lRho, lT);
    }
  }
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void SpinerEOSDependsRhoSie::ValuesAtReferenceState(
    Real &rho, Real &temp, Real &sie, Real &press, Real &cv, Real &bmod, Real &dpde,
    Real &dvdt, Indexer_t &&lambda) const {
  rho = rhoNormal_;
  temp = TNormal_;
  sie = sieNormal_;
  press = PNormal_;
  cv = CvNormal_;
  bmod = bModNormal_;
  dpde = dPdENormal_;
  dvdt = dVdTNormal_;
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoSie::interpRhoT_(
    const Real rho, const Real T, const DataBox &db, Indexer_t &&lambda) const {
  const Real lRho = toLog_(rho, lRhoOffset_);
  const Real lT = toLog_(T, lTOffset_);
  if (!variadic_utils::is_nullptr(lambda)) {
    IndexerUtils::Get<IndexableTypes::LogDensity>(lambda, Lambda::lRho) = lRho;
  }
  return db.interpToReal(lRho, lT);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoSie::interpRhoSie_(
    const Real rho, const Real sie, const DataBox &db, Indexer_t &&lambda) const {
  const Real lRho = toLog_(rho, lRhoOffset_);
  const Real lE = toLog_(sie, lEOffset_);
  if (!variadic_utils::is_nullptr(lambda)) {
    IndexerUtils::Get<IndexableTypes::LogDensity>(lambda, Lambda::lRho) = lRho;
  }
  return db.interpToReal(lRho, lE);
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real SpinerEOSDependsRhoSie::lRhoFromPlT_(
    const Real P, const Real lT, Indexer_t &&lambda) const {
  const RootFinding1D::RootCounts *pcounts =
      (memoryStatus_ == DataStatus::OnDevice) ? nullptr : &counts;
  RootFinding1D::Status status = RootFinding1D::Status::SUCCESS;
  Real lRho;
  Real dPdRhoMax = dPdRhoMax_.interpToReal(lT);
  Real PMax = PlRhoMax_.interpToReal(lT);
  if (dPdRhoMax > 0 && P > PMax) {
    Real rho = (P - PMax) / dPdRhoMax + rhoMax_;
    lRho = toLog_(rho, lRhoOffset_);
    if (pcounts != nullptr) {
      pcounts->increment(0);
    }
  } else {
    Real lRhoGuess = reproducible_ ? lRhoMin_ : 0.5 * (lRhoMin_ + lRhoMax_);
    if (!variadic_utils::is_nullptr(lambda)) {
      Real lRho_cache =
          IndexerUtils::Get<IndexableTypes::LogDensity>(lambda, Lambda::lRho);
      if ((lRhoMin_ <= lRho_cache) && (lRho_cache <= lRhoMax_)) {
        lRhoGuess = lRho_cache;
      }
    }
    const callable_interp::l_interp PFunc(dependsRhoT_.P, lT);
    status = ROOT_FINDER(PFunc, P, lRhoGuess, lRhoMin_, lRhoMax_, robust::EPS(),
                         robust::EPS(), lRho, pcounts);
    if (status != RootFinding1D::Status::SUCCESS) {
#if SPINER_EOS_VERBOSE
      std::stringstream errorMessage;
      errorMessage << "inverting P table for logRho failed\n"
                   << "matid     = " << matid_ << "\n"
                   << "lT        = " << lT << "\n"
                   << "P         = " << P << "\n"
                   << "lRhoGuess = " << lRhoGuess << std::endl;
      EOS_ERROR(errorMessage.str().c_str());
#endif // SPINER_EOS_VERBOSE
      lRho = reproducible_ ? lRhoMin_ : lRhoGuess;
    }
  }
  if (!variadic_utils::is_nullptr(lambda)) {
    IndexerUtils::Get<IndexableTypes::LogDensity>(lambda, Lambda::lRho) = lRho;
    if constexpr (variadic_utils::is_indexable_v<Indexer_t, IndexableTypes::RootStatus>) {
      lambda[IndexableTypes::RootStatus()] = static_cast<Real>(status);
    }
  }
  return lRho;
}

} // namespace singularity

#endif // SINGULARITY_USE_SPINER_WITH_HDF5
#endif // _SINGULARITY_EOS_EOS_EOS_SPINER_HPP_
