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

#ifndef _SINGULARITY_EOS_EOS_EOS_SPINER_RHO_SIE_HPP_
#define _SINGULARITY_EOS_EOS_EOS_SPINER_RHO_SIE_HPP_

#ifdef SINGULARITY_USE_SPINER_WITH_HDF5

#include <type_traits>

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
#include <singularity-eos/eos/eos_spiner_common.hpp>

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
  static constexpr int NGRIDS = 3;
  using Grid_t = Spiner::PiecewiseGrid1D<Real, NGRIDS>;
  using DataBox = Spiner::DataBox<Real, Grid_t>;

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
  int matid() const { return matid_; }
  PORTABLE_FORCEINLINE_FUNCTION Real lRhoOffset() const { return lRhoOffset_; }
  PORTABLE_FORCEINLINE_FUNCTION Real lTOffset() const { return lTOffset_; }
  PORTABLE_FORCEINLINE_FUNCTION Real lEOffset() const { return lEOffset_; }
  PORTABLE_FORCEINLINE_FUNCTION Real rhoMin() const {
    return fromLog_(lRhoMin_, lRhoOffset_);
  }
  PORTABLE_FORCEINLINE_FUNCTION Real rhoMax() const { return rhoMax_; }
  PORTABLE_FORCEINLINE_FUNCTION Real TMin() const {
    return fromLog_(sie_.range(0).min(), lTOffset_);
  }
  PORTABLE_FORCEINLINE_FUNCTION Real TMax() const {
    return fromLog_(sie_.range(0).max(), lTOffset_);
  }
  PORTABLE_FORCEINLINE_FUNCTION Real sieMin() const {
    return fromLog_(T_.range(0).min(), lEOffset_);
  }
  PORTABLE_FORCEINLINE_FUNCTION Real sieMax() const {
    return fromLog_(T_.range(0).max(), lEOffset_);
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

  constexpr static inline int nlambda() noexcept { return _n_lambda; }
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
                               hid_t lEGroup, hid_t coldGroupd);
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
  DataBox PCold_, sieCold_, bModCold_;
  DataBox dPdRhoCold_, dPdECold_, dTdRhoCold_, dTdECold_, dEdTCold_;
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
      &PlRhoMax_, &dPdRhoMax_, &PCold_, &sieCold_, &bModCold_, &dPdRhoCold_, &dPdECold_, \
      &dTdRhoCold_, &dTdECold_, &dEdTCold_,
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

inline SpinerEOSDependsRhoSie::SpinerEOSDependsRhoSie(const std::string &filename,
                                                      const std::string &materialName,
                                                      TableSplit split,
                                                      bool reproducibility_mode)
    : split_(split), reproducible_(reproducibility_mode),
      memoryStatus_(DataStatus::OnHost) {

  std::string matid_str;
  hid_t file, matGroup, lTGroup, lEGroup, coldGroup;
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
  coldGroup = H5Gopen(matGroup, SP5::Depends::coldCurve, H5P_DEFAULT);

  status +=
      H5LTget_attribute_int(file, materialName.c_str(), SP5::Material::matid, &matid_);
  matid_str = std::to_string(matid_);

  status += loadDataboxes_(matid_str, file, lTGroup, lEGroup, coldGroup);

  status += H5Gclose(lTGroup);
  status += H5Gclose(lEGroup);
  status += H5Gclose(matGroup);
  status += H5Fclose(file);
  status += H5Gclose(coldGroup);

  if (status != H5_SUCCESS) {
    EOS_ERROR("SpinerDependsRhoSie: HDF5 error\n");
  }
}

herr_t SpinerEOSDependsRhoSie::loadDataboxes_(const std::string &matid_str, hid_t file,
                                              hid_t lTGroup, hid_t lEGroup,
                                              hid_t coldGroup) {
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

  // cold curves
  status += PCold_.loadHDF(coldGroup, SP5::Fields::P);
  status += sieCold_.loadHDF(coldGroup, SP5::Fields::sie);
  status += bModCold_.loadHDF(coldGroup, SP5::Fields::bMod);
  status += dPdRhoCold_.loadHDF(coldGroup, SP5::Fields::dPdRho);

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
  Real lRho = toLog_(rho, lRhoOffset_);
  return sieCold_.interpToReal(lRho);
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

inline SpinerEOSDependsRhoSie::SpinerEOSDependsRhoSie(const std::string &filename,
                                                      int matid, TableSplit split,
                                                      bool reproducibility_mode)
    : matid_(matid), split_(split), reproducible_(reproducibility_mode),
      memoryStatus_(DataStatus::OnHost) {

  std::string matid_str = std::to_string(matid);
  hid_t file, matGroup, lTGroup, lEGroup, coldGroup;
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
  coldGroup = H5Gopen(matGroup, SP5::Depends::coldCurve, H5P_DEFAULT);

  status += loadDataboxes_(matid_str, file, lTGroup, lEGroup, coldGroup);

  status += H5Gclose(lTGroup);
  status += H5Gclose(lEGroup);
  status += H5Gclose(matGroup);
  status += H5Fclose(file);
  status += H5Fclose(coldGroup);

  if (status != H5_SUCCESS) {
    EOS_ERROR("SpinerDependsRhoSIE: HDF5 error\n");
  }
}

} // namespace singularity

#endif // SINGULARITY_USE_SPINER_WITH_HDF5
#endif // _SINGULARITY_EOS_EOS_EOS_SPINER_RHO_SIE_HPP_
