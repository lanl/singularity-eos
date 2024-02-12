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

#ifndef _SINGULARITY_EOS_EOS_EOS_STELLAR_COLLAPSE_HPP_
#define _SINGULARITY_EOS_EOS_EOS_STELLAR_COLLAPSE_HPP_
#include <type_traits>
#ifdef SINGULARITY_USE_SPINER_WITH_HDF5

// C++ includes
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// C includes
#include <cstdlib>
#include <hdf5.h>
#include <hdf5_hl.h>

// ports-of-call
#include <ports-of-call/portability.hpp>

// singularity-eos
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/fast-math/logs.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>
#include <singularity-eos/base/sp5/singularity_eos_sp5.hpp>
#include <singularity-eos/eos/eos_base.hpp>

// spiner
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/sp5.hpp>

#define STELLAR_COLLAPSE_EOS_VERBOSE (0)

namespace singularity {
using namespace eos_base;

// Note the Stellar Collapse tables have units of:
// 1. Ye (unitless)
// 2. log(MeV) for temperature
// 3. log(g/cm^3) for density.
//
// TODO(JMM): For now the bottom of the table is a floor and the top
// is linear extrapolation in log-log space. We should reconsider this
// and introduce extrapolation as needed.
class StellarCollapse : public EosBase<StellarCollapse> {
 public:
  using DataBox = Spiner::DataBox<Real>;
  using Grid_t = Spiner::RegularGrid1D<Real>;

  // A weakly typed index map for lambdas
  struct Lambda {
    enum Index { Ye = 0, lT = 1 };
  };

  // Generic functions provided by the base class. These contain
  // e.g. the vector overloads that use the scalar versions declared
  // here We explicitly list, rather than using the macro because we
  // overload some methods.
  using EosBase<StellarCollapse>::TemperatureFromDensityInternalEnergy;
  using EosBase<StellarCollapse>::InternalEnergyFromDensityTemperature;
  using EosBase<StellarCollapse>::PressureFromDensityTemperature;
  using EosBase<StellarCollapse>::PressureFromDensityInternalEnergy;
  using EosBase<StellarCollapse>::MinInternalEnergyFromDensity;
  using EosBase<StellarCollapse>::EntropyFromDensityTemperature;
  using EosBase<StellarCollapse>::EntropyFromDensityInternalEnergy;
  using EosBase<StellarCollapse>::SpecificHeatFromDensityTemperature;
  using EosBase<StellarCollapse>::SpecificHeatFromDensityInternalEnergy;
  using EosBase<StellarCollapse>::BulkModulusFromDensityTemperature;
  using EosBase<StellarCollapse>::BulkModulusFromDensityInternalEnergy;
  using EosBase<StellarCollapse>::GruneisenParamFromDensityTemperature;
  using EosBase<StellarCollapse>::GruneisenParamFromDensityInternalEnergy;
  using EosBase<StellarCollapse>::FillEos;

  inline StellarCollapse(const std::string &filename, bool use_sp5 = false,
                         bool filter_bmod = true);

  // Saves to an SP5 file
  inline void Save(const std::string &filename);

  PORTABLE_INLINE_FUNCTION
  StellarCollapse() : memoryStatus_(DataStatus::Deallocated) {}

  inline StellarCollapse GetOnDevice();

  PORTABLE_INLINE_FUNCTION
  Real TemperatureFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION
  Real InternalEnergyFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION
  Real PressureFromDensityTemperature(const Real rho, const Real temperature,
                                      Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION
  Real PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                         Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION
  Real MinInternalEnergyFromDensity(const Real rho, Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION
  Real EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                     Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION
  Real EntropyFromDensityInternalEnergy(const Real rho, const Real sie,
                                        Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION
  Real SpecificHeatFromDensityTemperature(const Real rho, const Real temperature,
                                          Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION
  Real SpecificHeatFromDensityInternalEnergy(const Real rho, const Real sie,
                                             Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION
  Real BulkModulusFromDensityTemperature(const Real rho, const Real temperature,
                                         Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION
  Real BulkModulusFromDensityInternalEnergy(const Real rho, const Real sie,
                                            Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION
  Real GruneisenParamFromDensityTemperature(const Real rho, const Real temperature,
                                            Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION
  Real GruneisenParamFromDensityInternalEnergy(const Real rho, const Real sie,
                                               Real *lambda = nullptr) const;
  PORTABLE_INLINE_FUNCTION
  void DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                            Real *lambda, Real &rho, Real &sie) const;

  // Properties of an NSE EOS
  PORTABLE_INLINE_FUNCTION
  void MassFractionsFromDensityTemperature(const Real rho, const Real temperature,
                                           Real &Xa, Real &Xh, Real &Xn, Real &Xp,
                                           Real &Abar, Real &Zbar,
                                           Real *lambda = nullptr) const;

  PORTABLE_INLINE_FUNCTION
  void FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
               const unsigned long output, Real *lambda = nullptr) const;

  PORTABLE_INLINE_FUNCTION
  void ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                              Real &bmod, Real &dpde, Real &dvdt,
                              Real *lambda = nullptr) const;
  // Generic functions provided by the base class. These contain e.g. the vector
  // overloads that use the scalar versions declared here
  static constexpr unsigned long PreferredInput() { return _preferred_input; }
  static inline unsigned long scratch_size(std::string method, unsigned int nelements) {
    return 0;
  }
  static inline unsigned long max_scratch_size(unsigned int nelements) { return 0; }
  PORTABLE_FORCEINLINE_FUNCTION Real rhoMin() const { return rho_(lRhoMin_); }
  PORTABLE_FORCEINLINE_FUNCTION Real rhoMax() const { return rho_(lRhoMax_); }
  PORTABLE_FORCEINLINE_FUNCTION Real TMin() const { return T_(lTMin_); }
  PORTABLE_FORCEINLINE_FUNCTION Real TMax() const { return T_(lTMax_); }
  PORTABLE_FORCEINLINE_FUNCTION Real YeMin() const { return YeMin_; }
  PORTABLE_FORCEINLINE_FUNCTION Real YeMax() const { return YeMax_; }
  PORTABLE_FORCEINLINE_FUNCTION Real sieMin() const { return sieMin_; }
  PORTABLE_FORCEINLINE_FUNCTION Real sieMax() const { return sieMax_; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    printf("StellarCollapse parameters:\n"
           "depends on log10(rho), log10(T), Ye\n"
           "lrho bounds = %.14e, %.14e\n"
           "lT bounds = %.14e, %.14e\n"
           "Ye bounds = %.14e, %.14e\n",
           lRhoMin_, lRhoMax_, lTMin_, lTMax_, YeMin_, YeMax_);
    return;
  }
  PORTABLE_FORCEINLINE_FUNCTION Real MinimumDensity() const { return rhoMin(); }
  PORTABLE_FORCEINLINE_FUNCTION Real MinimumTemperature() const { return TMin(); }
  PORTABLE_INLINE_FUNCTION
  int nlambda() const noexcept { return _n_lambda; }
  inline RootFinding1D::Status rootStatus() const { return status_; }
  RootFinding1D::RootCounts counts;
  inline void Finalize();
  static std::string EosType() { return std::string("StellarCollapse"); }
  static std::string EosPyType() { return EosType(); }

  // A utility function for working with 3D DataBoxes that converts
  // from log10 to fastlog 10. Mostly for unit testing.

  // TODO(JMM) Should this be in a utilities function somewhere? In
  // the math folder or something? 3D is pretty specific to Stellar
  // Collapse, so I think we can leave it here for now?
  inline static void dataBoxToFastLogs(DataBox &db, DataBox &scratch,
                                       bool dependent_var_log);

 private:
  inline void LoadFromSP5File_(const std::string &filename);
  inline void LoadFromStellarCollapseFile_(const std::string &filename, bool filter_bmod);
  inline int readSCInt_(const hid_t &file_id, const std::string &name);
  inline void readBounds_(const hid_t &file_id, const std::string &name, int size,
                          Real &lo, Real &hi);
  inline void readSCDset_(const hid_t &file_id, const std::string &name,
                          const Grid_t &Ye_grid, const Grid_t &lT_grid,
                          const Grid_t &lRho_grid, DataBox &db);

  inline void medianFilter_(DataBox &db);
  inline void medianFilter_(const DataBox &in, DataBox &out);
  inline void fillMedianBuffer_(Real buffer[], int width, int iY, int iT, int irho,
                                const DataBox &tab) const;
  inline Real findMedian_(Real buffer[], int size) const;
  inline void computeBulkModulus_();
  inline void computeColdAndHotCurves_();
  inline void setNormalValues_();

  PORTABLE_FORCEINLINE_FUNCTION void checkLambda_(Real *lambda) const noexcept {
    if (lambda == nullptr) {
      EOS_ERROR("StellarCollapse: lambda must contain Ye and 1 space for caching.\n");
    }
  }

  PORTABLE_FORCEINLINE_FUNCTION Real toLog_(const Real x,
                                            const Real offset) const noexcept {
    // StellarCollapse can't use fast logs, unless we re-grid onto the
    // "fast log grid"
    return FastMath::log10(std::abs(std::max(x, -offset) + offset) + robust::SMALL());
  }
  PORTABLE_FORCEINLINE_FUNCTION Real fromLog_(const Real lx,
                                              const Real offset) const noexcept {
    // StellarCollapse can't use fast logs, unless we re-grid onto the
    // "fast log grid"
    // return std::pow(10., lx) - offset;
    return FastMath::pow10(lx) - offset;
  }
  PORTABLE_FORCEINLINE_FUNCTION Real lRho_(const Real rho) const noexcept {
    Real out = toLog_(rho, lRhoOffset_);
    return out;
    // return out < lRhoMin_ ? lRhoMin_ : out;
  }
  PORTABLE_FORCEINLINE_FUNCTION Real lT_(const Real T) const noexcept {
    return toLog_(T, lTOffset_);
  }
  PORTABLE_FORCEINLINE_FUNCTION Real rho_(const Real lRho) const noexcept {
    Real rho = fromLog_(lRho, lRhoOffset_);
    return rho < 0 ? 0 : rho;
  }
  PORTABLE_FORCEINLINE_FUNCTION Real T_(const Real lT) const noexcept {
    return fromLog_(lT, lTOffset_);
  }
  PORTABLE_FORCEINLINE_FUNCTION Real le2e_(const Real le) const noexcept {
    return fromLog_(le, lEOffset_);
  }
  PORTABLE_FORCEINLINE_FUNCTION Real e2le_(const Real e) const noexcept {
    return toLog_(e, lEOffset_);
  }
  PORTABLE_FORCEINLINE_FUNCTION Real lP2P_(const Real lP) const noexcept {
    return fromLog_(lP, lPOffset_);
  }
  PORTABLE_FORCEINLINE_FUNCTION Real P2lP_(const Real P) const noexcept {
    return toLog_(P, lPOffset_);
  }
  PORTABLE_FORCEINLINE_FUNCTION Real lB2B_(const Real lB) const noexcept {
    return fromLog_(lB, lBOffset_);
  }
  PORTABLE_FORCEINLINE_FUNCTION Real B2lB_(const Real B) const noexcept {
    return toLog_(B, lBOffset_);
  }

  PORTABLE_INLINE_FUNCTION Real lTFromlRhoSie_(const Real lRho, const Real sie,
                                               Real *lambda) const noexcept;
  PORTABLE_INLINE_FUNCTION __attribute__((always_inline)) void
  getLogsFromRhoT_(const Real rho, const Real temp, Real *lambda, Real &lRho, Real &lT,
                   Real &Ye) const noexcept {
    checkLambda_(lambda);
    lRho = lRho_(rho);
    lT = lT_(temp);
    Ye = lambda[Lambda::Ye];
    lambda[Lambda::lT] = lT;
  }
  PORTABLE_INLINE_FUNCTION __attribute__((always_inline)) void
  getLogsFromRhoSie_(const Real rho, const Real sie, Real *lambda, Real &lRho, Real &lT,
                     Real &Ye) const noexcept {
    lRho = lRho_(rho);
    lT = lTFromlRhoSie_(lRho, sie, lambda);
    Ye = lambda[Lambda::Ye];
    return;
  }

  static constexpr const unsigned long _preferred_input =
      thermalqs::density | thermalqs::temperature;

  // Dependent variables
  DataBox lP_, lE_, dPdRho_, dPdE_, dEdT_, lBMod_;
  DataBox entropy_; // kb/baryon
  DataBox Xa_;      // mass fraction of alpha particles
  DataBox Xh_;      // mass fraction of heavy ions
  DataBox Xn_;      // mass fraction of neutrons
  DataBox Xp_;      // mass fraction of protons
  DataBox Abar_;    // Average atomic mass
  DataBox Zbar_;    // Average atomic number
  // Spiner::DataBox gamma_; // polytropic index. dlog(P)/dlog(rho).
  // dTdRho_, dTdE_, dEdRho_, dEdT_;

  // Bounds of dependent variables. Needed for root finding.
  DataBox eCold_, eHot_;

  // Independent variable bounds
  int numRho_, numT_, numYe_;
  Real lRhoMin_, lRhoMax_;
  Real lTMin_, lTMax_;
  Real YeMin_, YeMax_;
  Real sieMin_, sieMax_;

  static constexpr Real MeV2GK_ = 11.604525006;
  static constexpr Real GK2MeV_ = 1. / MeV2GK_;
  static constexpr Real MeV2K_ = 1.e9 * MeV2GK_;
  static constexpr Real K2MeV_ = 1. / MeV2K_;
  static constexpr Real TNormal_ = 5 * GK2MeV_; // Threshold of NSE
  static constexpr Real rhoNormal_ = 2.e12;     // 1./100'th of nuclear density
  static constexpr Real YeNormal_ = 0.3517;     // Beta equilibrium value
  Real sieNormal_, PNormal_, SNormal_;
  Real CvNormal_, bModNormal_, dPdENormal_, dVdTNormal_;

  // offsets must be non-negative
  Real lEOffset_;
  static constexpr Real lRhoOffset_ = 0.0; // TODO(JMM): Address if this ever changes
  static constexpr Real lTOffset_ = 0.0;
  static constexpr Real lPOffset_ = 0.0;
  static constexpr Real lBOffset_ = 0.0;

  // whereAmI_ and status_ used only for reporting. They are not thread-safe.
  mutable RootFinding1D::Status status_ = RootFinding1D::Status::SUCCESS;
  static constexpr const Real ROOT_THRESH = 1e-14; // TODO: experiment
  DataStatus memoryStatus_ = DataStatus::Deallocated;
  static constexpr const int _n_lambda = 2;
  static constexpr const char *_lambda_names[] = {"Ye", "log(T)"};

  // Stuff for median filter smoothing
  static constexpr Real DELTASMOOTH = 10.0;
  static constexpr int MF_W = 3;
  static constexpr int MF_S = (2 * MF_W + 1) * (2 * MF_W + 1) * (2 * MF_W + 1);
};

// ======================================================================
// Implementation details below
// ======================================================================

namespace callable_interp {

class LogT {
 public:
  using DataBox = Spiner::DataBox<Real>;
  PORTABLE_INLINE_FUNCTION
  LogT(const DataBox &field, const Real Ye, const Real lRho)
      : field_(field), Ye_(Ye), lRho_(lRho) {}
  PORTABLE_INLINE_FUNCTION Real operator()(const Real lT) const {
    return field_.interpToReal(Ye_, lT, lRho_);
  }

 private:
  const DataBox &field_;
  const Real Ye_, lRho_;
};

} // namespace callable_interp

// For some reason, the linker doesn't like this being a member field
// of StellarCollapse.  So we'll make it a global variable.
constexpr char METADATA_NAME[] = "Metadata";

inline StellarCollapse::StellarCollapse(const std::string &filename, bool use_sp5,
                                        bool filter_bmod) {
  if (use_sp5) {
    LoadFromSP5File_(filename);
  } else {
    LoadFromStellarCollapseFile_(filename, filter_bmod);
  }
  setNormalValues_();
}

// Saves to an SP5 file
inline void StellarCollapse::Save(const std::string &filename) {
  herr_t status = H5_SUCCESS;
  hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // Metadata
  hid_t metadata = H5Gcreate(file, METADATA_NAME, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status += H5LTset_attribute_string(file, METADATA_NAME, SP5::Offsets::messageName,
                                     SP5::Offsets::message);
  status +=
      H5LTset_attribute_double(file, METADATA_NAME, SP5::Offsets::sie, &lEOffset_, 1);
  H5Gclose(metadata);

  // Databoxes
  status += lP_.saveHDF(file, "logpress");
  status += lE_.saveHDF(file, "logenergy");
  status += dPdRho_.saveHDF(file, "dpdrhoe");
  status += dPdE_.saveHDF(file, "dpderho");
  status += dEdT_.saveHDF(file, "dedt");
  status += entropy_.saveHDF(file, "entropy");
  status += Xa_.saveHDF(file, "Xa");
  status += Xh_.saveHDF(file, "Xh");
  status += Xn_.saveHDF(file, "Xn");
  status += Xp_.saveHDF(file, "Xp");
  status += Abar_.saveHDF(file, "Abar");
  status += Zbar_.saveHDF(file, "Zbar");
  status += lBMod_.saveHDF(file, "logbulkmodulus");
  status += eCold_.saveHDF(file, "ecold");
  status += eHot_.saveHDF(file, "ehot");

  status += H5Fclose(file);
  if (status != H5_SUCCESS) {
    EOS_ERROR("[StellarCollapse::Save]: There was a problem with HDF5\n");
  }
}

inline StellarCollapse StellarCollapse::GetOnDevice() {
  StellarCollapse other;
  other.lP_ = Spiner::getOnDeviceDataBox<Real>(lP_);
  other.lE_ = Spiner::getOnDeviceDataBox<Real>(lE_);
  other.dPdRho_ = Spiner::getOnDeviceDataBox<Real>(dPdRho_);
  other.dPdE_ = Spiner::getOnDeviceDataBox<Real>(dPdE_);
  other.dEdT_ = Spiner::getOnDeviceDataBox<Real>(dEdT_);
  other.entropy_ = Spiner::getOnDeviceDataBox<Real>(entropy_);
  other.Xa_ = Spiner::getOnDeviceDataBox<Real>(Xa_);
  other.Xh_ = Spiner::getOnDeviceDataBox<Real>(Xh_);
  other.Xn_ = Spiner::getOnDeviceDataBox<Real>(Xn_);
  other.Xp_ = Spiner::getOnDeviceDataBox<Real>(Xp_);
  other.Abar_ = Spiner::getOnDeviceDataBox<Real>(Abar_);
  other.Zbar_ = Spiner::getOnDeviceDataBox<Real>(Zbar_);
  other.lBMod_ = Spiner::getOnDeviceDataBox<Real>(lBMod_);
  other.eCold_ = Spiner::getOnDeviceDataBox<Real>(eCold_);
  other.eHot_ = Spiner::getOnDeviceDataBox<Real>(eHot_);
  other.memoryStatus_ = DataStatus::OnDevice;
  other.numRho_ = numRho_;
  other.numT_ = numT_;
  other.numYe_ = numYe_;
  other.lTMin_ = lTMin_;
  other.lTMax_ = lTMax_;
  other.YeMin_ = YeMin_;
  other.YeMax_ = YeMax_;
  other.sieMin_ = sieMin_;
  other.sieMax_ = sieMax_;
  other.lEOffset_ = lEOffset_;
  other.sieNormal_ = sieNormal_;
  other.PNormal_ = PNormal_;
  other.SNormal_ = SNormal_;
  other.CvNormal_ = CvNormal_;
  other.bModNormal_ = bModNormal_;
  other.dPdENormal_ = dPdENormal_;
  other.dVdTNormal_ = dVdTNormal_;
  other.status_ = status_;
  return other;
}

inline void StellarCollapse::Finalize() {
  lP_.finalize();
  lE_.finalize();
  dPdRho_.finalize();
  dPdE_.finalize();
  dEdT_.finalize();
  entropy_.finalize();
  Xa_.finalize();
  Xh_.finalize();
  Xn_.finalize();
  Xp_.finalize();
  Abar_.finalize();
  Zbar_.finalize();
  lBMod_.finalize();
  eCold_.finalize();
  eHot_.finalize();
  memoryStatus_ = DataStatus::Deallocated;
}

PORTABLE_INLINE_FUNCTION
Real StellarCollapse::TemperatureFromDensityInternalEnergy(const Real rho, const Real sie,
                                                           Real *lambda) const {
  const Real lRho = lRho_(rho);
  const Real lT = lTFromlRhoSie_(lRho, sie, lambda);
  return T_(lT);
}

PORTABLE_INLINE_FUNCTION
Real StellarCollapse::InternalEnergyFromDensityTemperature(const Real rho,
                                                           const Real temp,
                                                           Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoT_(rho, temp, lambda, lRho, lT, Ye);
  const Real lE = lE_.interpToReal(Ye, lT, lRho);
  return le2e_(lE);
}

PORTABLE_INLINE_FUNCTION
Real StellarCollapse::PressureFromDensityTemperature(const Real rho,
                                                     const Real temperature,
                                                     Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoT_(rho, temperature, lambda, lRho, lT, Ye);
  const Real lP = lP_.interpToReal(Ye, lT, lRho);
  return lP2P_(lP);
}

PORTABLE_INLINE_FUNCTION
Real StellarCollapse::PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                                        Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoSie_(rho, sie, lambda, lRho, lT, Ye);
  const Real lP = lP_.interpToReal(Ye, lT, lRho);
  return lP2P_(lP);
}
PORTABLE_INLINE_FUNCTION
Real StellarCollapse::MinInternalEnergyFromDensity(const Real rho, Real *lambda) const {
  MinInternalEnergyIsNotEnabled("Stellar Collapse");
  return 0.0;
}

PORTABLE_INLINE_FUNCTION
Real StellarCollapse::EntropyFromDensityTemperature(const Real rho,
                                                    const Real temperature,
                                                    Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoT_(rho, temperature, lambda, lRho, lT, Ye);
  const Real entropy = entropy_.interpToReal(Ye, lT, lRho);
  return (entropy > robust::EPS() ? entropy : robust::EPS());
}

PORTABLE_INLINE_FUNCTION
Real StellarCollapse::EntropyFromDensityInternalEnergy(const Real rho, const Real sie,
                                                       Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoSie_(rho, sie, lambda, lRho, lT, Ye);
  const Real entropy = entropy_.interpToReal(Ye, lT, lRho);
  return (entropy > robust::EPS() ? entropy : robust::EPS());
}

PORTABLE_INLINE_FUNCTION
Real StellarCollapse::SpecificHeatFromDensityTemperature(const Real rho,
                                                         const Real temperature,
                                                         Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoT_(rho, temperature, lambda, lRho, lT, Ye);
  const Real Cv = dEdT_.interpToReal(Ye, lT, lRho);
  return (Cv > robust::EPS() ? Cv : robust::EPS());
}

PORTABLE_INLINE_FUNCTION
Real StellarCollapse::SpecificHeatFromDensityInternalEnergy(const Real rho,
                                                            const Real sie,
                                                            Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoSie_(rho, sie, lambda, lRho, lT, Ye);
  const Real Cv = dEdT_.interpToReal(Ye, lT, lRho);
  return (Cv > robust::EPS() ? Cv : robust::EPS());
}

PORTABLE_INLINE_FUNCTION
Real StellarCollapse::BulkModulusFromDensityTemperature(const Real rho,
                                                        const Real temperature,
                                                        Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoT_(rho, temperature, lambda, lRho, lT, Ye);
  const Real lbmod = lBMod_.interpToReal(Ye, lT, lRho);
  const Real bMod = lB2B_(lbmod);
  return bMod > robust::EPS() ? bMod : robust::EPS();
}

PORTABLE_INLINE_FUNCTION
Real StellarCollapse::GruneisenParamFromDensityTemperature(const Real rho,
                                                           const Real temp,
                                                           Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoT_(rho, temp, lambda, lRho, lT, Ye);
  const Real dpde = dPdE_.interpToReal(Ye, lT, lRho);
  const Real gm1 = std::abs(dpde) / (std::abs(rho) + robust::EPS());
  return gm1;
}

PORTABLE_INLINE_FUNCTION
Real StellarCollapse::BulkModulusFromDensityInternalEnergy(const Real rho, const Real sie,
                                                           Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoSie_(rho, sie, lambda, lRho, lT, Ye);
  const Real lbmod = lBMod_.interpToReal(Ye, lT, lRho);
  const Real bMod = lB2B_(lbmod);
  return bMod;
}

PORTABLE_INLINE_FUNCTION
Real StellarCollapse::GruneisenParamFromDensityInternalEnergy(const Real rho,
                                                              const Real sie,
                                                              Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoSie_(rho, sie, lambda, lRho, lT, Ye);
  const Real dpde = dPdE_.interpToReal(Ye, lT, lRho);
  const Real gm1 = std::abs(dpde) / (std::abs(rho) + robust::EPS());
  return gm1;
}

// TODO(JMM): Fill in this stub if we ever use this EOS in a PTE code.
PORTABLE_INLINE_FUNCTION
void StellarCollapse::DensityEnergyFromPressureTemperature(const Real press,
                                                           const Real temp, Real *lambda,
                                                           Real &rho, Real &sie) const {
  EOS_ERROR("StellarCollapse::DensityEnergyFromPRessureTemperature is a stub");
}

PORTABLE_INLINE_FUNCTION
void StellarCollapse::MassFractionsFromDensityTemperature(
    const Real rho, const Real temperature, Real &Xa, Real &Xh, Real &Xn, Real &Xp,
    Real &Abar, Real &Zbar, Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoT_(rho, temperature, lambda, lRho, lT, Ye);
  Xa = Xa_.interpToReal(Ye, lT, lRho);
  Xh = Xh_.interpToReal(Ye, lT, lRho);
  Xn = Xn_.interpToReal(Ye, lT, lRho);
  Xp = Xp_.interpToReal(Ye, lT, lRho);
  Abar = Abar_.interpToReal(Ye, lT, lRho);
  Zbar = Zbar_.interpToReal(Ye, lT, lRho);
}

PORTABLE_INLINE_FUNCTION
void StellarCollapse::FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv,
                              Real &bmod, const unsigned long output,
                              Real *lambda) const {
  Real lRho, lT, Ye;
  const unsigned long input = ~output;
  if (output == thermalqs::none) {
    UNDEFINED_ERROR;
  }
  if (output & thermalqs::density) {
    EOS_ERROR("StellarCollapse cannot output density at this time");
  }
  if (input & thermalqs::temperature) {
    getLogsFromRhoT_(rho, temp, lambda, lRho, lT, Ye);
  } else if (input & thermalqs::specific_internal_energy) {
    getLogsFromRhoSie_(rho, energy, lambda, lRho, lT, Ye);
  } else {
    UNDEFINED_ERROR;
  }
  if (output & thermalqs::specific_internal_energy) {
    const Real lE = lE_.interpToReal(Ye, lT, lRho);
    energy = le2e_(lE);
  }
  if (output & thermalqs::pressure) {
    const Real lP = lP_.interpToReal(Ye, lT, lRho);
    press = lP2P_(lP);
  }
  if (output & thermalqs::specific_heat) {
    const Real Cv = dEdT_.interpToReal(Ye, lT, lRho);
    cv = (Cv > robust::EPS() ? Cv : robust::EPS());
  }
  if (output & thermalqs::bulk_modulus) {
    const Real lbmod = lBMod_.interpToReal(Ye, lT, lRho);
    bmod = lB2B_(lbmod);
  }
}

PORTABLE_INLINE_FUNCTION
void StellarCollapse::ValuesAtReferenceState(Real &rho, Real &temp, Real &sie,
                                             Real &press, Real &cv, Real &bmod,
                                             Real &dpde, Real &dvdt, Real *lambda) const {
  rho = rhoNormal_;
  temp = TNormal_;
  sie = sieNormal_;
  press = PNormal_;
  cv = CvNormal_;
  bmod = bModNormal_;
  dpde = dPdENormal_;
  dvdt = dVdTNormal_;
  Real lT = lT_(temp);
  lambda[Lambda::Ye] = YeNormal_;
  lambda[Lambda::lT] = lT;
}

inline void StellarCollapse::LoadFromSP5File_(const std::string &filename) {
  herr_t status = H5_SUCCESS;

  hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Offsets
  hid_t metadata = H5Gopen(file, METADATA_NAME, H5P_DEFAULT);
  status += H5LTget_attribute_double(file, METADATA_NAME, SP5::Offsets::sie, &lEOffset_);
  status += H5Gclose(metadata);

  // Databoxes
  status += lP_.loadHDF(file, "logpress");
  status += lE_.loadHDF(file, "logenergy");
  status += dPdRho_.loadHDF(file, "dpdrhoe");
  status += dPdE_.loadHDF(file, "dpderho");
  status += dEdT_.loadHDF(file, "dedt");
  status += entropy_.loadHDF(file, "entropy");
  status += Xa_.loadHDF(file, "Xa");
  status += Xh_.loadHDF(file, "Xh");
  status += Xn_.loadHDF(file, "Xn");
  status += Xp_.loadHDF(file, "Xp");
  status += Abar_.loadHDF(file, "Abar");
  status += Zbar_.loadHDF(file, "Zbar");
  status += lBMod_.loadHDF(file, "logbulkmodulus");
  status += eCold_.loadHDF(file, "ecold");
  status += eHot_.loadHDF(file, "ehot");

  status += H5Fclose(file);
  if (status != H5_SUCCESS) {
    EOS_ERROR("[StellarCollapse::Load]: There was a problem with HDF5\n");
  }

  // bounds, etc.
  auto YeGrid = lP_.range(2);
  auto lTGrid = lP_.range(1);
  auto lRGrid = lP_.range(0);
  numRho_ = lRGrid.nPoints();
  numT_ = lTGrid.nPoints();
  numYe_ = YeGrid.nPoints();
  lRhoMin_ = lRGrid.min();
  lRhoMax_ = lRGrid.max();
  lTMin_ = lTGrid.min();
  lTMax_ = lTGrid.max();
  YeMin_ = YeGrid.min();
  YeMax_ = YeGrid.max();
  sieMin_ = eCold_.min();
  sieMax_ = eHot_.max();
}

// Read data directly from a stellar collapse eos file
inline void StellarCollapse::LoadFromStellarCollapseFile_(const std::string &filename,
                                                          bool filter_bmod) {
  // Start HDF5 stuff
  // ---------------------------------------------------------------------
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  herr_t status = H5_SUCCESS;

  numRho_ = readSCInt_(file_id, "pointsrho");
  numT_ = readSCInt_(file_id, "pointstemp");
  numYe_ = readSCInt_(file_id, "pointsye");

  status = H5LTread_dataset_double(file_id, "energy_shift", &lEOffset_);
  if (status < 0) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("An HDF5 error ocurred while reading energy_shift\n");
  }

  Real lRhoMin, lRhoMax, lTMin, lTMax;
  readBounds_(file_id, "logrho", numRho_, lRhoMin, lRhoMax);
  readBounds_(file_id, "logtemp", numT_, lTMin, lTMax);
  readBounds_(file_id, "ye", numYe_, YeMin_, YeMax_);

  // Convert temperature to MeV
  const Real logMev2K = std::log10(MeV2K_);
  lTMin += logMev2K;
  lTMax += logMev2K;

  // Generate grids for reading stellar collapse format tables
  Grid_t Ye_grid(YeMin_, YeMax_, numYe_);
  Grid_t lT_grid(lTMin, lTMax, numT_);
  Grid_t lRho_grid(lRhoMin, lRhoMax, numRho_);

  // Tables
  // fastest -> slowest:
  // Ye -> lT -> lRho
  // TODO(JMM): Might be worth re-ordering these.
  // For cache efficiency, lT should be fastest, not lRho.

  // Dependent vars in log space, so interpolation is log-log
  readSCDset_(file_id, "logpress", Ye_grid, lT_grid, lRho_grid, lP_);
  readSCDset_(file_id, "logenergy", Ye_grid, lT_grid, lRho_grid, lE_);

  // Dependent vars not in log space, so interpolation is log-linear
  readSCDset_(file_id, "dpdrhoe", Ye_grid, lT_grid, lRho_grid, dPdRho_);
  readSCDset_(file_id, "dpderho", Ye_grid, lT_grid, lRho_grid, dPdE_);
  readSCDset_(file_id, "dedt", Ye_grid, lT_grid, lRho_grid, dEdT_);

  // TODO(JMM): entropy, mass fractions, and average atomic mass and
  // numbers aren't exposed in eos_variant. So you need to pull the
  // type out.
  readSCDset_(file_id, "entropy", Ye_grid, lT_grid, lRho_grid, entropy_);
  readSCDset_(file_id, "Xa", Ye_grid, lT_grid, lRho_grid, Xa_);
  readSCDset_(file_id, "Xh", Ye_grid, lT_grid, lRho_grid, Xh_);
  readSCDset_(file_id, "Xn", Ye_grid, lT_grid, lRho_grid, Xn_);
  readSCDset_(file_id, "Xp", Ye_grid, lT_grid, lRho_grid, Xp_);
  readSCDset_(file_id, "Abar", Ye_grid, lT_grid, lRho_grid, Abar_);
  readSCDset_(file_id, "Zbar", Ye_grid, lT_grid, lRho_grid, Zbar_);

  H5Fclose(file_id);
  // -----------------------------------------------------------------------
  // End HDF5 stuff

  // TODO(BRR): in the future, if reading in chemical potentials, convert from MeV to K

  // Convert MeV to K in tabulated data
  for (int iY = 0; iY < numYe_; ++iY) {
    for (int iT = 0; iT < numT_; ++iT) {
      for (int irho = 0; irho < numRho_; ++irho) {
        dEdT_(iY, iT, irho) *= K2MeV_;
      }
    }
  }

  // Filter thermo derivs before reinterpolating
  if (filter_bmod) {
    medianFilter_(dPdRho_);
    medianFilter_(dPdE_);
    medianFilter_(dEdT_);
  }

  // Re-interpolate tables in case we want fast-log gridding
  DataBox scratch(numYe_, numT_, numRho_);
  // logged quantities
  dataBoxToFastLogs(lP_, scratch, true);
  dataBoxToFastLogs(lE_, scratch, true);
  // linear quantities
  dataBoxToFastLogs(dPdRho_, scratch, false);
  dataBoxToFastLogs(dPdE_, scratch, false);
  dataBoxToFastLogs(dEdT_, scratch, false);
  // non-standard quantities
  dataBoxToFastLogs(entropy_, scratch, false);
  dataBoxToFastLogs(Xa_, scratch, false);
  dataBoxToFastLogs(Xh_, scratch, false);
  dataBoxToFastLogs(Xn_, scratch, false);
  dataBoxToFastLogs(Xp_, scratch, false);
  dataBoxToFastLogs(Abar_, scratch, false);
  dataBoxToFastLogs(Zbar_, scratch, false);

  // Generate bounds
  Ye_grid = lP_.range(2);
  lT_grid = lP_.range(1);
  lRho_grid = lP_.range(0);
  YeMin_ = Ye_grid.min();
  YeMax_ = Ye_grid.max();
  lTMin_ = lT_grid.min();
  lTMax_ = lT_grid.max();
  lRhoMin_ = lRho_grid.min();
  lRhoMax_ = lRho_grid.max();

  // Finally, compute bulk modulus, hot and cold curves from
  // re-interpolated data.
  computeBulkModulus_();
  computeColdAndHotCurves_();
}

inline int StellarCollapse::readSCInt_(const hid_t &file_id, const std::string &name) {
  int data;
  hsize_t one = 1;
  hsize_t file_grid_dims[] = {1};
  hsize_t file_start[] = {0};
  hsize_t file_count[] = {1};
  hid_t filespace = H5Screate_simple(1, file_grid_dims, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  hid_t memspace = H5Screate_simple(1, &one, NULL);

  hid_t dset_id = H5Dopen(file_id, name.c_str(), H5P_DEFAULT);
  herr_t status =
      H5Dread(dset_id, H5T_STD_I32LE, memspace, filespace, H5P_DEFAULT, &data);
  H5Dclose(dset_id);
  H5Sclose(memspace);
  H5Sclose(filespace);
  if (status < 0) {
    std::string msg = "Failed to read dataset " + name + "\n";
    EOS_ERROR(msg.c_str());
  }
  return data;
}

// Read (inclusive) bounds on the independent variables Assumes
// uniform tables in log-space. The file requires we read in the
// array, but we only need the bounds themselves.
inline void StellarCollapse::readBounds_(const hid_t &file_id, const std::string &name,
                                         int size, Real &lo, Real &hi) {
  std::vector<Real> table(size);
  herr_t status = H5LTread_dataset_double(file_id, name.c_str(), table.data());
  if (status != H5_SUCCESS) {
    EOS_ERROR("An HDF5 error ocurred while reading bounds");
  }
  lo = table[0];
  hi = table[size - 1];
}

/* For some insane reason, HDF5 won't read into the pointer exposed by
 * `DataBox.data()` So we hack our way around this by using a std::vector and
 * then copying the elements. There's some evidence this is an HDF5 bug. See:
 * https://forum.hdfgroup.org/t/is-this-a-bug-in-hdf5-1-8-6/2211
 */
inline void StellarCollapse::readSCDset_(const hid_t &file_id, const std::string &name,
                                         const Grid_t &Ye_grid, const Grid_t &lT_grid,
                                         const Grid_t &lRho_grid, DataBox &db) {
  herr_t exists = H5LTfind_dataset(file_id, name.c_str());
  if (!exists) {
    std::string msg = "Tried to read dataset " + name + " but it doesn't exist\n";
    EOS_ERROR(msg.c_str());
  }
  db.resize(Ye_grid.nPoints(), lT_grid.nPoints(), lRho_grid.nPoints());
  herr_t status = H5LTread_dataset(file_id, name.c_str(), H5T_REAL, db.data());
  if (status != H5_SUCCESS) {
    std::string msg = "Tried to read dataset " + name + " but it failed exist\n";
    PORTABLE_ALWAYS_THROW_OR_ABORT(msg.c_str());
  }

  // bounds
  db.setRange(2, Ye_grid);
  db.setRange(1, lT_grid);
  db.setRange(0, lRho_grid);
}

// Reinterpolate tab from its original grid spacing to the one using
// the native log gridding for stellar collapse (usually fast logs).
// Scratch is used as a temporary storage buffer and is assumed to be
// the same shape as db
// Assume index 3 is linear, indexes 2 and 1 are logarithmic
inline void StellarCollapse::dataBoxToFastLogs(DataBox &db, DataBox &scratch,
                                               bool dependent_var_log) {
  auto log10toNQT = [](const Real x) { return FastMath::log10(std::pow(10, x)); };
  auto NQTtolog10 = [](const Real x) { return std::log10(FastMath::pow10(x)); };
  auto gridToNQT = [&](const Grid_t &g) {
    const Real l10min = g.min();
    const Real l10max = g.max();
    const Real lmin = log10toNQT(l10min);
    const Real lmax = log10toNQT(l10max);
    return Grid_t(lmin, lmax, g.nPoints());
  };

  auto &r2 = db.range(2);
  auto &r1 = db.range(1);
  auto &r0 = db.range(0);

  Grid_t newr1 = gridToNQT(r1);
  Grid_t newr0 = gridToNQT(r0);

  for (int i2 = 0; i2 < r2.nPoints(); ++i2) {
    Real x2 = r2.x(i2);
    for (int i1 = 0; i1 < newr1.nPoints(); ++i1) {
      Real lx1 = newr1.x(i1);
      Real l10x1 = NQTtolog10(lx1);
      for (int i0 = 0; i0 < newr0.nPoints(); ++i0) {
        Real lx0 = newr0.x(i0);
        Real l10x0 = NQTtolog10(lx0);
        Real val = db.interpToReal(x2, l10x1, l10x0);
        if (dependent_var_log) {
          val = log10toNQT(val);
        }
        scratch(i2, i1, i0) = val;
      }
    }
  }
  for (int i = 0; i < db.size(); ++i) {
    db(i) = scratch(i);
  }
  // range(2) is already ok
  db.setRange(1, newr1);
  db.setRange(0, newr0);
}

inline void StellarCollapse::medianFilter_(DataBox &db) {
  DataBox tmp;
  tmp.copy(db);
  medianFilter_(tmp, db);
  free(tmp);
}

inline void StellarCollapse::medianFilter_(const DataBox &in, DataBox &out) {
  Real buffer[MF_S];
  // filter, overwriting as needed
  for (int iY = MF_W; iY < numYe_ - MF_W; ++iY) {
    for (int iT = MF_W; iT < numT_ - MF_W; ++iT) {
      for (int irho = MF_W; irho < numRho_ - MF_W; ++irho) {
        out(iY, iT, irho) = in(iY, iT, irho);
        fillMedianBuffer_(buffer, MF_W, iY, iT, irho, in);
        Real point = in(iY, iT, irho);
        Real avg = findMedian_(buffer, MF_S);
        int bad = std::abs(avg - point) / std::abs(avg) > DELTASMOOTH;
        if (bad) out(iY, iT, irho) = avg;
      }
    }
  }
}

inline void StellarCollapse::fillMedianBuffer_(Real buffer[], int width, int iY, int iT,
                                               int irho, const DataBox &tab) const {
  int i = 0;
  for (int iWy = -width; iWy <= width; iWy++) {
    for (int iWt = -width; iWt <= width; iWt++) {
      for (int iWr = -width; iWr <= width; iWr++) {
        buffer[i++] = tab(iY + iWy, iT + iWt, irho + iWr);
      }
    }
  }
}

// Note modifies buffer
inline Real StellarCollapse::findMedian_(Real buffer[], int size) const {
  std::qsort(buffer, size, sizeof(Real), [](const void *ap, const void *bp) {
    double a = *((double *)ap);
    double b = *((double *)bp);
    return (a > b) - (a < b);
  });
  if (size % 2 == 0) {
    return 0.5 * (buffer[(size - 1) / 2] + buffer[size / 2]);
  }
  return buffer[size / 2];
}

// TODO(JMM): I tabulate bulk modulus on a log scale. Should I?
inline void StellarCollapse::computeBulkModulus_() {
  lBMod_.copyMetadata(lP_);
  for (int iY = 0; iY < numYe_; ++iY) {
    Real Ye = lBMod_.range(2).x(iY);
    for (int iT = 0; iT < numT_; ++iT) {
      Real lT = lBMod_.range(1).x(iT);
      for (int irho = 0; irho < numRho_; ++irho) {
        Real lRho = lBMod_.range(0).x(irho);
        Real rho = rho_(lRho);
        Real lP = lP_(iY, iT, irho);
        Real P = lP2P_(lP);
        Real PoR = robust::ratio(P, rho);
        // assume table is hardened
        Real bMod = rho * dPdRho_(iY, iT, irho) + PoR * dPdE_(iY, iT, irho);
        if (bMod < robust::EPS()) bMod = robust::EPS();
        lBMod_(iY, iT, irho) = B2lB_(bMod);
      }
    }
  }
}

inline void StellarCollapse::computeColdAndHotCurves_() {
  eCold_.resize(numYe_, numRho_);
  eHot_.resize(numYe_, numRho_);
  int iTCold = 0;
  int iTHot = numT_ - 1;
  for (int iY = 0; iY < numYe_; ++iY) {
    for (int irho = 0; irho < numRho_; ++irho) {
      Real lECold = lE_(iY, iTCold, irho);
      Real lEHot = lE_(iY, iTHot, irho);
      eCold_(iY, irho) = le2e_(lECold);
      eHot_(iY, irho) = le2e_(lEHot);
    }
  }
  sieMin_ = eCold_.min();
  sieMax_ = eHot_.max();
  eCold_.setRange(0, lRhoMin_, lRhoMax_, numRho_);
  eCold_.setRange(1, YeMin_, YeMax_, numYe_);
  eHot_.setRange(0, lRhoMin_, lRhoMax_, numRho_);
  eHot_.setRange(1, YeMin_, YeMax_, numYe_);
}

inline void StellarCollapse::setNormalValues_() {
  const Real lT = lT_(TNormal_);
  const Real lRho = lRho_(rhoNormal_);
  const Real Ye = YeNormal_;

  const Real lE = lE_.interpToReal(Ye, lT, lRho);
  sieNormal_ = le2e_(lE);

  const Real lP = lP_.interpToReal(Ye, lT, lRho);
  PNormal_ = lP2P_(lP);

  const Real entropy = entropy_.interpToReal(Ye, lT, lRho);
  SNormal_ = entropy;

  const Real Cv = dEdT_.interpToReal(Ye, lT, lRho);
  CvNormal_ = (Cv > robust::EPS() ? Cv : robust::EPS());

  const Real lB = lBMod_.interpToReal(Ye, lT, lRho);
  const Real bMod = lB2B_(lB);
  bModNormal_ = bMod > robust::EPS() ? bMod : robust::EPS();

  dPdENormal_ = dPdE_.interpToReal(Ye, lT, lRho);

  Real dPdR = dPdRho_.interpToReal(Ye, lT, lRho);
  dVdTNormal_ = dPdENormal_ * CvNormal_ / (rhoNormal_ * rhoNormal_ * dPdR);
}

PORTABLE_INLINE_FUNCTION
Real StellarCollapse::lTFromlRhoSie_(const Real lRho, const Real sie,
                                     Real *lambda) const noexcept {
  checkLambda_(lambda);
  RootFinding1D::Status status = RootFinding1D::Status::SUCCESS;
  using RootFinding1D::regula_falsi;
  Real lT;
  Real Ye = lambda[Lambda::Ye];
  Real lTGuess = lambda[Lambda::lT];

  const RootFinding1D::RootCounts *pcounts =
      (memoryStatus_ == DataStatus::OnDevice) ? nullptr : &counts;

  // If sie above hot curve or below cold curve, force it onto the table.
  // TODO(JMM): Rethink this as needed.
  if (sie <= eCold_.interpToReal(Ye, lRho)) {
    lT = lTGuess = lTMin_;
    if (pcounts != nullptr) {
      pcounts->increment(0);
    }
  } else if (sie >= eHot_.interpToReal(Ye, lRho)) {
    lT = lTGuess = lTMax_;
    if (pcounts != nullptr) {
      pcounts->increment(0);
    }
  } else {
    // if the guess isn't in the bounds, bound it
    if (!(lTMin_ <= lTGuess && lTGuess <= lTMax_)) {
      lTGuess = 0.5 * (lTMin_ + lTMax_);
    }
    // Get log(sie)
    Real lE = e2le_(sie);
    const callable_interp::LogT lEFunc(lE_, Ye, lRho);
    status = regula_falsi(lEFunc, lE, lTGuess, lTMin_, lTMax_, ROOT_THRESH, ROOT_THRESH,
                          lT, pcounts);
    if (status != RootFinding1D::Status::SUCCESS) {
#if STELLAR_COLLAPSE_EOS_VERBOSE
      std::stringstream errorMessage;
      errorMessage << "Inverting log(sie) table for log(T) failed\n"
                   << "Ye      = " << Ye << "\n"
                   << "lRho    = " << lRho << "\n"
                   << "sie     = " << sie << "\n"
                   << "lE      = " << lE << "\n"
                   << "lTGuess = " << lTGuess << std::endl;
      EOS_ERROR(errorMessage.str().c_str());
#endif // STELLAR_COLLAPSE_EOS_VERBOSE
      lT = lTGuess;
    }
  }
  if (memoryStatus_ != DataStatus::OnDevice) {
    status_ = status;
  }
  lambda[Lambda::lT] = lT;
  return lT;
}
} // namespace singularity

#endif // SINGULARITY_USE_SPINER_WITH_HDF5
#endif // _SINGULARITY_EOS_EOS_EOS_STELLAR_COLLAPSE_HPP_
