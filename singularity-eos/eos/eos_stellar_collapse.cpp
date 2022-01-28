//------------------------------------------------------------------------------
// © 2021. Triad National Security, LLC. All rights reserved.  This
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

#ifdef SPINER_USE_HDF

// C++ includes
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// C includes
#include <cstdlib>
#include <hdf5.h>
#include <hdf5_hl.h>

// singularity includes
#include <singularity-eos/eos/eos.hpp>

#include <ports-of-call/portability.hpp>
#include <root-finding-1d/root_finding.hpp>
#include <sp5/singularity_eos_sp5.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/sp5.hpp>

#define STELLAR_COLLAPSE_EOS_VERBOSE (0)

namespace callable_interp {

class LogT {
 public:
  PORTABLE_INLINE_FUNCTION
  LogT(const Spiner::DataBox &field, const Real Ye, const Real lRho)
      : field_(field), Ye_(Ye), lRho_(lRho) {}
  PORTABLE_INLINE_FUNCTION Real operator()(const Real lT) const {
    return field_.interpToReal(Ye_, lT, lRho_);
  }

 private:
  const Spiner::DataBox &field_;
  const Real Ye_, lRho_;
};

} // namespace callable_interp

namespace singularity {

// For some reason, the linker doesn't like this being a member field
// of StellarCollapse.  So we'll make it a global variable.
constexpr char METADATA_NAME[] = "Metadata";

StellarCollapse::StellarCollapse(const std::string &filename, bool use_sp5,
                                 bool filter_bmod)
    : filename_(filename.c_str()) {

  if (use_sp5) {
    LoadFromSP5File_(filename);
  } else {
    LoadFromStellarCollapseFile_(filename);
    if (filter_bmod) {
      medianFilter_(dPdRho_); // needed if pulling the data
      medianFilter_(dPdE_);   // directly from Stellar Collapse tables
      medianFilter_(dEdT_);
    }
    computeBulkModulus_();
    computeColdAndHotCurves_();
  }
  setNormalValues_();
}

// Saves to an SP5 file
void StellarCollapse::Save(const std::string &filename) {
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

StellarCollapse StellarCollapse::GetOnDevice() {
  StellarCollapse other;
  other.lP_ = Spiner::getOnDeviceDataBox(lP_);
  other.lE_ = Spiner::getOnDeviceDataBox(lE_);
  other.dPdRho_ = Spiner::getOnDeviceDataBox(dPdRho_);
  other.dPdE_ = Spiner::getOnDeviceDataBox(dPdE_);
  other.dEdT_ = Spiner::getOnDeviceDataBox(dEdT_);
  other.entropy_ = Spiner::getOnDeviceDataBox(entropy_);
  other.Xa_ = Spiner::getOnDeviceDataBox(Xa_);
  other.Xh_ = Spiner::getOnDeviceDataBox(Xh_);
  other.Xn_ = Spiner::getOnDeviceDataBox(Xn_);
  other.Xp_ = Spiner::getOnDeviceDataBox(Xp_);
  other.Abar_ = Spiner::getOnDeviceDataBox(Abar_);
  other.Zbar_ = Spiner::getOnDeviceDataBox(Zbar_);
  other.lBMod_ = Spiner::getOnDeviceDataBox(lBMod_);
  other.eCold_ = Spiner::getOnDeviceDataBox(eCold_);
  other.eHot_ = Spiner::getOnDeviceDataBox(eHot_);
  other.memoryStatus_ = DataStatus::OnDevice;
  other.numRho_ = numRho_;
  other.numT_ = numT_;
  other.numYe_ = numYe_;
  other.lTMin_ = lTMin_;
  other.lTMax_ = lTMax_;
  other.YeMin_ = YeMin_;
  other.YeMax_ = YeMax_;
  other.lEOffset_ = lEOffset_;
  other.sieNormal_ = sieNormal_;
  other.PNormal_ = PNormal_;
  other.CvNormal_ = CvNormal_;
  other.bModNormal_ = bModNormal_;
  other.dPdENormal_ = dPdENormal_;
  other.dVdTNormal_ = dVdTNormal_;
  other.status_ = status_;
  return other;
}

void StellarCollapse::Finalize() {
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

PORTABLE_FUNCTION
Real StellarCollapse::TemperatureFromDensityInternalEnergy(const Real rho, const Real sie,
                                                           Real *lambda) const {
  const Real lRho = lRho_(rho);
  const Real lT = lTFromlRhoSie_(lRho, sie, lambda);
  return T_(lT);
}

PORTABLE_FUNCTION
Real StellarCollapse::InternalEnergyFromDensityTemperature(const Real rho,
                                                           const Real temp,
                                                           Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoT_(rho, temp, lambda, lRho, lT, Ye);
  const Real lE = lE_.interpToReal(Ye, lT, lRho);
  return le2e_(lE);
}

PORTABLE_FUNCTION
Real StellarCollapse::PressureFromDensityTemperature(const Real rho,
                                                     const Real temperature,
                                                     Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoT_(rho, temperature, lambda, lRho, lT, Ye);
  const Real lP = lP_.interpToReal(Ye, lT, lRho);
  return lP2P_(lP);
}

PORTABLE_FUNCTION
Real StellarCollapse::PressureFromDensityInternalEnergy(const Real rho, const Real sie,
                                                        Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoSie_(rho, sie, lambda, lRho, lT, Ye);
  const Real lP = lP_.interpToReal(Ye, lT, lRho);
  return lP2P_(lP);
}

PORTABLE_FUNCTION
Real StellarCollapse::SpecificHeatFromDensityTemperature(const Real rho,
                                                         const Real temperature,
                                                         Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoT_(rho, temperature, lambda, lRho, lT, Ye);
  const Real Cv = dEdT_.interpToReal(Ye, lT, lRho);
  return (Cv > EPS ? Cv : EPS);
}

PORTABLE_FUNCTION
Real StellarCollapse::SpecificHeatFromDensityInternalEnergy(const Real rho,
                                                            const Real sie,
                                                            Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoSie_(rho, sie, lambda, lRho, lT, Ye);
  const Real Cv = dEdT_.interpToReal(Ye, lT, lRho);
  return (Cv > EPS ? Cv : EPS);
}

PORTABLE_FUNCTION
Real StellarCollapse::BulkModulusFromDensityTemperature(const Real rho,
                                                        const Real temperature,
                                                        Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoT_(rho, temperature, lambda, lRho, lT, Ye);
  const Real lbmod = lBMod_.interpToReal(Ye, lT, lRho);
  const Real bMod = lB2B_(lbmod);
  return bMod > EPS ? bMod : EPS;
}

PORTABLE_FUNCTION
Real StellarCollapse::GruneisenParamFromDensityTemperature(const Real rho,
                                                           const Real temp,
                                                           Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoT_(rho, temp, lambda, lRho, lT, Ye);
  const Real dpde = dPdE_.interpToReal(Ye, lT, lRho);
  const Real gm1 = std::abs(dpde) / (std::abs(rho) + EPS);
  return gm1;
}

PORTABLE_FUNCTION
Real StellarCollapse::BulkModulusFromDensityInternalEnergy(const Real rho, const Real sie,
                                                           Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoSie_(rho, sie, lambda, lRho, lT, Ye);
  const Real lbmod = lBMod_.interpToReal(Ye, lT, lRho);
  const Real bMod = lB2B_(lbmod);
  return bMod;
}

PORTABLE_FUNCTION
Real StellarCollapse::GruneisenParamFromDensityInternalEnergy(const Real rho,
                                                              const Real sie,
                                                              Real *lambda) const {
  Real lRho, lT, Ye;
  getLogsFromRhoSie_(rho, sie, lambda, lRho, lT, Ye);
  const Real dpde = dPdE_.interpToReal(Ye, lT, lRho);
  const Real gm1 = std::abs(dpde) / (std::abs(rho) + EPS);
  return gm1;
}

// TODO(JMM): Fill in this stub if we ever use this EOS in a PTE code.
PORTABLE_FUNCTION
void StellarCollapse::DensityEnergyFromPressureTemperature(const Real press,
                                                           const Real temp, Real *lambda,
                                                           Real &rho, Real &sie) const {
  EOS_ERROR("StellarCollapse::DensityEnergyFromPRessureTemperature is a stub");
}

PORTABLE_FUNCTION
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
    cv = (Cv > EPS ? Cv : EPS);
  }
  if (output & thermalqs::bulk_modulus) {
    const Real lbmod = lBMod_.interpToReal(Ye, lT, lRho);
    bmod = lB2B_(lbmod);
  }
}

PORTABLE_FUNCTION
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

void StellarCollapse::LoadFromSP5File_(const std::string &filename) {
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
void StellarCollapse::LoadFromStellarCollapseFile_(const std::string &filename) {
  // Open the file.
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  herr_t status = H5_SUCCESS;

  numRho_ = readSCInt_(file_id, "pointsrho");
  numT_ = readSCInt_(file_id, "pointstemp");
  numYe_ = readSCInt_(file_id, "pointsye");

  status = H5LTread_dataset_double(file_id, "energy_shift", &lEOffset_);
  if (status < 0) {
    EOS_ERROR("An HDF5 error ocurred while reading energy_shift\n");
  }

  // Bounds
  readBounds_(file_id, "logrho", numRho_, lRhoMin_, lRhoMax_);
  readBounds_(file_id, "logtemp", numT_, lTMin_, lTMax_);
  const Real logMev2K = std::log10(MeV2K_);
  lTMin_ += logMev2K;
  lTMax_ += logMev2K;
  readBounds_(file_id, "ye", numYe_, YeMin_, YeMax_);

  // Tables
  // fastest -> slowest:
  // Ye -> lT -> lRho
  // TODO(JMM): Might be worth re-ordering these.
  // For cache efficiency, lT should be fastest, not lRho.
  readSCDset_(file_id, "logpress", lP_);
  readSCDset_(file_id, "logenergy", lE_);
  readSCDset_(file_id, "dpdrhoe", dPdRho_);
  readSCDset_(file_id, "dpderho", dPdE_);
  readSCDset_(file_id, "dedt", dEdT_);

  // TODO(JMM): entropy, mass fractions, and average atomic mass and
  // numbers aren't exposed in eos_variant. So you need to pull the
  // type out.
  readSCDset_(file_id, "entropy", entropy_);
  readSCDset_(file_id, "Xa", Xa_);
  readSCDset_(file_id, "Xh", Xh_);
  readSCDset_(file_id, "Xn", Xn_);
  readSCDset_(file_id, "Xp", Xp_);
  readSCDset_(file_id, "Abar", Abar_);
  readSCDset_(file_id, "Zbar", Zbar_);

  // TODO(BRR): in the future, if reading in chemical potentials, convert from MeV to K

  // Convert MeV to K in tabulated data
  for (int iY = 0; iY < numYe_; ++iY) {
    for (int iT = 0; iT < numT_; ++iT) {
      for (int irho = 0; irho < numRho_; ++irho) {
        dEdT_(iY, iT, irho) *= K2MeV_;
      }
    }
  }

  H5Fclose(file_id);
}

int StellarCollapse::readSCInt_(const hid_t &file_id, const std::string &name) {
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
void StellarCollapse::readBounds_(const hid_t &file_id, const std::string &name, int size,
                                  Real &lo, Real &hi) {
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
void StellarCollapse::readSCDset_(const hid_t &file_id, const std::string &name,
                                  Spiner::DataBox &db) {
  herr_t exists = H5LTfind_dataset(file_id, name.c_str());
  if (!exists) {
    std::string msg = "Tried to read dataset " + name + " but it doesn't exist\n";
    EOS_ERROR(msg.c_str());
  }
  db.resize(numYe_, numT_, numRho_);
  herr_t status = H5LTread_dataset(file_id, name.c_str(), H5T_REAL, db.data());
  if (status != H5_SUCCESS) {
    std::string msg = "Tried to read dataset " + name + " but it failed exist\n";
    EOS_ERROR(msg.c_str());
  }

  // bounds
  db.setRange(2, YeMin_, YeMax_, numYe_);
  db.setRange(1, lTMin_, lTMax_, numT_);
  db.setRange(0, lRhoMin_, lRhoMax_, numRho_);
}

void StellarCollapse::medianFilter_(Spiner::DataBox &db) {
  Spiner::DataBox tmp;
  tmp.copy(db);
  medianFilter_(tmp, db);
}

void StellarCollapse::medianFilter_(const Spiner::DataBox &in, Spiner::DataBox &out) {
  Real buffer[MF_S];
  // filter, overwriting as needed
  for (int iY = MF_W; iY < numYe_ - MF_W; ++iY) {
    for (int iT = MF_W; iT < numT_ - MF_W; ++iT) {
      for (int irho = MF_W; irho < numRho_ - MF_W; ++irho) {
        out(iY, iT, irho) = in(iY, iT, irho);
        fillMedianBuffer_(buffer, MF_W, iY, iT, irho, in);
        Real point = in(iY, iT, irho);
        Real avg = findMedian_(buffer, MF_S);
        int bad = std::abs(avg - point) / std::abs(avg) > EPSSMOOTH;
        if (bad) out(iY, iT, irho) = avg;
      }
    }
  }
}

void StellarCollapse::fillMedianBuffer_(Real buffer[], int width, int iY, int iT,
                                        int irho, const Spiner::DataBox &tab) const {
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
Real StellarCollapse::findMedian_(Real buffer[], int size) const {
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
void StellarCollapse::computeBulkModulus_() {
  lBMod_.copyMetadata(lP_);
  for (int iY = 0; iY < numYe_; ++iY) {
    Real Ye = lBMod_.range(2).x(iY);
    for (int iT = 0; iT < numT_; ++iT) {
      Real lT = lBMod_.range(1).x(iT);
      for (int irho = 0; irho < numRho_; ++irho) {
        Real lRho = lBMod_.range(0).x(irho);
        Real rho = rho_(lRho);
        Real lP = lP_(iY, iT, irho);
        Real lPoR = lP - lRho;
        Real PoR = fromLog_(lPoR, 0.0);
        // assume table is hardened
        Real bMod = rho * dPdRho_(iY, iT, irho) + PoR * dPdE_(iY, iT, irho);
        if (bMod < EPS) bMod = EPS;
        lBMod_(iY, iT, irho) = B2lB_(bMod);
      }
    }
  }
}

void StellarCollapse::computeColdAndHotCurves_() {
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

void StellarCollapse::setNormalValues_() {
  const Real lT = lT_(TNormal_);
  const Real lRho = lRho_(rhoNormal_);
  const Real Ye = YeNormal_;

  const Real lE = lE_.interpToReal(Ye, lT, lRho);
  sieNormal_ = le2e_(lE);

  const Real lP = lP_.interpToReal(Ye, lT, lRho);
  PNormal_ = lP2P_(lP);

  const Real Cv = dEdT_.interpToReal(Ye, lT, lRho);
  CvNormal_ = (Cv > EPS ? Cv : EPS);

  const Real lB = lBMod_.interpToReal(Ye, lT, lRho);
  const Real bMod = lB2B_(lB);
  bModNormal_ = bMod > EPS ? bMod : EPS;

  dPdENormal_ = dPdE_.interpToReal(Ye, lT, lRho);

  Real dPdR = dPdRho_.interpToReal(Ye, lT, lRho);
  dVdTNormal_ = dPdENormal_ * CvNormal_ / (rhoNormal_ * rhoNormal_ * dPdR);
}

PORTABLE_FUNCTION
Real StellarCollapse::lTFromlRhoSie_(const Real lRho, const Real sie,
                                     Real *lambda) const noexcept {
  checkLambda_(lambda);
  RootFinding1D::Status status = RootFinding1D::Status::SUCCESS;
  using RootFinding1D::findRoot;
  Real lT;
  Real Ye = lambda[Lambda::Ye];
  Real lTGuess = lambda[Lambda::lT];

  // If sie above hot curve or below cold curve, force it onto the table.
  // TODO(JMM): Rethink this as needed.
  if (sie <= eCold_.interpToReal(Ye, lRho)) {
    lT = lTGuess = lTMin_;
    counts.increment(0);
  } else if (sie >= eHot_.interpToReal(Ye, lRho)) {
    lT = lTGuess = lTMax_;
    counts.increment(0);
  } else {
    // if the guess isn't in the bounds, bound it
    if (!(lTMin_ <= lTGuess && lTGuess <= lTMax_)) {
      lTGuess = 0.5 * (lTMin_ + lTMax_);
    }
    // Get log(sie)
    Real lE = e2le_(sie);
    const callable_interp::LogT lEFunc(lE_, Ye, lRho);
    status = findRoot(lEFunc, lE, lTGuess, lTMin_, lTMax_, ROOT_THRESH, ROOT_THRESH, lT,
                      counts);
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
  status_ = status;
  lambda[Lambda::lT] = lT;
  return lT;
}
} // namespace singularity

#endif // SPINER_USE_HDF
