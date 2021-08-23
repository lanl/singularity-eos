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

#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
// #include <iostream> // debug
// #include <stdio.h> // debug

#include <hdf5.h>
#include <hdf5_hl.h>

#include <eos/eos.hpp>
#include <sp5/singularity_eos_sp5.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/sp5.hpp>
#include <sp5/singularity_eos_sp5.hpp>
#include <root-finding-1d/root_finding.hpp>
#include <ports-of-call/portability.hpp>

#define SPINER_EOS_VERBOSE (0)

namespace singularity {

// TODO: we're using log-linear interpolation, not log-log
// this may be suboptimal. We may want a way to do some variables
// in log-log. In particular, it might be good to do:
// pressure, energy, and bulk modulus in log-log.
// ~JMM

// replace lambdas with callable
namespace callable_interp {
  
  class l_interp {
  private:
    const Spiner::DataBox& field;
    const Real fixed;
  public:
    PORTABLE_INLINE_FUNCTION
    l_interp(const Spiner::DataBox& field_, const Real fixed_)
      : field{field_},
        fixed{fixed_}
    { }
    
    PORTABLE_INLINE_FUNCTION Real operator()(const Real x) const {
      return field.interpToReal(x, fixed);
    }
  };
  
  class r_interp {
  private:
    const Spiner::DataBox& field;
    const Real fixed;
  public:
    PORTABLE_INLINE_FUNCTION
    r_interp(const Spiner::DataBox& field_, const Real fixed_)
      : field{field_},
        fixed{fixed_}
    { }
    
    PORTABLE_INLINE_FUNCTION Real operator()(const Real x) const {
      return field.interpToReal(fixed, x);
    }
  };

  class prod_interp_1d {
  private:
    const Spiner::DataBox& field1, field2;
    const Real r;
  public:
    PORTABLE_INLINE_FUNCTION
    prod_interp_1d(const Spiner::DataBox& field1_,
                   const Spiner::DataBox& field2_,
                   const Real r_)
      : field1{field1_}, field2{field2_}, r{r_}
    { }
    PORTABLE_INLINE_FUNCTION Real operator()(const Real x) const {
      return field1.interpToReal(x)*field2.interpToReal(x)*r*x;
    }
  };

  class interp {
  private:
    const Spiner::DataBox& field;
  public:
    PORTABLE_INLINE_FUNCTION
    interp(const Spiner::DataBox& field_) : field(field_) { }
    PORTABLE_INLINE_FUNCTION Real operator()(const Real x) const {
      return field.interpToReal(x);
    }
  };
} // namespace callable_interp

SpinerEOSDependsRhoT::SpinerEOSDependsRhoT(const std::string& filename,
                                           int matid,
                                           bool reproducibility_mode)
  : filename_(filename.c_str())
  , matid_(matid)
  , reproducible_(reproducibility_mode)
  , status_(RootFinding1D::Status::SUCCESS)
  , memoryStatus_(DataStatus::OnHost) {

  std::string matid_str = std::to_string(matid);
  hsize_t dims;
  size_t size;
  H5T_class_t tclass;
  hid_t file, matGroup, lTGroup, coldGroup;
  herr_t status = H5_SUCCESS;
  std::vector<char> materialName;

  file = H5Fopen(filename.c_str(),
                 H5F_ACC_RDONLY, H5P_DEFAULT);
  matGroup = H5Gopen(file, matid_str.c_str(), H5P_DEFAULT);
  lTGroup = H5Gopen(matGroup,   SP5::Depends::logRhoLogT,   H5P_DEFAULT);
  coldGroup = H5Gopen(matGroup, SP5::Depends::coldCurve,    H5P_DEFAULT);

  // This works because chars are 1 byte. ~ JMM
  status += H5LTget_attribute_info(file, matid_str.c_str(),
                                   SP5::Material::name,
                                   &dims,&tclass,&size);
  materialName.resize((int)size+1);
  status += H5LTget_attribute_string(file, matid_str.c_str(),
                                     SP5::Material::name,
                                     materialName.data());
  materialName_ = materialName.data();

  status += loadDataboxes_(matid_str, file, lTGroup, coldGroup);

  status += H5Gclose(lTGroup);
  status += H5Gclose(coldGroup);
  status += H5Gclose(matGroup);
  status += H5Fclose(file);

  if (status != H5_SUCCESS) {
    EOS_ERROR("SpinerDependsRHoT: HDF5 error\n"); //TODO: make this better
  }
}

SpinerEOSDependsRhoT::SpinerEOSDependsRhoT(const std::string& filename,
                                           const std::string& materialName,
                                           bool reproducibility_mode)
  : filename_(filename.c_str())
  , materialName_(materialName.c_str())
  , reproducible_(reproducibility_mode)
  , status_(RootFinding1D::Status::SUCCESS)
  , memoryStatus_(DataStatus::OnHost) {

  std::string matid_str;
  hid_t file, matGroup, lTGroup, coldGroup;
  herr_t status = H5_SUCCESS;

  file = H5Fopen(filename.c_str(),
                 H5F_ACC_RDONLY, H5P_DEFAULT);
  matGroup = H5Gopen(file, materialName.c_str(), H5P_DEFAULT);
  lTGroup = H5Gopen(matGroup,   SP5::Depends::logRhoLogT,   H5P_DEFAULT);
  coldGroup = H5Gopen(matGroup, SP5::Depends::coldCurve,    H5P_DEFAULT);

  status += H5LTget_attribute_int(file, materialName.c_str(),
                                  SP5::Material::matid,
                                  &matid_);
  matid_str = std::to_string(matid_);

  status += loadDataboxes_(matid_str, file, lTGroup, coldGroup);

  status += H5Gclose(lTGroup);
  status += H5Gclose(coldGroup);
  status += H5Gclose(matGroup);
  status += H5Fclose(file);

  if (status != H5_SUCCESS) {
    EOS_ERROR("SpinerDependsRhoT: HDF5 error\n");
  }
}

SpinerEOSDependsRhoT SpinerEOSDependsRhoT::GetOnDevice() {
  SpinerEOSDependsRhoT other;
  other.P_ = Spiner::getOnDeviceDataBox(P_);
  other.sie_ = Spiner::getOnDeviceDataBox(sie_);
  other.bMod_ = Spiner::getOnDeviceDataBox(bMod_);
  other.dPdRho_ = Spiner::getOnDeviceDataBox(dPdRho_);
  other.dPdE_ = Spiner::getOnDeviceDataBox(dPdE_);
  other.dTdRho_ = Spiner::getOnDeviceDataBox(dTdRho_);
  other.dTdE_ = Spiner::getOnDeviceDataBox(dTdE_);
  other.dEdRho_ = Spiner::getOnDeviceDataBox(dEdRho_);
  other.dEdT_ = Spiner::getOnDeviceDataBox(dEdT_);
  other.PMax_ = Spiner::getOnDeviceDataBox(PMax_);
  other.sielTMax_ = Spiner::getOnDeviceDataBox(sielTMax_);
  other.dEdTMax_ = Spiner::getOnDeviceDataBox(dEdTMax_);
  other.gm1Max_ = Spiner::getOnDeviceDataBox(gm1Max_);
  other.PCold_ = Spiner::getOnDeviceDataBox(PCold_);
  other.sieCold_ = Spiner::getOnDeviceDataBox(sieCold_);
  other.bModCold_ = Spiner::getOnDeviceDataBox(bModCold_);
  other.dPdRhoCold_ = Spiner::getOnDeviceDataBox(dPdRhoCold_);
  other.dPdECold_ = Spiner::getOnDeviceDataBox(dPdECold_);
  other.dTdRhoCold_ = Spiner::getOnDeviceDataBox(dTdRhoCold_);
  other.dTdECold_ = Spiner::getOnDeviceDataBox(dTdECold_);
  other.dEdTCold_ = Spiner::getOnDeviceDataBox(dEdTCold_);
  other.lTColdCrit_ = Spiner::getOnDeviceDataBox(lTColdCrit_);
  other.lRhoMin_ = lRhoMin_;
  other.lRhoMax_ = lRhoMax_;
  other.rhoMax_ = rhoMax_;
  other.lRhoMinSearch_ = lRhoMinSearch_;
  other.lTMin_ = lTMin_;
  other.lTMax_ = lTMax_;
  other.TMax_ = TMax_;
  other.lRhoOffset_ = lRhoOffset_;
  other.lTOffset_ = lTOffset_;
  other.rhoNormal_ = rhoNormal_;
  other.TNormal_ = TNormal_;
  other.sieNormal_ = sieNormal_;
  other.PNormal_ = PNormal_;
  other.CvNormal_ = CvNormal_;
  other.bModNormal_ = bModNormal_;
  other.dPdENormal_ = dPdENormal_;
  other.dVdTNormal_ = dVdTNormal_;
  other.numRho_ = numRho_;
  other.numT_ = numT_;
  other.matid_ = matid_;
  other.reproducible_ = reproducible_;
  other.status_ = status_;
  other.memoryStatus_ = DataStatus::OnDevice;
  return other;
}

void SpinerEOSDependsRhoT::Finalize() {
  P_.finalize();
  sie_.finalize();
  bMod_.finalize();
  dPdRho_.finalize();
  dPdE_.finalize();
  dTdRho_.finalize();
  dTdE_.finalize();
  dEdRho_.finalize();
  dEdT_.finalize();
  PMax_.finalize();
  sielTMax_.finalize();
  dEdTMax_.finalize();
  gm1Max_.finalize();
  PCold_.finalize();
  sieCold_.finalize();
  bModCold_.finalize();
  dPdRhoCold_.finalize();
  dPdECold_.finalize();
  dTdRhoCold_.finalize();
  dEdTCold_.finalize();
  dTdECold_.finalize();
  lTColdCrit_.finalize();
  memoryStatus_ = DataStatus::Deallocated;
}

herr_t SpinerEOSDependsRhoT::loadDataboxes_(const std::string& matid_str,
                                            hid_t file, hid_t lTGroup,
                                            hid_t coldGroup) {
  herr_t status = H5_SUCCESS;

  // offsets
  status += H5LTget_attribute_double(file, matid_str.c_str(),
                                     SP5::Offsets::rho,
                                     &lRhoOffset_);
  status += H5LTget_attribute_double(file, matid_str.c_str(),
                                     SP5::Offsets::T,
                                     &lTOffset_);
  lRhoOffset_ = std::abs(lRhoOffset_);
  lTOffset_   = std::abs(lTOffset_);
  // normal density
  status += H5LTget_attribute_double(file, matid_str.c_str(),
                                     SP5::Material::normalDensity,
                                     &rhoNormal_);
  rhoNormal_ = std::abs(rhoNormal_);

  // tables
  status += P_.loadHDF(lTGroup,      SP5::Fields::P);
  status += sie_.loadHDF(lTGroup,    SP5::Fields::sie);
  status += bMod_.loadHDF(lTGroup,   SP5::Fields::bMod);
  status += dPdRho_.loadHDF(lTGroup, SP5::Fields::dPdRho);
  status += dPdE_.loadHDF(lTGroup,   SP5::Fields::dPdE);
  status += dTdRho_.loadHDF(lTGroup, SP5::Fields::dTdRho);
  status += dTdE_.loadHDF(lTGroup,   SP5::Fields::dTdE);
  status += dEdRho_.loadHDF(lTGroup, SP5::Fields::dEdRho);
  status += dEdT_.loadHDF(lTGroup,   SP5::Fields::dEdT);

  // cold curves
  status += PCold_.loadHDF(coldGroup,      SP5::Fields::P);
  status += sieCold_.loadHDF(coldGroup,    SP5::Fields::sie);
  status += bModCold_.loadHDF(coldGroup,   SP5::Fields::bMod);
  status += dPdRhoCold_.loadHDF(coldGroup, SP5::Fields::dPdRho);

  numRho_ = bMod_.dim(2);
  numT_   = bMod_.dim(1);

  // set bounds
  lRhoMin_ = P_.range(1).min();
  lRhoMax_ = P_.range(1).max();
  rhoMax_  = fromLog_(lRhoMax_,lRhoOffset_);
  lTMin_   = P_.range(0).min();
  lTMax_   = P_.range(0).max();
  TMax_    = fromLog_(lTMax_,lTOffset_);

  Real rhoMin = fromLog_(lRhoMin_,lRhoOffset_);
  Real rhoMinSearch = std::max(rhoMin,
                               std::max(std::abs(EPS)*10,
                                        std::abs(EPS*rhoMin)));
  lRhoMinSearch_ = toLog_(rhoMinSearch,lRhoOffset_);

  // bulk modulus can be wrong in the tables. Use FLAG's approach to
  // fix the table.
  fixBulkModulus_();

  // find critical temperature Tcrit(rho)
  // where sie(rho,Tcrit(rho)) = sieCold(rho)
  setlTColdCrit_();

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
    bModCold_(j)   = bMod_.interpToReal(lRho,lT);
    dPdECold_(j)   = dPdE_.interpToReal(lRho,lT);
    dTdRhoCold_(j) = dTdRho_.interpToReal(lRho,lT);
    dTdECold_(j)   = dTdE_.interpToReal(lRho,lT);
    dEdTCold_(j)   = dEdT_.interpToReal(lRho,lT);
  }
  
  // major vs. minor axes change, so this must be done by hand
  PMax_.resize(numRho_);
  PMax_.setRange(0, bMod_.range(1));
  dEdTMax_.copyMetadata(PMax_);
  gm1Max_.copyMetadata(PMax_);
  sielTMax_.copyMetadata(PMax_);
  for (int j = 0; j < numRho_; j++) {
    Real rho = PMax_.range(0).x(j);
    PMax_(j) = P_(j,numT_-1);
    dEdTMax_(j) = dEdT_(j, numT_-1);
    gm1Max_(j) = dPdE_(j, numT_-1)/(rho + EPS); // max gruneisen
    sielTMax_(j) = sie_(j, numT_-1);
  }

  // reference state
  Real lRhoNormal = lRho_(rhoNormal_);
  // if rho normal not on the table, set it to the middle
  if ( !(lRhoMin_ < lRhoNormal && lRhoNormal < lRhoMax_) ) {
    lRhoNormal = 0.5*(lRhoMin_ + lRhoMax_);
    rhoNormal_ = rho_(lRhoNormal);
  }
  // Same for temperature. Use room temperature if it's available
  TNormal_ = ROOM_TEMPERATURE;
  Real lTNormal = lT_(TNormal_);
  if ( !(lTMin_ < lTNormal && lTNormal < lTMax_) ) {
    lTNormal = 0.5*(lTMin_ + lTMax_);
    TNormal_ = T_(lTNormal);
  }
  sieNormal_ = sie_.interpToReal(lRhoNormal,lTNormal);
  PNormal_ = P_.interpToReal(lRhoNormal,lTNormal);
  CvNormal_ = dEdT_.interpToReal(lRhoNormal,lTNormal);
  bModNormal_ = bMod_.interpToReal(lRhoNormal,lTNormal);
  dPdENormal_ = dPdE_.interpToReal(lRhoNormal,lTNormal);
  Real dPdR = dPdRho_.interpToReal(lRhoNormal,lTNormal);
  dVdTNormal_ = dPdENormal_*CvNormal_/(rhoNormal_*rhoNormal_*dPdR);
  
  return status;
}

void SpinerEOSDependsRhoT::fixBulkModulus_() {
  // assumes all databoxes are the same size
  // TODO: do we need to smooth this data with a median filter
  // or something like that?
  for (int j = 0; j < numRho_; j++) {
    Real lRho = bMod_.range(1).x(j);
    Real rho = rho_(lRho);
    for (int i = 0; i < numT_; i++) {
      Real lT = bMod_.range(0).x(i);
      Real press = P_.interpToReal(lRho,lT);
      Real DPDR = dPdRho_.interpToReal(lRho,lT);
      Real DPDE = dPdE_.interpToReal(lRho,lT);
      Real DEDR = dEdRho_.interpToReal(lRho,lT);
      Real DTDE = dTdE_.interpToReal(lRho,lT);
      Real bMod;
      if (DPDE > 0.0 && rho > 0.0) {
        bMod = rho*DPDR + DPDE*(press/rho - rho*DEDR);
      } else if (rho > 0.0) {
        bMod = std::max(rho*DPDR,0.0);
      } else {
        bMod = 0.0;
      }
      bMod_(j,i) = std::max(bMod,std::abs(EPS));
    }
  }
}

void SpinerEOSDependsRhoT::setlTColdCrit_() {
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
    for (int i = 0; i < numT_-1; i++) {
      Real sign_curr = sie_(j,i) - sieCold;
      Real sign_next = sie_(j,i+1) - sieCold;
      if ( sign_curr < 0 && sign_next > 0 ) {
        crossings.push_back(i);
        directions.push_back(1);
        last_pos_crossing = crossings.size()-1;
      } else if ( sign_curr > 0 && sign_next < 0 ) {
        crossings.push_back(i);
        directions.push_back(-1);
      }
    }
    // react accordingly
    // TODO(JMM): If we have multiple zero crossings,
    // there's no good uniquely defined solution.
    // We choose the highest-temperature crossing,
    // which results in the least code. But this may be wrong.
    if (last_pos_crossing < 0) { // off the grid
      lTColdCrit_(j) = lTMin_;
    } else { // at least one pos crossing. Use last one.
      const callable_interp::r_interp sieFunc(sie_, lRho);
      Real lT;
      int ilast = crossings[last_pos_crossing];
      Real lTlower = bMod_.range(0).x(ilast);
      Real lTupper = bMod_.range(0).x(ilast + 1);
      Real lTGuess = 0.5*(lTlower + lTupper);
      auto status = findRoot(sieFunc, sieCold, lTGuess,
                             lTlower, lTupper,
                             ROOT_THRESH, ROOT_THRESH,
                             lT, counts);
      if (status != RootFinding1D::Status::SUCCESS) {
        lT = lTGuess;
      }
      lTColdCrit_(j) = lT;
      // lTColdCrit_(j) = lTMin_;
    }
    lTColdCrit_(j) = std::max(lTMin_,lTColdCrit_(j));
  }
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoT::TemperatureFromDensityInternalEnergy(const Real rho,
                                                                const Real sie,
                                                                Real *lambda) const {
  TableStatus whereAmI;
  const Real lRho = lRho_(rho);
  const Real lT = lTFromlRhoSie_(lRho,sie,whereAmI,lambda);
  return T_(lT);
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoT::InternalEnergyFromDensityTemperature(const Real rho,
                                                                const Real temperature,
                                                                Real *lambda) const {
  Real lRho, lT;
  getLogsRhoT_(rho, temperature, lRho, lT, lambda);
  TableStatus whereAmI = getLocDependsRhoT_(lRho,lT);
  return sieFromlRhoTlT_(lRho,temperature,lT,whereAmI);
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoT::PressureFromDensityTemperature(const Real rho,
                                                          const Real temperature,
                                                          Real *lambda) const {
  Real lRho, lT;
  getLogsRhoT_(rho, temperature, lRho, lT, lambda);
  TableStatus whereAmI = getLocDependsRhoT_(lRho,lT);
  return PFromRholRhoTlT_(rho,lRho,temperature,lT,whereAmI);
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoT::PressureFromDensityInternalEnergy(const Real rho,
                                                             const Real sie,
                                                             Real *lambda) const {
  TableStatus whereAmI;
  Real lRho = lRho_(rho);
  Real lT = lTFromlRhoSie_(lRho,sie,whereAmI,lambda);
  Real P;
  if ( whereAmI == TableStatus::OffBottom ) { // cold curve
    P = PCold_.interpToReal(lRho);
  } else if ( whereAmI == TableStatus::OffTop ) { // ideal gas
    const Real gm1 = gm1Max_.interpToReal(lRho);
    P = gm1*rho*sie;
  } else { // on table
    P = P_.interpToReal(lRho,lT);
  }
  return P;
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoT::SpecificHeatFromDensityTemperature(const Real rho,
                                                              const Real temperature,
                                                              Real *lambda) const {
  Real lRho, lT;
  getLogsRhoT_(rho, temperature, lRho, lT, lambda);
  TableStatus whereAmI = getLocDependsRhoT_(lRho,lT);
  return CvFromlRholT_(lRho,lT,whereAmI);
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoT::SpecificHeatFromDensityInternalEnergy(const Real rho,
                                                                 const Real sie,
                                                                 Real *lambda) const {
  TableStatus whereAmI;
  Real Cv;
  const Real lRho = lRho_(rho);
  const Real lT = lTFromlRhoSie_(lRho,sie,whereAmI,lambda);
  if ( whereAmI == TableStatus::OffBottom ) { // cold curve
    // on cold curve. Currently, we assume constant extrapolation.
    // TODO(JMM): Do something more sophisticated
    Cv = dEdTCold_.interpToReal(lRho);
  } else if ( whereAmI == TableStatus::OffTop ) { // ideal gas
    Cv = dEdTMax_.interpToReal(lRho); // Cv assumed constant in T
  } else { // on table
    Cv = dEdT_.interpToReal(lRho,lT);
  }
  return Cv > EPS ? Cv : EPS;
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoT::BulkModulusFromDensityTemperature(const Real rho,
                                                             const Real temperature,
                                                             Real *lambda) const 
{
  Real lRho, lT;
  getLogsRhoT_(rho, temperature, lRho, lT, lambda);
  TableStatus whereAmI = getLocDependsRhoT_(lRho,lT);
  return bModFromRholRhoTlT_(rho,lRho,temperature,lT,whereAmI);
}

PORTABLE_FUNCTION
Real
SpinerEOSDependsRhoT::GruneisenParamFromDensityTemperature(const Real rho,
                                                           const Real temp,
                                                           Real *lambda) const 
{
  Real lRho, lT, gm1;
  getLogsRhoT_(rho, temp, lRho, lT, lambda);

  TableStatus whereAmI = getLocDependsRhoT_(lRho, lT);
  if ( whereAmI == TableStatus::OffBottom ) {
    // use cold curves
    Real dpde = dPdECold_.interpToReal(lRho);
    gm1 = std::abs(dpde)/(std::abs(rho) + EPS);
  } else if ( whereAmI == TableStatus::OffTop ) {
    gm1 = gm1Max_.interpToReal(lRho);
  } else { // on table
    const Real dpde = dPdE_.interpToReal(lRho,lT);
    gm1 = std::abs(dpde)/(std::abs(rho) + EPS);
  }
  return gm1;
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoT::BulkModulusFromDensityInternalEnergy(const Real rho,
                                                                const Real sie,
                                                                Real *lambda) const {
  TableStatus whereAmI;
  Real bMod;
  const Real lRho = lRho_(rho);
  const Real lT = lTFromlRhoSie_(lRho, sie, whereAmI, lambda);
  if ( whereAmI == TableStatus::OffBottom) {
    bMod = bModCold_.interpToReal(lRho);
  } else if ( whereAmI == TableStatus::OffTop ) {
    const Real gm1 = gm1Max_.interpToReal(lRho);
    bMod = (gm1 + 1)*gm1*rho*sie;
  } else { // on table
    bMod = bMod_.interpToReal(rho,sie);
  }
  return bMod;
}

PORTABLE_FUNCTION
Real
SpinerEOSDependsRhoT::GruneisenParamFromDensityInternalEnergy(const Real rho,
                                                              const Real sie,
                                                              Real *lambda) const {
  TableStatus whereAmI;
  Real gm1;
  const Real lRho = lRho_(rho);
  const Real lT = lTFromlRhoSie_(lRho, sie, whereAmI, lambda);
  if ( whereAmI == TableStatus::OffBottom ) {
    Real dpde = dPdECold_.interpToReal(lRho);
    gm1 = std::abs(dpde)/(std::abs(rho) + EPS);
  } else if ( whereAmI == TableStatus::OffTop ) {
    gm1 = gm1Max_.interpToReal(lRho);
  } else {
    const Real dpde = dPdE_.interpToReal(lRho,lT);
    gm1 = std::abs(dpde)/(std::abs(rho) + EPS);
  }
  return gm1;
}

// TODO(JMM): This would be faster with hand-tuned code
PORTABLE_FUNCTION
void SpinerEOSDependsRhoT::DensityEnergyFromPressureTemperature(const Real press,
                                                                const Real temp,
                                                                Real *lambda,
                                                                Real& rho,
                                                                Real& sie) const {
  TableStatus whereAmI;
  Real lT = lT_(temp);
  Real lRho = lRhoFromPlT_(press, lT, whereAmI, lambda);
  rho = rho_(lRho);
  sie = InternalEnergyFromDensityTemperature(rho, temp, lambda);
}

PORTABLE_FUNCTION
void SpinerEOSDependsRhoT::PTofRE(const Real rho, const Real sie,
                                  Real * lambda, Real& press, Real& temp,
                                  Real & dpdr, Real & dpde,
                                  Real & dtdr, Real & dtde) const {
  TableStatus whereAmI;
  const Real lRho = lRho_(rho);
  const Real lT = lTFromlRhoSie_(lRho, sie, whereAmI, lambda);
  temp = T_(lT);
  if ( whereAmI == TableStatus::OffBottom ) {
    press = PCold_.interpToReal(lRho);
    dpdr = dPdRhoCold_.interpToReal(lRho);
    dpde = dPdECold_.interpToReal(lRho);
    dtdr = dTdRhoCold_.interpToReal(lRho);
    dtde = dTdECold_.interpToReal(lRho);
  } else if ( whereAmI == TableStatus::OffTop ) {
    const Real Cv = dEdTMax_.interpToReal(lRho);
    const Real gm1 = gm1Max_.interpToReal(lRho);
    const Real e0 = sielTMax_.interpToReal(lRho);
    const Real e = e0 + Cv*(temp - TMax_);
    press = gm1*rho*e;
    dpdr =  gm1*e;
    dpde = gm1*rho;
    dtdr = -press/(gm1*Cv*rho*rho);
    dtde = 1./Cv;
  } else { // on table
    press = P_.interpToReal(lRho,lT);
    dpdr = dPdRho_.interpToReal(lRho,lT);
    dpde = dPdE_.interpToReal(lRho,lT);
    dtdr = dTdRho_.interpToReal(lRho,lT);
    dtde = dTdE_.interpToReal(lRho,lT);
  }
}

PORTABLE_FUNCTION
void SpinerEOSDependsRhoT::FillEos(Real& rho, Real& temp,
                                   Real& energy, Real& press,
                                   Real& cv, Real& bmod,
                                   const unsigned long output,
                                   Real *lambda) const {
  Real lRho, lT;
  TableStatus whereAmI;
  const unsigned long input = ~output;
  if (output == thermalqs::none) {
    UNDEFINED_ERROR;
  }
  if (output&thermalqs::density) {
    if (input&thermalqs::pressure && input&thermalqs::temperature) {
      lT = lT_(temp);
      lRho = lRhoFromPlT_(press,lT,whereAmI,lambda);
      rho = rho_(lRho);
    } else {
      UNDEFINED_ERROR;
    }
  } else {
    lRho = lRho_(rho);
  }
  if (output&thermalqs::temperature) {
    if (input&thermalqs::density && input&thermalqs::specific_internal_energy) {
      lT = lTFromlRhoSie_(lRho, energy, whereAmI, lambda);
    } else if (input&thermalqs::density && input&thermalqs::pressure) {
      lT = lTFromlRhoP_(lRho, press, whereAmI, lambda);
     }
    else {
      UNDEFINED_ERROR;
    }
    temp = T_(lT);
  } else {
    lT = lT_(temp);
  }
  whereAmI = getLocDependsRhoT_(lRho,lT);
  if (output&thermalqs::specific_internal_energy) {
    energy = sieFromlRhoTlT_(lRho,temp,lT,whereAmI);
  }
  if (output&thermalqs::pressure) {
    press = PFromRholRhoTlT_(rho,lRho,temp,lT,whereAmI);
  }
  if (output&thermalqs::specific_heat) {
    cv = CvFromlRholT_(lRho,lT,whereAmI);
  }
  if (output&thermalqs::bulk_modulus) {
    bmod = bModFromRholRhoTlT_(rho,lRho,temp,lT,whereAmI);
  }
  if (lambda != nullptr) {
    lambda[Lambda::lRho] = lRho;
    lambda[Lambda::lT]   = lT;
  }
}

PORTABLE_FUNCTION
void SpinerEOSDependsRhoT::ValuesAtReferenceState(Real& rho, Real& temp,
                                                  Real& sie, Real& press,
                                                  Real& cv, Real& bmod,
                                                  Real& dpde, Real& dvdt,
                                                  Real *lambda) const {
  rho = rhoNormal_;
  temp = TNormal_;
  sie = sieNormal_;
  press = PNormal_;
  cv = CvNormal_;
  bmod = bModNormal_;
  dpde = dPdENormal_;
  dvdt = dVdTNormal_;
}

PORTABLE_FUNCTION
void SpinerEOSDependsRhoT::getLogsRhoT_(const Real rho,
                                        const Real temperature,
                                        Real& lRho, Real& lT,
                                        Real* lambda) const {
  lRho = lRho_(rho);
  lT = lT_(temperature);
  if (lambda != nullptr) {
    lambda[Lambda::lRho] = lRho;
    lambda[Lambda::lT]   = lT;
  }
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoT::lRhoFromPlT_(const Real P,
                                        const Real lT,
                                        TableStatus& whereAmI,
                                        Real* lambda) const {
  using RootFinding1D::findRoot;
  RootFinding1D::Status status = RootFinding1D::Status::SUCCESS;

  Real lRho;
  Real lRhoGuess = reproducible_ ? lRhoMax_ : 0.5*(lRhoMin_ + lRhoMax_);
  // Real lRhoGuess = lRhoMin_ + 0.9*(lRhoMax_ - lRhoMin_);
  if (lambda != nullptr
      && lRhoMin_ <= lambda[Lambda::lRho]
      && lambda[Lambda::lRho] <= lRhoMax_) {
    lRhoGuess = lambda[Lambda::lRho];
  }
  if ( lT <= lTMin_ ) { // cold curve
    whereAmI = TableStatus::OffBottom;
    const callable_interp::interp PFunc(PCold_);
    status = findRoot(PFunc,
                      P, lRhoGuess,
                      //lRhoMin_, lRhoMax_,
                      lRhoMinSearch_, lRhoMax_,
                      ROOT_THRESH, ROOT_THRESH,
                      lRho, counts);
  } else if ( lT >= lTMax_ ) { // ideal gas
    whereAmI = TableStatus::OffTop;
    const callable_interp::prod_interp_1d PFunc(gm1Max_,dEdTMax_,lT);
    status = findRoot(PFunc,
                      P, lRhoGuess,
                      //lRhoMin_, lRhoMax_,
                      lRhoMinSearch_, lRhoMax_,
                      ROOT_THRESH, ROOT_THRESH,
                      lRho, counts);
  } else { // on table
    whereAmI = TableStatus::OnTable;
    const callable_interp::l_interp PFunc(P_, lT);
    status = findRoot(PFunc,
                      P, lRhoGuess,
                      //lRhoMin_, lRhoMax_,
                      lRhoMinSearch_, lRhoMax_,
                      ROOT_THRESH, ROOT_THRESH,
                      lRho, counts);
  }
  if (status_ != RootFinding1D::Status::SUCCESS) {
#if EPINER_EOS_VERBOSE
    std::stringstream errorMessage;
    errorMessage << "inverting P table for logRho failed\n"
                 << "matid     = " << matid_ << "\n"
                 << "lT        = " << lT << "\n"
                 << "P         = " << P << "\n"
                 << "lRhoGuess = " << lRhoGuess
                 << std::endl;
    EOS_ERROR(errorMessage.str().c_str());
#endif //SPINER_EOS_VERBOSE
    lRho = reproducible_ ? lRhoMax_ : lRhoGuess;
  }
  if (lambda != nullptr) {
    lambda[Lambda::lRho]   = lRho;
    lambda[Lambda::lT]     = lT;
  }
  status_ = status;
  whereAmI_ = whereAmI;
  return lRho;
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoT::lTFromlRhoSie_(const Real lRho,
                                          const Real sie,
                                          TableStatus& whereAmI,
                                          Real* lambda) const {

  RootFinding1D::Status status = RootFinding1D::Status::SUCCESS;
  using RootFinding1D::findRoot;
  Real lT;

  whereAmI = getLocDependsRhoSie_(lRho,sie);
  if ( whereAmI == TableStatus::OffBottom ) {
    // On the cold curve. No unique temperature.
    // return the minimum temperature in the table.
    lT = lTMin_; 
    counts.increment(0);
  } else if ( whereAmI == TableStatus::OffTop ) { // Assume ideal gas
    const Real Cv = dEdTMax_.interpToReal(lRho);
    const Real e0 = sielTMax_.interpToReal(lRho);
    const Real T = TMax_ + (sie - e0)/(Cv+EPS);
    lT = lT_(T);
    counts.increment(0);
  } else {
    Real lTGuess = reproducible_ ? lTMin_ : 0.5*(lTMin_ + lTMax_);
    if (lambda != nullptr
        && lTMin_ <= lambda[Lambda::lT] && lambda[Lambda::lT] <= lTMax_) {
      lTGuess = lambda[Lambda::lT];
    }
    const callable_interp::r_interp sieFunc(sie_, lRho);
    status = findRoot(sieFunc,
                      sie, lTGuess,
                      lTMin_, lTMax_,
                      ROOT_THRESH, ROOT_THRESH,
                      lT, counts);

    if (status != RootFinding1D::Status::SUCCESS) {
#if SPINER_EOS_VERBOSE
      std::stringstream errorMessage;
      errorMessage << "inverting sie table for logT failed\n"
                   << "matid   = " << matid_ << "\n"
                   << "lRho    = " << lRho << "\n"
                   << "sie     = " << sie << "\n"
                   << "lTGuess = " << lTGuess
                   << std::endl;
      EOS_ERROR(errorMessage.str().c_str());
#endif // SPINER_EOS_VERBOSE
      lT = reproducible_ ? lTMin_ : lTGuess;
    }
  }
  if (lambda != nullptr) {
    lambda[Lambda::lRho]   = lRho;
    lambda[Lambda::lT]     = lT;
  }
  status_ = status;
  whereAmI_ = whereAmI;
  return lT;
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoT::lTFromlRhoP_(const Real lRho,
                                        const Real press,
                                        TableStatus& whereAmI,
                                        Real* lambda) const {
  RootFinding1D::Status status = RootFinding1D::Status::SUCCESS;
  using RootFinding1D::findRoot;
  Real lT, lTGuess;

  // Assumes P is monotone in T
  const Real PCold = PCold_.interpToReal(lRho);
  const Real PMax = PMax_.interpToReal(lRho);
  if ( press <= PCold ) {
    whereAmI = TableStatus::OffBottom;
    lT = lTMin_;
    counts.increment(0);
  } else if ( press >= PMax ) {
    whereAmI = TableStatus::OffTop;
    lT = lTMax_;
    counts.increment(0);
  } else {
    whereAmI = TableStatus::OnTable;
    if (lambda != nullptr
        && lTMin_ <= lambda[Lambda::lT]
        && lambda[Lambda::lT] <= lTMax_) {
      lTGuess = lambda[Lambda::lT];
    } else {
      lTGuess = 0.5*(lTMin_ + lTMax_);
    }
    const callable_interp::r_interp PFunc(P_, lRho);
    status = findRoot(PFunc,
                      press, lTGuess,
                      lTMin_, lTMax_,
                      ROOT_THRESH, ROOT_THRESH,
                      lT, counts);
    if (status != RootFinding1D::Status::SUCCESS) {
#if SPINER_EOS_VERBOSE
      std::stringstream errorMessage;
      errorMessage << "inverting P table for logT failed\n"
                   << "matid   = " << matid_ << "\n"
                   << "lRho    = " << lRho << "\n"
                   << "P       = " << press << "\n"
                   << "lTGuess = " << lTGuess
                   << std::endl;
      EOS_ERROR(errorMessage.str().c_str());
#endif // SPINER_EOS_VERBOSE
      lT = reproducible_ ? lTMin_ : lTGuess;
    }
  }
  if (lambda != nullptr) {
    lambda[Lambda::lRho]   = lRho;
    lambda[Lambda::lT]     = lT;
  }
  status_ = status;
  whereAmI_ = whereAmI;
  return lT;
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoT::sieFromlRhoTlT_(const Real lRho,
                                           const Real T, const Real lT,
                                           const TableStatus& whereAmI) const {
  Real sie;
  if ( whereAmI == TableStatus::OffBottom ) { // cold curve
    sie = sieCold_.interpToReal(lRho);
  } else if ( whereAmI == TableStatus::OffTop ) { // ideal gas
    const Real Cv = dEdTMax_.interpToReal(lRho);
    const Real e0 = sielTMax_.interpToReal(lRho);
    sie = e0 + Cv*(T - TMax_);
  } else { // on table
    sie = sie_.interpToReal(lRho, lT);
  }
  return sie;
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoT::PFromRholRhoTlT_(const Real rho, const Real lRho,
                                            const Real T, const Real lT,
                                            const TableStatus& whereAmI) const {
  Real P;
  if ( whereAmI == TableStatus::OffBottom ) { // cold curve
    P =  PCold_.interpToReal(lRho);
  } else if ( whereAmI == TableStatus::OffTop ) { // ideal gas
    const Real Cv = dEdTMax_.interpToReal(lRho);
    const Real gm1 = gm1Max_.interpToReal(lRho);
    const Real e0 = sielTMax_.interpToReal(lRho);
    const Real e = e0 + Cv*(T - TMax_);
    P = gm1*rho*e;
  } else { // if ( whereAmI == TableStatus::OnTable) {
    P = P_.interpToReal(lRho, lT);
  }
  return P;
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoT::CvFromlRholT_(const Real lRho, const Real lT,
                                         const TableStatus& whereAmI) const {
  Real Cv;
  if ( whereAmI == TableStatus::OffBottom ) { // cold curve
    // on cold curve. Currently, we assume constant extrapolation.
    // TODO(JMM): Do something more sophisticated
    Cv = dEdTCold_.interpToReal(lRho);
  } else if ( whereAmI == TableStatus::OffTop ) { // ideal gas
    Cv = dEdTMax_.interpToReal(lRho); // Cv assumed constant in T
  } else { // on table
    Cv = dEdT_.interpToReal(lRho,lT);
  }
  return Cv > EPS ? Cv : EPS;
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoT::bModFromRholRhoTlT_(const Real rho, const Real lRho,
                                               const Real T, const Real lT,
                                               const TableStatus& whereAmI) const {
  Real bMod;
  if ( whereAmI == TableStatus::OffBottom ) {
    bMod = bModCold_.interpToReal(lRho);
  } else if ( whereAmI == TableStatus::OffTop ) {
    const Real Cv = dEdTMax_.interpToReal(lRho);
    const Real gm1 = gm1Max_.interpToReal(lRho);
    const Real e0 = sielTMax_.interpToReal(lRho);
    const Real e = e0 + Cv*(T - TMax_);
    bMod = (gm1 + 1)*gm1*rho*e;
  } else { // on table
    bMod = bMod_.interpToReal(lRho,lT);
  }
  return bMod > EPS ? bMod : EPS;
}

PORTABLE_FUNCTION TableStatus
SpinerEOSDependsRhoT::getLocDependsRhoSie_(const Real lRho,
                                           const Real sie) const {
  TableStatus whereAmI;
  if ( sie >= sielTMax_.interpToReal(lRho) ) {
    whereAmI = TableStatus::OffTop;
  } else if ( sie <= sieCold_.interpToReal(lRho) ) {
    whereAmI = TableStatus::OffBottom;
  } else {
    whereAmI = TableStatus::OnTable;
  }
  whereAmI_ = whereAmI;
  return whereAmI;
}

PORTABLE_FUNCTION TableStatus
SpinerEOSDependsRhoT::getLocDependsRhoT_(const Real lRho, const Real lT) const {
  TableStatus whereAmI;
  if ( lT <= lTMin_ ) whereAmI = TableStatus::OffBottom;
  else if ( lT >= lTMax_ ) whereAmI = TableStatus::OffTop;
  else whereAmI = TableStatus::OnTable;
  whereAmI_ = whereAmI;
  return whereAmI;
}

SpinerEOSDependsRhoSie::SpinerEOSDependsRhoSie(const std::string& filename,
                                               int matid,
                                               bool reproducibility_mode)
  : filename_(filename.c_str())
  , matid_(matid)
  , reproducible_(reproducibility_mode)
  , status_(RootFinding1D::Status::SUCCESS)
  , memoryStatus_(DataStatus::OnHost) {

  std::string matid_str = std::to_string(matid);
  hsize_t dims;
  size_t size;
  H5T_class_t tclass;
  hid_t file, matGroup, lTGroup, lEGroup;
  herr_t status = H5_SUCCESS;
  std::vector<char> materialName;

  file = H5Fopen(filename.c_str(),
                 H5F_ACC_RDONLY, H5P_DEFAULT);
  matGroup = H5Gopen(file, matid_str.c_str(), H5P_DEFAULT);
  lTGroup = H5Gopen(matGroup, SP5::Depends::logRhoLogT,   H5P_DEFAULT);
  lEGroup = H5Gopen(matGroup, SP5::Depends::logRhoLogSie, H5P_DEFAULT);

  // This works because chars are 1 byte. ~ JMM
  status += H5LTget_attribute_info(file, matid_str.c_str(),
                                   SP5::Material::name,
                                   &dims,&tclass,&size);
  materialName.resize((int)size+1);
  status += H5LTget_attribute_string(file, matid_str.c_str(),
                                     SP5::Material::name,
                                     materialName.data());
  materialName_ = materialName.data();

  status += loadDataboxes_(matid_str, file, lTGroup, lEGroup);

  status += H5Gclose(lTGroup);
  status += H5Gclose(lEGroup);
  status += H5Gclose(matGroup);
  status += H5Fclose(file);
  
  if (status != H5_SUCCESS) {
    EOS_ERROR("SpinerDependsTHoSIE: HDF5 error\n");
  }
}

SpinerEOSDependsRhoSie::SpinerEOSDependsRhoSie(const std::string& filename,
                                               const std::string& materialName,
                                               bool reproducibility_mode)
  : filename_(filename.c_str())
  , materialName_(materialName.c_str())
  , reproducible_(reproducibility_mode)
  , status_(RootFinding1D::Status::SUCCESS)
  , memoryStatus_(DataStatus::OnHost) {

  std::string matid_str;
  hid_t file, matGroup, lTGroup, lEGroup;
  herr_t status = H5_SUCCESS;

  file = H5Fopen(filename.c_str(),
                 H5F_ACC_RDONLY, H5P_DEFAULT);
  matGroup = H5Gopen(file, materialName.c_str(), H5P_DEFAULT);
  lTGroup = H5Gopen(matGroup, SP5::Depends::logRhoLogT,   H5P_DEFAULT);
  lEGroup = H5Gopen(matGroup, SP5::Depends::logRhoLogSie, H5P_DEFAULT);

  status += H5LTget_attribute_int(file, materialName.c_str(),
                                  SP5::Material::matid,
                                  &matid_);
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

herr_t SpinerEOSDependsRhoSie::loadDataboxes_(const std::string& matid_str,
                                              hid_t file,
                                              hid_t lTGroup,
                                              hid_t lEGroup) {
  herr_t status = H5_SUCCESS;

  // offsets
  status += H5LTget_attribute_double(file, matid_str.c_str(),
                                     SP5::Offsets::rho,
                                     &lRhoOffset_);
  status += H5LTget_attribute_double(file, matid_str.c_str(),
                                     SP5::Offsets::T,
                                     &lTOffset_);
  status += H5LTget_attribute_double(file, matid_str.c_str(),
                                     SP5::Offsets::sie,
                                     &lEOffset_);
  lRhoOffset_ = std::abs(lRhoOffset_);
  lTOffset_   = std::abs(lTOffset_);
  lEOffset_   = std::abs(lEOffset_);

  // normal density
  status += H5LTget_attribute_double(file, matid_str.c_str(),
                                     SP5::Material::normalDensity,
                                     &rhoNormal_);
  rhoNormal_ = std::abs(rhoNormal_);

  // sometimes independent variables
  status += sie_.loadHDF(lTGroup, SP5::Fields::sie);
  status += T_.loadHDF(lEGroup,   SP5::Fields::T);

  // dependent variables
  // depends on rho and T
  status += dependsRhoT_.P.loadHDF(lTGroup,       SP5::Fields::P);
  status += dependsRhoT_.bMod.loadHDF(lTGroup,    SP5::Fields::bMod);
  status += dependsRhoT_.dPdRho.loadHDF(lTGroup,  SP5::Fields::dPdRho);
  status += dependsRhoT_.dPdE.loadHDF(lTGroup,    SP5::Fields::dPdE);
  status += dependsRhoT_.dTdRho.loadHDF(lTGroup,  SP5::Fields::dTdRho);
  status += dependsRhoT_.dTdE.loadHDF(lTGroup,    SP5::Fields::dTdE);
  status += dependsRhoT_.dEdRho.loadHDF(lTGroup,  SP5::Fields::dEdRho);
  // depends on rho and e
  status += dependsRhoSie_.P.loadHDF(lEGroup,       SP5::Fields::P);
  status += dependsRhoSie_.bMod.loadHDF(lEGroup,    SP5::Fields::bMod);
  status += dependsRhoSie_.dPdRho.loadHDF(lEGroup,  SP5::Fields::dPdRho);
  status += dependsRhoSie_.dPdE.loadHDF(lEGroup,    SP5::Fields::dPdE);
  status += dependsRhoSie_.dTdRho.loadHDF(lEGroup,  SP5::Fields::dTdRho);
  status += dependsRhoSie_.dTdE.loadHDF(lEGroup,    SP5::Fields::dTdE);
  status += dependsRhoSie_.dEdRho.loadHDF(lEGroup,  SP5::Fields::dEdRho);
  

  // Fix up bulk modulus
  calcBMod_(dependsRhoT_);
  calcBMod_(dependsRhoSie_);

  // Metadata for root finding extrapolation
  numRho_ = sie_.dim(2);
  lRhoMin_ = sie_.range(1).min();
  lRhoMax_ = sie_.range(1).max();
  rhoMax_  = fromLog_(lRhoMax_,lRhoOffset_);

  // slice to maximum of rho
  PlRhoMax_  = dependsRhoT_.P.slice(numRho_-1);
  dPdRhoMax_ = dependsRhoT_.dPdRho.slice(numRho_-1);

  // reference state
  Real lRhoNormal = toLog_(rhoNormal_,lRhoOffset_);
  // if rho normal not on the table, set it to the middle
  if ( !(lRhoMin_ < lRhoNormal && lRhoNormal < lRhoMax_) ) {
    lRhoNormal = 0.5*(lRhoMin_ + lRhoMax_);
    rhoNormal_ = fromLog_(lRhoNormal,lRhoOffset_);
  }
  // Same for temperature. Use room temperature if it's available
  TNormal_ = ROOM_TEMPERATURE;
  Real lTNormal = toLog_(TNormal_,lTOffset_);
  Real lTMin = sie_.range(0).min();
  Real lTMax = sie_.range(0).max();
  if ( !(lTMin < lTNormal && lTNormal < lTMax) ) {
    lTNormal = 0.5*(lTMin + lTMax);
    TNormal_ = fromLog_(lTNormal,lTOffset_);
  }
  sieNormal_ = sie_.interpToReal(lRhoNormal,lTNormal);
  PNormal_ = dependsRhoT_.P.interpToReal(lRhoNormal,lTNormal);
  CvNormal_ = 1./dependsRhoT_.dTdE.interpToReal(lRhoNormal,lTNormal);
  bModNormal_ = dependsRhoT_.bMod.interpToReal(lRhoNormal,lTNormal);
  dPdENormal_ = dependsRhoT_.dPdE.interpToReal(lRhoNormal,lTNormal);
  Real dPdR = dependsRhoT_.dPdRho.interpToReal(lRhoNormal,lTNormal);
  dVdTNormal_ = dPdENormal_*CvNormal_/(rhoNormal_*rhoNormal_*dPdR);

  return status;
}

void SpinerEOSDependsRhoSie::calcBMod_(SP5Tables& tables) {
  for (int j = 0; j < tables.bMod.dim(2); j++) {
    Real lRho = tables.bMod.range(1).x(j);
    Real rho = fromLog_(lRho,lRhoOffset_);
    for (int i = 0; i < tables.bMod.dim(1); i++) {
      Real press = tables.P(j,i);
      Real DPDR  = tables.dPdRho(j,i);
      Real DPDE  = tables.dPdE(j,i);
      Real DEDR  = tables.dEdRho(j,i);
      Real DTDE  = tables.dTdE(j,i);
      Real bMod;
      if (DPDE > 0.0 && rho > 0.0) {
        bMod = rho*DPDR + DPDE*(press/rho - rho*DEDR);
      } else if (rho > 0.0) {
        bMod = std::max(rho*DPDR,0.0);
      } else {
        bMod = 0.0;
      }
      tables.bMod(j,i) = std::max(bMod,std::abs(Spiner::EPS));
    }
  }
}

SpinerEOSDependsRhoSie SpinerEOSDependsRhoSie::GetOnDevice() {
  SpinerEOSDependsRhoSie other;
  using Spiner::getOnDeviceDataBox;
  other.sie_ = getOnDeviceDataBox(sie_);
  other.T_   = getOnDeviceDataBox(T_);
  other.dependsRhoT_.P        = getOnDeviceDataBox(dependsRhoT_.P);
  other.dependsRhoT_.bMod     = getOnDeviceDataBox(dependsRhoT_.bMod);
  other.dependsRhoT_.dPdRho   = getOnDeviceDataBox(dependsRhoT_.dPdRho);
  other.dependsRhoT_.dPdE     = getOnDeviceDataBox(dependsRhoT_.dPdE);
  other.dependsRhoT_.dTdRho   = getOnDeviceDataBox(dependsRhoT_.dTdRho);
  other.dependsRhoT_.dTdE     = getOnDeviceDataBox(dependsRhoT_.dTdE);
  other.dependsRhoT_.dEdRho   = getOnDeviceDataBox(dependsRhoT_.dEdRho);
  other.dependsRhoSie_.P      = getOnDeviceDataBox(dependsRhoSie_.P);
  other.dependsRhoSie_.bMod   = getOnDeviceDataBox(dependsRhoSie_.bMod);
  other.dependsRhoSie_.dPdRho = getOnDeviceDataBox(dependsRhoSie_.dPdRho);
  other.dependsRhoSie_.dPdE   = getOnDeviceDataBox(dependsRhoSie_.dPdE);
  other.dependsRhoSie_.dTdRho = getOnDeviceDataBox(dependsRhoSie_.dTdRho);
  other.dependsRhoSie_.dTdE   = getOnDeviceDataBox(dependsRhoSie_.dTdE);
  other.dependsRhoSie_.dEdRho = getOnDeviceDataBox(dependsRhoSie_.dEdRho);
  other.numRho_       = numRho_;
  other.lRhoMin_      = lRhoMin_;
  other.lRhoMax_      = lRhoMax_;
  other.rhoMax_       = rhoMax_;
  other.PlRhoMax_     = getOnDeviceDataBox(PlRhoMax_);
  other.dPdRhoMax_    = getOnDeviceDataBox(dPdRhoMax_);
  other.lRhoOffset_   = lRhoOffset_;
  other.lTOffset_     = lTOffset_;
  other.lEOffset_     = lEOffset_;
  other.rhoNormal_    = rhoNormal_;
  other.TNormal_      = TNormal_;
  other.sieNormal_    = sieNormal_;
  other.PNormal_      = PNormal_;
  other.CvNormal_     = CvNormal_;
  other.bModNormal_   = bModNormal_;
  other.dPdENormal_   = dPdENormal_;
  other.dVdTNormal_   = dVdTNormal_;
  other.matid_        = matid_;
  other.reproducible_ = reproducible_;
  other.status_       = status_;
  other.memoryStatus_ = DataStatus::OnDevice;
  return other;
}

void SpinerEOSDependsRhoSie::Finalize() {
  sie_.finalize();
  T_.finalize();
  dependsRhoT_.P.finalize();
  dependsRhoT_.bMod.finalize();
  dependsRhoT_.dPdRho.finalize();
  dependsRhoT_.dPdE.finalize();
  dependsRhoT_.dTdRho.finalize();
  dependsRhoT_.dTdE.finalize();
  dependsRhoT_.dEdRho.finalize();
  dependsRhoSie_.P.finalize();
  dependsRhoSie_.bMod.finalize();
  dependsRhoSie_.dPdRho.finalize();
  dependsRhoSie_.dPdE.finalize();
  dependsRhoSie_.dTdRho.finalize();
  dependsRhoSie_.dTdE.finalize();
  dependsRhoSie_.dEdRho.finalize();
  if (memoryStatus_ == DataStatus::OnDevice) { // these are slices on host
    PlRhoMax_.finalize();
    dPdRhoMax_.finalize();
  }
  memoryStatus_ = DataStatus::Deallocated;
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoSie::TemperatureFromDensityInternalEnergy(const Real rho,
                                                                  const Real sie,
                                                                  Real *lambda) const {
  return interpRhoSie_(rho,sie,T_,lambda);
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoSie::InternalEnergyFromDensityTemperature(const Real rho,
                                                                  const Real T,
                                                                  Real *lambda) const {
  return interpRhoT_(rho,T,sie_,lambda);
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoSie::PressureFromDensityTemperature(const Real rho,
                                                            const Real T,
                                                            Real *lambda) const {
  return interpRhoT_(rho,T,dependsRhoT_.P,lambda);
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoSie::PressureFromDensityInternalEnergy(const Real rho,
                                                               const Real sie,
                                                               Real *lambda) const {
  return interpRhoSie_(rho,sie,dependsRhoSie_.P,lambda);
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoSie::SpecificHeatFromDensityTemperature(const Real rho,
                                                                const Real T,
                                                                Real *lambda) const {
  return 1./interpRhoT_(rho,T,dependsRhoT_.dTdE,lambda);
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoSie::SpecificHeatFromDensityInternalEnergy(const Real rho,
                                                                   const Real sie,
                                                                   Real *lambda) const {
  return 1./interpRhoSie_(rho,sie,dependsRhoSie_.dTdE,lambda);
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoSie::BulkModulusFromDensityTemperature(const Real rho,
                                                               const Real T,
                                                               Real *lambda) const {
  return interpRhoT_(rho,T,dependsRhoT_.bMod,lambda);
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoSie::BulkModulusFromDensityInternalEnergy(const Real rho,
                                                                  const Real sie,
                                                                  Real *lambda) const {
  return interpRhoSie_(rho,sie,dependsRhoSie_.bMod,lambda);
}

PORTABLE_FUNCTION
Real
SpinerEOSDependsRhoSie::GruneisenParamFromDensityTemperature(const Real rho,
                                                             const Real T,
                                                             Real *lambda) 
                                                             const {
  const Real dpde = interpRhoT_(rho, T, dependsRhoT_.dPdE, lambda);
  return dpde/rho;
}

PORTABLE_FUNCTION
Real
SpinerEOSDependsRhoSie::GruneisenParamFromDensityInternalEnergy(const Real rho,
                                                                const Real sie,
                                                                Real *lambda)
                                                                const {
  const Real lRho = toLog_(rho,lRhoOffset_);
  const Real lE = toLog_(sie,lEOffset_);
  const Real dpde = dependsRhoSie_.dPdE.interpToReal(lRho,lE);
  return dpde/rho;
}

PORTABLE_FUNCTION
void SpinerEOSDependsRhoSie::DensityEnergyFromPressureTemperature(const Real press,
                                                                  const Real temp,
                                                                  Real *lambda,
                                                                  Real& rho, Real& sie) const {
  Real lT = toLog_(temp,lTOffset_);
  Real lRho = lRhoFromPlT_(press,lT,lambda);
  rho = fromLog_(lRho,lRhoOffset_);
  sie = sie_.interpToReal(lRho,lT);
}

PORTABLE_FUNCTION
void SpinerEOSDependsRhoSie::PTofRE(const Real rho, const Real sie,
                                    Real * lambda,
                                    Real& press, Real& temp,
                                    Real & dpdr, Real & dpde,
                                    Real & dtdr, Real & dtde) const {
  Real lRho = toLog_(rho,lRhoOffset_);
  Real lE = toLog_(sie,lEOffset_);
  press = dependsRhoSie_.P.interpToReal(lRho,lE);
  temp = T_.interpToReal(lRho,lE);
  dpdr = dependsRhoSie_.dPdRho.interpToReal(lRho,lE);
  dpde = dependsRhoSie_.dPdE.interpToReal(lRho,lE);
  dtdr = dependsRhoSie_.dTdRho.interpToReal(lRho,lE);
  dtde = dependsRhoSie_.dTdE.interpToReal(lRho,lE);
}

PORTABLE_FUNCTION
void SpinerEOSDependsRhoSie::FillEos(Real& rho, Real& temp,
                                     Real& energy, Real& press,
                                     Real& cv, Real& bmod,
                                     const unsigned long output,
                                     Real *lambda) const {
  Real lRho, lT, lE;
  if (output == thermalqs::none) {
    UNDEFINED_ERROR;
  }
  if (output&thermalqs::temperature && output&thermalqs::specific_internal_energy) {
    UNDEFINED_ERROR;
  }
  if (output&thermalqs::density) {
    if (!(output&thermalqs::pressure || output&thermalqs::temperature)) {
      lT = toLog_(temp,lTOffset_);
      lRho = lRhoFromPlT_(press,lT,lambda);
      rho = fromLog_(lRho,lRhoOffset_);
    } else {
      UNDEFINED_ERROR;
    } 
  } else {
    lRho = toLog_(rho,lRhoOffset_);
    if (lambda != nullptr) *lambda = lRho;
  }
  if (output&thermalqs::temperature) {
    lE = toLog_(energy,lEOffset_);
    temp = T_.interpToReal(lRho,lE);
    if (output&thermalqs::pressure) {
      press = dependsRhoSie_.P.interpToReal(lRho,lE);
    }
    if (output&thermalqs::specific_heat) {
      cv = 1./dependsRhoSie_.dTdE.interpToReal(lRho,lE);
    }
    if (output&thermalqs::bulk_modulus) {
      bmod = dependsRhoSie_.bMod.interpToReal(lRho,lE);
    }
  }
  if (output&thermalqs::specific_internal_energy) {
    lT = toLog_(temp,lTOffset_);
    energy = sie_.interpToReal(lRho,lT);
    if (output&thermalqs::pressure) {
      press = dependsRhoT_.P.interpToReal(lRho,lT);
    }
    if (output&thermalqs::specific_heat) {
      cv = 1./dependsRhoT_.dTdE.interpToReal(lRho,lT);
    }
    if (output&thermalqs::bulk_modulus) {
      bmod = dependsRhoT_.bMod.interpToReal(lRho,lT);
    }
  }
}

PORTABLE_FUNCTION
void SpinerEOSDependsRhoSie::ValuesAtReferenceState(Real& rho, Real& temp,
                                                    Real& sie, Real& press,
                                                    Real& cv, Real& bmod,
                                                    Real& dpde, Real& dvdt,
                                                    Real *lambda) const {
  rho = rhoNormal_;
  temp = TNormal_;
  sie = sieNormal_;
  press = PNormal_;
  cv = CvNormal_;
  bmod = bModNormal_;
  dpde = dPdENormal_;
  dvdt = dVdTNormal_;
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoSie::interpRhoT_(const Real rho, const Real T,
                                         const Spiner::DataBox& db,
                                         Real *lambda) const {
  const Real lRho = toLog_(rho,lRhoOffset_);
  const Real lT = toLog_(T,lTOffset_);
  if (lambda != nullptr) { *lambda = lRho; }
  return db.interpToReal(lRho, lT);
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoSie::interpRhoSie_(const Real rho, const Real sie,
                                           const Spiner::DataBox& db,
                                           Real *lambda) const {
  const Real lRho = toLog_(rho,lRhoOffset_);
  const Real lE = toLog_(sie,lEOffset_);
  if (lambda != nullptr) { *lambda = lRho; }
  return db.interpToReal(lRho, lE);
}

PORTABLE_FUNCTION
Real SpinerEOSDependsRhoSie::lRhoFromPlT_(const Real P,
                                          const Real lT,
                                          Real* lambda) const {
  using RootFinding1D::findRoot;

  Real lRho;
  Real dPdRhoMax = dPdRhoMax_.interpToReal(lT);
  Real PMax = PlRhoMax_.interpToReal(lT);
  if ( dPdRhoMax > 0 && P > PMax) {
    Real rho = (P - PMax)/dPdRhoMax + rhoMax_;
    lRho = toLog_(rho,lRhoOffset_);
    counts.increment(0);
  } else {
    Real lRhoGuess = reproducible_ ? lRhoMin_ : 0.5*(lRhoMin_ + lRhoMax_);
    if (lambda != nullptr && lRhoMin_ <= *lambda && *lambda <= lRhoMax_) {
      lRhoGuess = *lambda;
    }
    const callable_interp::l_interp PFunc(dependsRhoT_.P, lT);
    status_ = findRoot(PFunc,
                       P, lRhoGuess,
                       lRhoMin_, lRhoMax_,
                       EPS, EPS,
                       lRho, counts);

    if (status_ != RootFinding1D::Status::SUCCESS) {
#if EPINER_EOS_VERBOSE
      std::stringstream errorMessage;
      errorMessage << "inverting P table for logRho failed\n"
                   << "matid     = " << matid_ << "\n"
                   << "lT        = " << lT << "\n"
                   << "P         = " << P << "\n"
                   << "lRhoGuess = " << lRhoGuess
                   << std::endl;
      EOS_ERROR(errorMessage.str().c_str());
#endif //SPINER_EOS_VERBOSE
      lRho = reproducible_ ? lRhoMin_ : lRhoGuess;
    }
  }
  if (lambda != nullptr) *lambda = lRho;
  return lRho;
}

} // namespace singularity

#endif
