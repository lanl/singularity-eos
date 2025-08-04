//======================================================================
// sesame2spiner tool for converting eospac to spiner
// Author: Jonah Miller (jonahm@lanl.gov)
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
//======================================================================

#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <hdf5.h>
#include <hdf5_hl.h>

#ifndef SPINER_USE_HDF
#error "HDF5 must be enabled"
#endif // SPINER_USE_HDF

#include <eospac-wrapper/eospac_wrapper.hpp>
#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/fast-math/logs.hpp>
#include <singularity-eos/base/sp5/singularity_eos_sp5.hpp>
#include <singularity-eos/eos/eos_spiner_sie_transforms.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/sp5.hpp>

#include "generate_files.hpp"
#include "io_eospac.hpp"
#include "parse_cli.hpp"
#include "parser.hpp"

using namespace EospacWrapper;

herr_t saveMaterial(hid_t loc, const SesameMetadata &metadata, const Bounds &lRhoBounds,
                    const Bounds &lTBounds, const Bounds &leBounds,
                    const std::string &name, const bool addSubtables,
                    Verbosity eospacWarn) {

  const int matid = metadata.matid;
  std::string sMatid = std::to_string(matid);

  herr_t status = 0;
  hid_t matGroup, lTGroup, leGroup, coldGroup;

  matGroup = H5Gcreate(loc, sMatid.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status += H5Lcreate_soft(sMatid.c_str(), loc, name.c_str(), H5P_DEFAULT, H5P_DEFAULT);

  // Dependent variables metadata
  status += H5LTset_attribute_string(loc, sMatid.c_str(), SP5::Offsets::messageName,
                                     SP5::Offsets::message);
  status += H5LTset_attribute_double(loc, sMatid.c_str(), SP5::Offsets::rho,
                                     &lRhoBounds.offset, 1);
  status +=
      H5LTset_attribute_double(loc, sMatid.c_str(), SP5::Offsets::T, &lTBounds.offset, 1);
  status += H5LTset_attribute_double(loc, sMatid.c_str(), SP5::Offsets::sie,
                                     &leBounds.offset, 1);

  // Material metadata
  status +=
      H5LTset_attribute_double(loc, sMatid.c_str(), SP5::Material::exchangeCoefficient,
                               &metadata.exchangeCoefficient, 1);
  status += H5LTset_attribute_double(loc, sMatid.c_str(), SP5::Material::meanAtomicMass,
                                     &metadata.meanAtomicMass, 1);
  status += H5LTset_attribute_double(loc, sMatid.c_str(), SP5::Material::meanAtomicNumber,
                                     &metadata.meanAtomicNumber, 1);
  status += H5LTset_attribute_double(loc, sMatid.c_str(), SP5::Material::solidBulkModulus,
                                     &metadata.solidBulkModulus, 1);
  status += H5LTset_attribute_double(loc, sMatid.c_str(), SP5::Material::normalDensity,
                                     &metadata.normalDensity, 1);
  status += H5LTset_attribute_string(loc, sMatid.c_str(), SP5::Material::comments,
                                     metadata.comments.c_str());
  status += H5LTset_attribute_int(loc, sMatid.c_str(), SP5::Material::matid,
                                  &metadata.matid, 1);
  status += H5LTset_attribute_string(loc, sMatid.c_str(), SP5::Material::name,
                                     metadata.name.c_str());

  lTGroup = H5Gcreate(matGroup, SP5::Depends::logRhoLogT, H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);

  leGroup = H5Gcreate(matGroup, SP5::Depends::logRhoLogSie, H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
  coldGroup =
      H5Gcreate(matGroup, SP5::Depends::coldCurve, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  status += saveTablesRhoSie(leGroup, matid, TableSplit::Total, lRhoBounds, leBounds,
                             eospacWarn);
  status +=
      saveTablesRhoT(lTGroup, matid, TableSplit::Total, lRhoBounds, lTBounds, eospacWarn);
  {
    DataBox P, sie, dPdRho, dEdRho, bMod, mask, transitionMask;
    eosColdCurves(matid, lRhoBounds, P, sie, dPdRho, dEdRho, bMod, mask, eospacWarn);
    // currently unused
    // eosColdCurveMask(matid, lRhoBounds, leBounds.grid.nPoints(), sie, transitionMask,
    //                  eospacWarn);

    status += P.saveHDF(coldGroup, SP5::Fields::P);
    status += sie.saveHDF(coldGroup, SP5::Fields::sie);
    status += bMod.saveHDF(coldGroup, SP5::Fields::bMod);
    status += dPdRho.saveHDF(coldGroup, SP5::Fields::dPdRho);
    status += dEdRho.saveHDF(coldGroup, SP5::Fields::dEdRho);
    // currently unused
    // status += mask.saveHDF(coldGroup, SP5::Fields::mask);
    // status += transitionMask.saveHDF(coldGroup, SP5::Fields::transitionMask);
  }

  if (addSubtables) {
    int i = 0;
    std::vector<TableSplit> splits = {TableSplit::ElectronOnly, TableSplit::IonCold};
    std::vector<std::string> grpnames = {SP5::SubTable::electronOnly,
                                         SP5::SubTable::ionCold};
    for (auto split : splits) {
      std::string grpname = grpnames[i++];
      {
        hid_t grp =
            H5Gcreate(leGroup, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status += saveTablesRhoSie(grp, matid, split, lRhoBounds, leBounds, eospacWarn);
        status += H5Gclose(grp);
      }
      {
        hid_t grp =
            H5Gcreate(lTGroup, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status += saveTablesRhoT(grp, matid, split, lRhoBounds, lTBounds, eospacWarn);
        status += H5Gclose(grp);
      }
    }
  }

  status += H5Gclose(leGroup);
  status += H5Gclose(lTGroup);
  status += H5Gclose(coldGroup);
  status += H5Gclose(matGroup);

  return status;
}

herr_t saveAllMaterials(const std::string &savename,
                        const std::vector<std::string> &filenames, bool printMetadata,
                        Verbosity eospacWarn) {
  std::vector<Params> params;
  std::vector<int> matids;
  std::unordered_map<std::string, int> used_names;
  std::unordered_set<int> used_matids;
  SesameMetadata metadata;
  hid_t file;
  herr_t status = H5_SUCCESS;

  for (auto const &filename : filenames) {
    Params p(filename);
    if (!p.Contains("matid")) {
      std::cerr << "Material file " << filename << "is missing matid.\n"
                << "Example input files:\n"
                << EXAMPLESTRING << std::endl;
      std::exit(1);
    }
    matids.push_back(p.Get<int>("matid"));
    params.push_back(p);
  }

  std::cout << "Saving to file " << savename << std::endl;
  file = H5Fcreate(savename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // singularity version
  H5LTset_attribute_string(file, "/", "singularity_version", SINGULARITY_VERSION);
  // log type. 0 for true, 1 for NQT1, 2 for NQT2, -1 for single precision true
  int log_type = singularity::FastMath::Settings::log_type;
  H5LTset_attribute_int(file, "/", SP5::logType, &log_type, 1);

  std::cout << "Processing " << matids.size() << " materials..." << std::endl;

  for (size_t i = 0; i < matids.size(); i++) {
    int matid = matids[i];
    if (used_matids.count(matid) > 0) {
      std::cerr << "...Duplicate matid " << matid << " detected. Skipping." << std::endl;
      continue;
    }
    used_matids.insert(matid);

    std::cout << "..." << matid << std::endl;

    eosGetMetadata(matid, metadata, Verbosity::Debug);
    if (printMetadata) std::cout << metadata << std::endl;

    std::string name = params[i].Get("name", metadata.name);
    if (name == "-1" || name == "") {
      std::string new_name = "material_" + std::to_string(i);
      std::cerr << "...WARNING: no reasonable name found. "
                << "Using a default name: " << new_name << std::endl;
      name = new_name;
    }
    if (used_names.count(name) > 0) {
      used_names[name] += 1;
      std::string new_name = name + "_" + std::to_string(used_names[name]);
      std::cerr << "...WARNING: Name " << name << " already used. "
                << "Using name: " << new_name << std::endl;
      name = new_name;
    } else {
      used_names[name] = 1;
    }

    Bounds lRhoBounds, lTBounds, leBounds;
    getMatBounds(i, matid, metadata, params[i], lRhoBounds, lTBounds, leBounds);

    if (eospacWarn == Verbosity::Debug) {
      std::cout << "bounds for log(rho), log(T), log(sie) are:\n"
                << lRhoBounds << lTBounds << leBounds << std::endl;
    }

    const bool add_subtables = params[i].Get("ionization", false);
    if (eospacWarn == Verbosity::Debug) {
      std::cout << "Adding subtables for partial ionization? " << add_subtables
                << std::endl;
    }

    status += saveMaterial(file, metadata, lRhoBounds, lTBounds, leBounds, name,
                           add_subtables, eospacWarn);
    if (status != H5_SUCCESS) {
      std::cerr << "WARNING: problem with HDf5" << std::endl;
    }
  }

  std::cout << "Cleaning up." << std::endl;
  status += H5Fclose(file);
  if (status != H5_SUCCESS) {
    std::cerr << "WARNING: problem with HDf5" << std::endl;
  }
  return status;
}

herr_t saveTablesRhoSie(hid_t loc, int matid, TableSplit split, const Bounds &lRhoBounds,
                        const Bounds &leBounds, Verbosity eospacWarn) {
  herr_t status = 0;
  DataBox P, T, bMod, dPdRho, dPdE, dTdRho, dTdE, dEdRho, sie_shift, mask;
  eosDataOfRhoSie(matid, split, lRhoBounds, leBounds, P, T, bMod, dPdRho, dPdE, dTdRho,
                  dTdE, dEdRho, sie_shift, mask, eospacWarn);
  status += P.saveHDF(loc, SP5::Fields::P);
  status += T.saveHDF(loc, SP5::Fields::T);
  status += bMod.saveHDF(loc, SP5::Fields::bMod);
  status += dPdRho.saveHDF(loc, SP5::Fields::dPdRho);
  status += dPdE.saveHDF(loc, SP5::Fields::dPdE);
  status += dTdRho.saveHDF(loc, SP5::Fields::dTdRho);
  status += dTdE.saveHDF(loc, SP5::Fields::dTdE);
  status += dEdRho.saveHDF(loc, SP5::Fields::dEdRho);
  //status += sie_shift.saveHDF(loc, SP5::Fields::sie);
  // currently unused
  // status += mask.saveHDF(loc, SP5::Fields::mask);
  return status;
}

herr_t saveTablesRhoT(hid_t loc, int matid, TableSplit split, const Bounds &lRhoBounds,
                      const Bounds &lTBounds, Verbosity eospacWarn) {
  herr_t status = 0;
  DataBox P, sie, bMod, dPdRho, dPdE, dTdRho, dTdE, dEdRho, dEdT, mask;
  eosDataOfRhoT(matid, split, lRhoBounds, lTBounds, P, sie, bMod, dPdRho, dPdE, dTdRho,
                dTdE, dEdRho, dEdT, mask, eospacWarn);
  status += P.saveHDF(loc, SP5::Fields::P);
  status += sie.saveHDF(loc, SP5::Fields::sie);
  status += bMod.saveHDF(loc, SP5::Fields::bMod);
  status += dPdRho.saveHDF(loc, SP5::Fields::dPdRho);
  status += dPdE.saveHDF(loc, SP5::Fields::dPdE);
  status += dTdRho.saveHDF(loc, SP5::Fields::dTdRho);
  status += dTdE.saveHDF(loc, SP5::Fields::dTdE);
  status += dEdRho.saveHDF(loc, SP5::Fields::dEdRho);
  status += dEdT.saveHDF(loc, SP5::Fields::dEdT);
  // Currently unused
  // status += mask.saveHDF(loc, SP5::Fields::mask);
  return status;
}

void getMatBounds(int i, int matid, const SesameMetadata &metadata, const Params &params,
                  Bounds &lRhoBounds, Bounds &lTBounds, Bounds &leBounds) {
 
  //move to top of function so all if and else statements have access
  using namespace singularity;
  TransformDataContainer data(matid, Verbosity::Quiet);
  AllTransform<TransformDataContainer> shift(data);

  // The "epsilon" shifts here are required to avoid eospac
  // extrapolation errors at table bounds
  constexpr Real TINY = std::numeric_limits<Real>::epsilon();
  auto TinyShift = [=](Real val, int sign) {
    Real shift = std::abs(std::min(10 * val * TINY, TINY));
    return val + sign * shift;
  };

  Real rhoMin = params.Get("rhomin", TinyShift(metadata.rhoMin, 1));
  Real rhoMax = params.Get("rhomax", metadata.rhoMax);
  Real TMin = params.Get("Tmin", TinyShift(metadata.TMin, 1));
  Real TMax = params.Get("Tmax", metadata.TMax);
  Real sieMin = params.Get("siemin", TinyShift(metadata.sieMin, 1));
  Real sieMax = params.Get("siemax", metadata.sieMax);

  checkValInMatBounds(matid, "rhoMin", rhoMin, metadata.rhoMin, metadata.rhoMax);
  checkValInMatBounds(matid, "rhoMax", rhoMax, metadata.rhoMin, metadata.rhoMax);
  checkValInMatBounds(matid, "TMin", TMin, metadata.TMin, metadata.TMax);
  checkValInMatBounds(matid, "TMax", TMax, metadata.TMin, metadata.TMax);
  checkValInMatBounds(matid, "sieMin", sieMin, metadata.sieMin, metadata.sieMax);
  checkValInMatBounds(matid, "sieMax", sieMax, metadata.sieMin, metadata.sieMax);

  Real shrinklRhoBounds = params.Get("shrinklRhoBounds", 0.);
  Real shrinklTBounds = params.Get("shrinklTBounds", 0.);
  Real shrinkleBounds = params.Get("shrinkleBounds", 0.);

  shrinklRhoBounds = std::min(1., std::max(shrinklRhoBounds, 0.));
  shrinklTBounds = std::min(1., std::max(shrinklTBounds, 0.));
  shrinkleBounds = std::min(1., std::max(shrinkleBounds, 0.));

  if (shrinklRhoBounds > 0 && (params.Contains("rhomin") || params.Contains("rhomax"))) {
    std::cerr << "WARNING [" << matid << "]: "
              << "shrinklRhoBounds > 0 and rhomin or rhomax set" << std::endl;
  }
  if (shrinklTBounds > 0 && (params.Contains("Tmin") || params.Contains("Tmax"))) {
    std::cerr << "WARNING [" << matid << "]: "
              << "shrinklTBounds > 0 and Tmin or Tmax set" << std::endl;
  }
  if (shrinkleBounds > 0 && (params.Contains("siemin") || params.Contains("siemax"))) {
    std::cerr << "WARNING [" << matid << "]: "
              << "shrinkleBounds > 0 and siemin or siemax set" << std::endl;
  }

  int ppdRho = params.Get("numrho/decade", PPD_DEFAULT_RHO);
  int numRhoDefault = Bounds::getNumPointsFromPPD(rhoMin, rhoMax, ppdRho);

  int ppdT = params.Get("numT/decade", PPD_DEFAULT_T);
  int numTDefault = Bounds::getNumPointsFromPPD(TMin, TMax, ppdT);

  int ppdSie = params.Get("numSie/decade", PPD_DEFAULT_T);
  int numSieDefault = Bounds::getNumPointsFromPPD(sieMin, sieMax, ppdSie);

  int numRho = params.Get("numrho", numRhoDefault);
  int numT = params.Get("numT", numTDefault);
  int numSie = params.Get("numsie", numSieDefault);

  constexpr Real TAnchor = 298.15;
  Real rhoAnchor = params.Get("rho_fine_center", metadata.normalDensity);
  if (std::isnan(rhoAnchor) || rhoAnchor <= 0 || rhoAnchor > 1e8) {
    std::cerr << "WARNING [" << matid << "] "
              << "normal density ill defined. Setting it to a sensible default."
              << std::endl;
    rhoAnchor = 1;
  }

  // Piecewise grids stuff
  const bool piecewiseRho = params.Get("piecewiseRho", true);
  const bool piecewiseT = params.Get("piecewiseT", true);
  const bool piecewiseSie = params.Get("piecewiseSie", true);

  const Real ppd_factor_rho_lo =
      params.Get("rhoCoarseFactorLo", COARSE_FACTOR_DEFAULT_RHO_LO);
  const Real ppd_factor_rho_hi =
      params.Get("rhoCoarseFactorHi", COARSE_FACTOR_DEFAULT_RHO_HI);
  const Real ppd_factor_T = params.Get("TCoarseFactor", COARSE_FACTOR_DEFAULT_T);
  const Real ppd_factor_sie = params.Get("sieCoarseFactor", COARSE_FACTOR_DEFAULT_T);
  const Real rho_fine_diameter =
      params.Get("rhoFineDiameterDecades", RHO_FINE_DIAMETER_DEFAULT);
  const Real TSplitPoint = params.Get("TSplitPoint", T_SPLIT_POINT_DEFAULT);
  const Real rho_fine_center = rhoAnchor;

  // These override the rho center/diameter settings
  Real rho_fine_min = params.Get("rhoFineMin", -1);
  Real rho_fine_max = params.Get("rhoFineMax", -1);
  if (rho_fine_min * rho_fine_max < 0) {
    std::cerr << "WARNING [" << matid << "]: "
              << "Either rhoFineMin or rhoFineMax is set while the other is still unset."
              << " Both must be set to be sensible. Ignoring." << std::endl;
    rho_fine_min = rho_fine_max = -1;
  }

  // Forces density and temperature to be in a region where an offset
  // is not needed. This improves resolution at low densities and
  // temperatures.

  // Extrapolation and other resolution tricks will be explored in the
  // future.
  if (rhoMin < STRICTLY_POS_MIN_RHO) rhoMin = STRICTLY_POS_MIN_RHO;
  if (TMin < STRICTLY_POS_MIN_T) TMin = STRICTLY_POS_MIN_T;

  if (piecewiseRho) {
    if (rho_fine_min > 0) {
      lRhoBounds = Bounds(Bounds::ThreeGrids(), rhoMin, rhoMax, rho_fine_center,
                          rho_fine_min, rho_fine_max, ppdRho, ppd_factor_rho_lo,
                          ppd_factor_rho_hi, true, shrinklRhoBounds);
    } else {
      lRhoBounds =
          Bounds(Bounds::ThreeGrids(), rhoMin, rhoMax, rho_fine_center, rho_fine_diameter,
                 ppdRho, ppd_factor_rho_lo, ppd_factor_rho_hi, true, shrinklRhoBounds);
    }
  } else {
    lRhoBounds = Bounds(rhoMin, rhoMax, numRho, true, shrinklRhoBounds, rhoAnchor);
  }
  if (piecewiseT) {
    lTBounds = Bounds(Bounds::TwoGrids(), TMin, TMax, TAnchor, TSplitPoint, ppdT,
                      ppd_factor_T, true, shrinklTBounds);
  } else {
    lTBounds = Bounds(TMin, TMax, numT, true, shrinklTBounds, TAnchor);
  }
  if (piecewiseSie) {
    // compute temperature as a reasonable anchor point
    constexpr int NT = 1;
    constexpr EOS_INTEGER nXYPairs = 2;
    EOS_INTEGER tableHandle[NT];
    EOS_INTEGER tableType[NT] = {EOS_Ut_DT};
    EOS_REAL rho[2], T[2], sie[2], dx[2], dy[2], sie_transformed[2];
    {
      eosSafeLoad(NT, matid, tableType, tableHandle, {"EOS_Ut_DT"}, Verbosity::Quiet);
      EOS_INTEGER eospacEofRT = tableHandle[0];
      rho[0] = rho[1] = densityToSesame(rhoAnchor);
      T[0] = temperatureToSesame(TAnchor);
      T[1] = temperatureToSesame(TSplitPoint);
      eosSafeInterpolate(&eospacEofRT, nXYPairs, rho, T, sie, dx, dy, "EofRT",
                         Verbosity::Quiet);

      //std::vector<Real> sie_transformed(nXYPairs);
      for (int i = 0; i < nXYPairs; ++i) {
        Real sie_phys = sieToSesame(sie[i]);
        Real rho_phys = densityToSesame(rho[i]);
        sie_transformed[i] = shift.transform(sie_phys, rho_phys);
    }
     
      eosSafeDestroy(NT, tableHandle, Verbosity::Quiet);
    }

    const Real sieAnchor = sie_transformed[0];
    const Real sieSplitPoint = sie_transformed[1]; //sie_tranformed when testing transformations
    leBounds = Bounds(Bounds::TwoGrids(), sieMin, sieMax, sieAnchor, sieSplitPoint,
                      ppdSie, ppd_factor_sie, true, shrinkleBounds);
  } else {


     //shift these as well for bounds creation in following else statemnt?
     sieMin = shift.transform(sieMin, rhoMIn);
     sieMax = shift.transform(sieMax, rhoMax);
     numsie = shift.transform(numSie, numRho);
     


    leBounds = Bounds(sieMin, sieMax, numSie, true, shrinkleBounds);
  }

  std::cout << "lRho bounds are\n"
            << lRhoBounds << "lT bounds are\n"
            << lTBounds << "lSie bounds are \n"
            << leBounds << std::endl;

  return;
}

bool checkValInMatBounds(int matid, const std::string &name, Real val, Real vmin,
                         Real vmax) {
  if (val < vmin || val > vmax) {
    std::cerr << "WARNING [" << matid << "]: " << name
              << " out of sesame table bounds. Consider changing this.\n"
              << "\t" << name << ", [bounds] = " << val << ", [" << vmin << ", " << vmax
              << "]" << std::endl;
    return false;
  }
  return true;
}
