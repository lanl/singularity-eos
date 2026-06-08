//======================================================================
// sesame2spiner tool for converting eospac to spiner
// Author: Jonah Miller (jonahm@lanl.gov)
// © 2021-2026. Triad National Security, LLC. All rights reserved.  This
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
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/sp5.hpp>

#include "generate_files.hpp"
#include "io_eospac.hpp"
#include "parse_cli.hpp"
#include "parser.hpp"

using namespace EospacWrapper;
using singularity::spiner_table_builder::constructRhoBounds;
using singularity::spiner_table_builder::constructTBounds;
using singularity::spiner_table_builder::SpinerTableGridParams;

herr_t saveMaterial(hid_t loc, const SesameMetadata &metadata, const Bounds &lRhoBounds,
                    const Bounds &lTBounds, const Bounds &leBounds,
                    const std::string &name, const bool addSubtables,
                    Verbosity eospacWarn) {

  const int matid = metadata.matid;
  std::string sMatid = std::to_string(matid);

  herr_t status = 0;
  hid_t matGroup, lTGroup, leGroup, coldGroup, mfGroup;

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

  {
    DataBox mf, mask;
    Bounds nphBounds;
    std::string phase_names;

    if (eosMassFraction(matid, lRhoBounds, lTBounds, nphBounds, mf, mask, phase_names,
                        eospacWarn)) {
      mfGroup = H5Gcreate(matGroup, SP5::Depends::massFrac, H5P_DEFAULT, H5P_DEFAULT,
                          H5P_DEFAULT);

      status += mf.saveHDF(mfGroup, SP5::Fields::massFrac);
      const int numphases = mf.dim(1);
      status += H5LTset_attribute_int(mfGroup, ".", "numphases", &numphases, 1);
      status +=
          H5LTset_attribute_string(mfGroup, ".", "phase names", phase_names.c_str());
      status += H5Gclose(mfGroup);
    }
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
    AddMaterials(params, matids, filename);
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
  DataBox P, T, bMod, dPdRho, dPdE, dTdRho, dTdE, dEdRho, mask;
  eosDataOfRhoSie(matid, split, lRhoBounds, leBounds, P, T, bMod, dPdRho, dPdE, dTdRho,
                  dTdE, dEdRho, mask, eospacWarn);
  status += P.saveHDF(loc, SP5::Fields::P);
  status += T.saveHDF(loc, SP5::Fields::T);
  status += bMod.saveHDF(loc, SP5::Fields::bMod);
  status += dPdRho.saveHDF(loc, SP5::Fields::dPdRho);
  status += dPdE.saveHDF(loc, SP5::Fields::dPdE);
  status += dTdRho.saveHDF(loc, SP5::Fields::dTdRho);
  status += dTdE.saveHDF(loc, SP5::Fields::dTdE);
  status += dEdRho.saveHDF(loc, SP5::Fields::dEdRho);
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

// Convert string-based Params to structured SpinerTableGridParams
// This bridges sesame2spiner's parameter file format with the SpinerEOS constructor API
SpinerTableGridParams paramsToGridParams(int matid, const SesameMetadata &metadata,
                                         const Params &params) {
  SpinerTableGridParams gridParams;

  // The "epsilon" shifts here are required to avoid eospac
  // extrapolation errors at table bounds
  constexpr Real TINY = std::numeric_limits<Real>::epsilon();
  auto TinyShift = [=](Real val, int sign) {
    Real shift = std::abs(std::min(10 * val * TINY, TINY));
    return val + sign * shift;
  };

  // Extract bounds from params or metadata
  gridParams.rhoMin = params.Get("rhomin", TinyShift(metadata.rhoMin, 1));
  gridParams.rhoMax = params.Get("rhomax", metadata.rhoMax);
  gridParams.TMin = params.Get("Tmin", TinyShift(metadata.TMin, 1));
  gridParams.TMax = params.Get("Tmax", metadata.TMax);
  gridParams.sieMin = params.Get("siemin", TinyShift(metadata.sieMin, 1));
  gridParams.sieMax = params.Get("siemax", metadata.sieMax);

  // Grid resolution controls
  gridParams.numRhoPerDecade = params.Get("numrho/decade", PPD_DEFAULT_RHO);
  gridParams.numTPerDecade = params.Get("numT/decade", PPD_DEFAULT_T);
  gridParams.numSiePerDecade = params.Get("numSie/decade", PPD_DEFAULT_T);

  // Allow direct specification of grid points (overrides per-decade)
  int numRhoDefault = Bounds::getNumPointsFromPPD(gridParams.rhoMin, gridParams.rhoMax,
                                                  gridParams.numRhoPerDecade);
  int numTDefault = Bounds::getNumPointsFromPPD(gridParams.TMin, gridParams.TMax,
                                                gridParams.numTPerDecade);
  int numSieDefault = Bounds::getNumPointsFromPPD(gridParams.sieMin, gridParams.sieMax,
                                                  gridParams.numSiePerDecade);
  gridParams.numRho = params.Get("numrho", numRhoDefault);
  gridParams.numT = params.Get("numT", numTDefault);
  gridParams.numSie = params.Get("numsie", numSieDefault);

  // Shrink bounds controls
  gridParams.shrinklRhoBounds =
      std::min(1., std::max(params.Get("shrinklRhoBounds", 0.), 0.));
  gridParams.shrinklTBounds =
      std::min(1., std::max(params.Get("shrinklTBounds", 0.), 0.));
  gridParams.shrinkleBounds =
      std::min(1., std::max(params.Get("shrinkleBounds", 0.), 0.));

  // Warnings for inconsistent settings
  if (gridParams.shrinklRhoBounds > 0 &&
      (params.Contains("rhomin") || params.Contains("rhomax"))) {
    std::cerr << "WARNING [" << matid << "]: "
              << "shrinklRhoBounds > 0 and rhomin or rhomax set" << std::endl;
  }
  if (gridParams.shrinklTBounds > 0 &&
      (params.Contains("Tmin") || params.Contains("Tmax"))) {
    std::cerr << "WARNING [" << matid << "]: "
              << "shrinklTBounds > 0 and Tmin or Tmax set" << std::endl;
  }
  if (gridParams.shrinkleBounds > 0 &&
      (params.Contains("siemin") || params.Contains("siemax"))) {
    std::cerr << "WARNING [" << matid << "]: "
              << "shrinkleBounds > 0 and siemin or siemax set" << std::endl;
  }

  // Material properties
  gridParams.matid = matid;
  gridParams.rhoNormal = params.Get("rho_fine_center", metadata.normalDensity);
  if (std::isnan(gridParams.rhoNormal) || gridParams.rhoNormal <= 0 ||
      gridParams.rhoNormal > 1e8) {
    std::cerr << "WARNING [" << matid << "] "
              << "normal density ill defined. Setting it to a sensible default."
              << std::endl;
    gridParams.rhoNormal = 1.0;
  }

  // Piecewise grid controls
  gridParams.piecewiseRho = params.Get("piecewiseRho", true);
  gridParams.piecewiseT = params.Get("piecewiseT", true);
  gridParams.piecewiseSie = params.Get("piecewiseSie", true);

  gridParams.rhoCoarseFactorLo =
      params.Get("rhoCoarseFactorLo", COARSE_FACTOR_DEFAULT_RHO_LO);
  gridParams.rhoCoarseFactorHi =
      params.Get("rhoCoarseFactorHi", COARSE_FACTOR_DEFAULT_RHO_HI);
  gridParams.TCoarseFactor = params.Get("TCoarseFactor", COARSE_FACTOR_DEFAULT_T);
  gridParams.sieCoarseFactor = params.Get("sieCoarseFactor", COARSE_FACTOR_DEFAULT_T);
  gridParams.rhoFineDiameterDecades =
      params.Get("rhoFineDiameterDecades", RHO_FINE_DIAMETER_DEFAULT);
  gridParams.TSplitPoint = params.Get("TSplitPoint", T_SPLIT_POINT_DEFAULT);

  // Optional fine grid bounds override
  Real rho_fine_min = params.Get("rhoFineMin", -1.0);
  Real rho_fine_max = params.Get("rhoFineMax", -1.0);
  if (rho_fine_min * rho_fine_max < 0) {
    std::cerr << "WARNING [" << matid << "]: "
              << "Either rhoFineMin or rhoFineMax is set while the other is still unset."
              << " Both must be set to be sensible. Ignoring." << std::endl;
    rho_fine_min = rho_fine_max = -1.0;
  }
  gridParams.rhoFineMin = rho_fine_min;
  gridParams.rhoFineMax = rho_fine_max;

  // Note: Abar and Zbar are not in Params, they come from metadata in sesame2spiner
  // Offsets are auto-computed by constructRhoBounds/constructTBounds, not in Params

  return gridParams;
}

void getMatBounds(int i, int matid, const SesameMetadata &metadata, const Params &params,
                  Bounds &lRhoBounds, Bounds &lTBounds, Bounds &leBounds) {

  // Convert string-based Params to structured grid parameters
  SpinerTableGridParams gridParams = paramsToGridParams(matid, metadata, params);

  // Validate that requested bounds are within metadata bounds
  checkValInMatBounds(matid, "rhoMin", gridParams.rhoMin, metadata.rhoMin,
                      metadata.rhoMax);
  checkValInMatBounds(matid, "rhoMax", gridParams.rhoMax, metadata.rhoMin,
                      metadata.rhoMax);
  checkValInMatBounds(matid, "TMin", gridParams.TMin, metadata.TMin, metadata.TMax);
  checkValInMatBounds(matid, "TMax", gridParams.TMax, metadata.TMin, metadata.TMax);
  checkValInMatBounds(matid, "sieMin", gridParams.sieMin, metadata.sieMin,
                      metadata.sieMax);
  checkValInMatBounds(matid, "sieMax", gridParams.sieMax, metadata.sieMin,
                      metadata.sieMax);

  // Use shared grid construction utilities for rho and T grids
  constructRhoBounds(gridParams, lRhoBounds);
  constructTBounds(gridParams, lTBounds);

  // Sie grid construction: requires EOSPAC queries for anchor points
  // Cannot use constructSieBounds because it needs a source EOS, not EOSPAC
  if (gridParams.piecewiseSie) {
    // Compute sie anchor points at (rhoNormal, TAnchor) and (rhoNormal, TSplitPoint)
    constexpr Real TAnchor = 298.15; // Room temperature
    constexpr int NT = 1;
    constexpr EOS_INTEGER nXYPairs = 2;
    EOS_INTEGER tableHandle[NT];
    EOS_INTEGER tableType[NT] = {EOS_Ut_DT};
    EOS_REAL rho[2], T[2], sie[2], dx[2], dy[2];
    {
      eosSafeLoad(NT, matid, tableType, tableHandle, {"EOS_Ut_DT"}, Verbosity::Quiet);
      EOS_INTEGER eospacEofRT = tableHandle[0];
      rho[0] = rho[1] = densityToSesame(gridParams.rhoNormal);
      T[0] = temperatureToSesame(TAnchor);
      T[1] = temperatureToSesame(gridParams.TSplitPoint);
      eosSafeInterpolate(&eospacEofRT, nXYPairs, rho, T, sie, dx, dy, "EofRT",
                         Verbosity::Quiet);
      eosSafeDestroy(NT, tableHandle, Verbosity::Quiet);
    }
    const Real sieAnchor = sie[0];
    const Real sieSplitPoint = sie[1];
    leBounds = Bounds(Bounds::TwoGrids(), gridParams.sieMin, gridParams.sieMax, sieAnchor,
                      sieSplitPoint, gridParams.numSiePerDecade,
                      gridParams.sieCoarseFactor, true, gridParams.shrinkleBounds);
  } else {
    // Uniform sie grid
    int numSie = gridParams.numSie;
    if (numSie <= 0) {
      numSie = Bounds::getNumPointsFromPPD(gridParams.sieMin, gridParams.sieMax,
                                           gridParams.numSiePerDecade);
    }
    leBounds = Bounds(gridParams.sieMin, gridParams.sieMax, numSie, true,
                      gridParams.shrinkleBounds);
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
