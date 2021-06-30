//======================================================================
// sesame2spiner tool for converting eospac to spiner
// Author: Jonah Miller (jonahm@lanl.gov)
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
//======================================================================

#include <string>
#include <cmath>
#include <cstdlib>

#include <hdf5.h>
#include <hdf5_hl.h>
#include <nlohmann/json.hpp>

#ifndef SPINER_USE_HDF
#error "HDF5 must be enabled"
#endif // SPINER_USE_HDF

#include <spiner/ports-of-call/portability.hpp>
#include <sp5/singularity_eos_sp5.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/sp5.hpp>

#include "io_eospac.hpp"
#include "generate_files.hpp"
#include "parser.hpp"

herr_t saveMaterial(hid_t loc,
		    const SesameMetadata& metadata,
		    const Bounds& lRhoBounds,
		    const Bounds& lTBounds,
		    const Bounds& leBounds,
		    const std::string& name,
		    Verbosity eospacWarn) {

  const int matid = metadata.matid;
  std::string sMatid = std::to_string(matid);

  herr_t status = 0;
  hid_t matGroup, lTGroup, leGroup, coldGroup;

  matGroup = H5Gcreate(loc, sMatid.c_str(),
		       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status += H5Lcreate_soft(sMatid.c_str(),
			   loc,
			   name.c_str(),
			   H5P_DEFAULT, H5P_DEFAULT);
  
  // Dependent variables metadata
  status += H5LTset_attribute_string(loc, sMatid.c_str(),
				     SP5::Offsets::messageName,
				     SP5::Offsets::message);
  status += H5LTset_attribute_double(loc, sMatid.c_str(),
				     SP5::Offsets::rho,
				     &lRhoBounds.offset,1);
  status += H5LTset_attribute_double(loc, sMatid.c_str(),
				     SP5::Offsets::T,
				     &lTBounds.offset,1);
  status += H5LTset_attribute_double(loc, sMatid.c_str(),
				     SP5::Offsets::sie,
				     &leBounds.offset,1);

  // Material metadata
  status += H5LTset_attribute_double(loc, sMatid.c_str(),
				     SP5::Material::exchangeCoefficient,
				     &metadata.exchangeCoefficient, 1);
  status += H5LTset_attribute_double(loc, sMatid.c_str(),
				     SP5::Material::meanAtomicMass,
				     &metadata.meanAtomicMass, 1);
  status += H5LTset_attribute_double(loc, sMatid.c_str(),
				     SP5::Material::meanAtomicNumber,
				     &metadata.meanAtomicNumber, 1);
  status += H5LTset_attribute_double(loc, sMatid.c_str(),
				     SP5::Material::solidBulkModulus,
				     &metadata.solidBulkModulus, 1);
  status += H5LTset_attribute_double(loc, sMatid.c_str(),
				     SP5::Material::normalDensity,
				     &metadata.normalDensity, 1);
  status += H5LTset_attribute_string(loc, sMatid.c_str(),
				     SP5::Material::comments,
				     metadata.comments.c_str());
  status += H5LTset_attribute_int(loc, sMatid.c_str(),
				  SP5::Material::matid,
				  &metadata.matid, 1);
  status += H5LTset_attribute_string(loc, sMatid.c_str(),
				     SP5::Material::name,
				     metadata.name.c_str());

  lTGroup = H5Gcreate(matGroup, SP5::Depends::logRhoLogT,
		      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  leGroup = H5Gcreate(matGroup, SP5::Depends::logRhoLogSie,
		      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  coldGroup = H5Gcreate(matGroup, SP5::Depends::coldCurve,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  {
    DataBox P, sie, T, bMod, dPdRho, dPdE, dTdRho, dTdE, dEdRho, mask;
    eosDataOfRhoSie(matid, lRhoBounds, leBounds,
                    P, T, bMod, dPdRho, dPdE, dTdRho, dTdE, dEdRho,
                    mask,
                    eospacWarn);
    status += P.saveHDF(leGroup,      SP5::Fields::P);
    status += T.saveHDF(leGroup,      SP5::Fields::T);
    status += bMod.saveHDF(leGroup,   SP5::Fields::bMod);
    status += dPdRho.saveHDF(leGroup, SP5::Fields::dPdRho);
    status += dPdE.saveHDF(leGroup,   SP5::Fields::dPdE);
    status += dTdRho.saveHDF(leGroup, SP5::Fields::dTdRho);
    status += dTdE.saveHDF(leGroup,   SP5::Fields::dTdE);
    status += dEdRho.saveHDF(leGroup, SP5::Fields::dEdRho);
    status += mask.saveHDF(leGroup,   SP5::Fields::mask);
  }

  {
    DataBox P, sie, T, bMod, dPdRho, dPdE, dTdRho, dTdE, dEdRho, dEdT, mask;
    eosDataOfRhoT(matid, lRhoBounds, lTBounds,
                  P, sie, bMod, dPdRho, dPdE, dTdRho, dTdE, dEdRho, dEdT,
                  mask,
                  eospacWarn);
    status += P.saveHDF(lTGroup,      SP5::Fields::P);
    status += sie.saveHDF(lTGroup,    SP5::Fields::sie);
    status += bMod.saveHDF(lTGroup,   SP5::Fields::bMod);
    status += dPdRho.saveHDF(lTGroup, SP5::Fields::dPdRho);
    status += dPdE.saveHDF(lTGroup,   SP5::Fields::dPdE);
    status += dTdRho.saveHDF(lTGroup, SP5::Fields::dTdRho);
    status += dTdE.saveHDF(lTGroup,   SP5::Fields::dTdE);
    status += dEdRho.saveHDF(lTGroup, SP5::Fields::dEdRho);
    status += dEdT.saveHDF(lTGroup,   SP5::Fields::dEdT);
    status += mask.saveHDF(lTGroup,   SP5::Fields::mask);
  }
  {
    DataBox P, sie, dPdRho, dEdRho, bMod, mask, transitionMask;
    eosColdCurves(matid, lRhoBounds,
                  P, sie, dPdRho, dEdRho, bMod, mask,
                  eospacWarn);
    eosColdCurveMask(matid, lRhoBounds, leBounds.grid.nPoints(),
                     sie, transitionMask,
                     eospacWarn);

    status += P.saveHDF(coldGroup,              SP5::Fields::P);
    status += sie.saveHDF(coldGroup,            SP5::Fields::sie);
    status += bMod.saveHDF(coldGroup,           SP5::Fields::bMod);
    status += dPdRho.saveHDF(coldGroup,         SP5::Fields::dPdRho);
    status += dEdRho.saveHDF(coldGroup,         SP5::Fields::dEdRho);
    status += mask.saveHDF(coldGroup,           SP5::Fields::mask);
    status += transitionMask.saveHDF(coldGroup, SP5::Fields::transitionMask);
  }

  status += H5Gclose(leGroup);
  status += H5Gclose(lTGroup);
  status += H5Gclose(coldGroup);
  status += H5Gclose(matGroup);

  return status;
}

herr_t saveAllMaterials(json& params,
			bool printMetadata,
			Verbosity eospacWarn) {
  std::vector<int> matids;
  SesameMetadata metadata;
  hid_t file;
  herr_t status = H5_SUCCESS;

  std::string savename = params.value("savename", SP5::defaultSesFileName);
  
  if ( params["materials"].is_null() ) {
    std::cerr << "Input file must have materials field.\n"
	      << "Example input file:\n"
	      << EXAMPLESTRING
	      << std::endl;
    std::exit(1);
  }
  
  for (auto & matParams : params["materials"]) {
    matids.push_back(matParams["matid"]);
  }
  if (matids.size() < 1) {
    std::cerr << "No materials in input file.\n"
	      << "Example input file:\n"
	      << EXAMPLESTRING
	      << std::endl;
    std::exit(1);
  }
  
  std::cout << "Saving to file " << savename << std::endl;
  file = H5Fcreate(savename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  
  std::cout << "Processing " << matids.size() << " materials..." << std::endl;
  
  for (size_t i = 0; i < matids.size(); i++) {
    int matid = matids[i];
    
    std::cout << "..." << matid << std::endl;
    
    eosGetMetadata(matid, metadata, Verbosity::Debug);
    if (printMetadata) std::cout << metadata << std::endl;
    std::string name = params["materials"][i].value("name",metadata.name); 

    Bounds lRhoBounds, lTBounds, leBounds;
    getMatBounds(i,matid,metadata,params,lRhoBounds,lTBounds,leBounds);

    if (eospacWarn == Verbosity::Debug) {
      std::cout << "bounds for log(rho), log(T), log(sie) are:\n"
		<< lRhoBounds
		<< lTBounds
		<< leBounds
		<< std::endl;
    }
    
    status += saveMaterial(file, metadata,
			   lRhoBounds, lTBounds, leBounds,
			   name, eospacWarn);
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

void getMatBounds(int i,
		  int matid,
		  const SesameMetadata& metadata,
		  json& params,
		  Bounds& lRhoBounds,
		  Bounds& lTBounds,
		  Bounds& leBounds) {
  Real rhoMin   = params["materials"][i].value("rhomin",metadata.rhoMin);
  Real rhoMax   = params["materials"][i].value("rhomax",metadata.rhoMax);
  Real TMin     = params["materials"][i].value("Tmin",metadata.TMin);
  Real TMax     = params["materials"][i].value("Tmax",metadata.TMax);
  Real sieMin   = params["materials"][i].value("siemin",metadata.sieMin);
  Real sieMax   = params["materials"][i].value("siemax",metadata.sieMax);

  checkValInMatBounds(matid, "rhoMin", rhoMin, metadata.rhoMin, metadata.rhoMax);
  checkValInMatBounds(matid, "rhoMax", rhoMax, metadata.rhoMin, metadata.rhoMax);
  checkValInMatBounds(matid, "TMin",   TMin,   metadata.TMin,   metadata.TMax);
  checkValInMatBounds(matid, "TMax",   TMax,   metadata.TMin,   metadata.TMax);
  checkValInMatBounds(matid, "sieMin", sieMin, metadata.sieMin, metadata.sieMax);
  checkValInMatBounds(matid, "sieMax", sieMax, metadata.sieMin, metadata.sieMax);

  Real shrinklRhoBounds = params["materials"][i].value("shrinklRhoBounds", 0.0);
  Real shrinklTBounds   = params["materials"][i].value("shrinklTBounds",   0.0);
  Real shrinkleBounds   = params["materials"][i].value("shrinkleBounds",   0.0);

  shrinklRhoBounds = std::min(1., std::max(shrinklRhoBounds, 0.));
  shrinklTBounds   = std::min(1., std::max(shrinklTBounds,   0.));
  shrinkleBounds   = std::min(1., std::max(shrinkleBounds,   0.));
  
  if (shrinklRhoBounds > 0
      && !(params["materials"][i]["rhomin"].is_null()
	   && params["materials"][i]["rhomax"].is_null())) {
    std::cerr << "WARNING [" << matid << "]: "
	      << "shrinklRhoBounds > 0 and rhomin or rhomax set"
	      << std::endl;
  }
  if (shrinklTBounds > 0
	&& !(params["materials"][i]["Tmin"].is_null()
	     && params["materials"][i]["Tmax"].is_null())) {
    std::cerr << "WARNING [" << matid << "]: "
	      << "shrinklTBounds > 0 and Tmin or Tmax set"
	      << std::endl;
  }
  if (shrinkleBounds > 0
      && !(params["materials"][i]["siemin"].is_null()
	   && params["materials"][i]["siemax"].is_null())) {
    std::cerr << "WARNING [" << matid << "]: "
	      << "shrinkleBounds > 0 and siemin or siemax set"
		<< std::endl;
  }

  int ppdRho = params["materials"][i].value("numrho/decade", PPD_DEFAULT);
  int numRhoDefault = getNumPointsFromPPD(rhoMin,rhoMax,ppdRho);

  int ppdT = params["materials"][i].value("numT/decade",   PPD_DEFAULT);
  int numTDefault = getNumPointsFromPPD(TMin,TMax,ppdT);

  int ppdSie = params["materials"][i].value("numSie/decade", PPD_DEFAULT);
  int numSieDefault = getNumPointsFromPPD(sieMin,sieMax,ppdSie);

  int numRho = params["materials"][i].value("numrho", numRhoDefault);
  int numT   = params["materials"][i].value("numT",     numTDefault);
  int numSie = params["materials"][i].value("numsie", numSieDefault);

  Real rhoAnchor = metadata.normalDensity;
  Real TAnchor = 298.15;

  lRhoBounds = Bounds(rhoMin,rhoMax,numRho,true,shrinklRhoBounds,rhoAnchor);
  lTBounds   = Bounds(TMin,TMax,numT,true,shrinklTBounds,TAnchor);
  leBounds   = Bounds(sieMin,sieMax,numSie,true,shrinkleBounds);
  
  return;
}

bool checkValInMatBounds(int matid,
                         const std::string& name,
                         Real val, Real vmin, Real vmax) {
  if ( val < vmin || val > vmax ) {
    std::cerr << "WARNING [" << matid << "]: "
              << name << " out of sesame table bounds. Consider changing this.\n"
              << "\t" << name << ", [bounds] = " << val
              << ", [" << vmin << ", " << vmax << "]"
              << std::endl;
    return false;
  }
  return true;
}

int getNumPointsFromPPD(Real min, Real max, int ppd) {
  Bounds b(min,max,3,true);
  Real ndecades = b.grid.max() - b.grid.min();
  return static_cast<int>(std::ceil(ppd*ndecades));
}
