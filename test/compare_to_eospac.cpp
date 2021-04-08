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

/*
  Comparison for SpinerEOS vs. EOSPAC
  Authors: Jonah Miller, Chad Meyers, Josh Dolence, Daniel Holiday
 */

// TODO there's a memory leak/double free somewhere

#ifdef SPINER_USE_HDF
#ifdef SINGULARITY_USE_EOSPAC

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <chrono>

#include "hdf5.h"
#include "hdf5_hl.h"

// #include "json.hpp"

#include "eos_Interface.h" // eospac API

#include "../utils/ports-of-call/portability.hpp"
#include "../utils/sp5/singularity_eos_sp5.hpp"
#include "../utils/spiner/spiner_types.hpp"
#include "../utils/spiner/sp5.hpp"
#include "../utils/spiner/databox.hpp"
#include "../utils/spiner/interpolation.hpp"

#include "../utils/sesame2spiner/io_eospac.hpp"

#include "portability.hpp"

#include "../eos/eos.hpp"

using namespace singularity;

// using nlohmann::json;
using duration = std::chrono::microseconds;
constexpr char diffFileName[] = "diffs.sp5";
constexpr int NTIMES = 10;

// #define singularity_normalize_cold_point

int main(int argc, char* argv[]) {

  herr_t status = H5_SUCCESS;

  #ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
  #endif
  {

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " filename [nFine]"
              << std::endl;
    std::exit(1);
  }
  std::string filename = argv[1];
  int nFine = argc >= 3 ? std::atoi(argv[2]) : 512;
  int nFineRho = nFine + 1;
  int nFineT = nFine - 1;
  std::string sp5name = filename;

  std::cout << "Comparison to EOSPAC" << std::endl;

  if (status != H5_SUCCESS) {
    std::cerr << "H5 error. Aborting" << std::endl;
    std::exit(1);
  }

  // params = fileToJson(filename); // reset to avoid horrible json magic
  hid_t file = H5Fcreate(diffFileName,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

  /*
  std::vector<int> matids;
  for (auto & matParams : params["materials"]) {
    matids.push_back(matParams["matid"]);
  }
  */
  std::vector<int> matids{5030};

  constexpr int NT = 3;
  EOS_INTEGER nXYPairs = (nFineRho)*(nFineT);
  DataBox xVals(nFineRho,nFineT);
  DataBox yVals(nFineRho,nFineT);
  DataBox vars(nFineRho,nFineT);
  DataBox dx(nFineRho,nFineT);
  DataBox dy(nFineRho,nFineT);
  DataBox pressEOSPAC(nFineRho,nFineT);
  DataBox pressDiff_h(nFineRho, nFineT);
  DataBox pressDiff_d(nFineRho, nFineT);

  #ifdef PORTABILITY_STRATEGY_KOKKOS
  using RView = Kokkos::View<Real*>;
  RView pressSpinerDView("pressSpinerDevice", nFineRho*nFineT);
  typename RView::HostMirror pressSpinerHView
    = Kokkos::create_mirror_view(pressSpinerDView);
  DataBox pressSpiner_d(pressSpinerDView.data(),nFineRho,nFineT);
  DataBox pressSpiner_hm(pressSpinerHView.data(),nFineRho,nFineT);
  #else
  DataBox pressSpiner_d(nFineRho,nFineT);
  DataBox pressSpiner_hm = pressSpiner_d.slice(2,0,nFineRho);
  #endif
  
  DataBox pressSpiner_h(nFineRho,nFineT);
  DataBox tempSpiner_h(nFineRho,nFineT);
  DataBox tempSpinerE_h(nFineRho,nFineT);
  DataBox tempEOSPAC(nFineRho,nFineT);
  DataBox tempDiff_h(nFineRho,nFineT);
  DataBox tempDiff_d(nFineRho,nFineT);
  DataBox tempDiffE_h(nFineRho,nFineT);
  DataBox tempDiffE_d(nFineRho,nFineT);

  #ifdef PORTABILITY_STRATEGY_KOKKOS
  RView tempSpinerDView("tempSpinerDevice", nFineRho*nFineT);
  RView tempSpinerDViewE("tempSpinerSieDevice", nFineRho*nFineT);
  RView::HostMirror tempSpinerHView
    = Kokkos::create_mirror_view(tempSpinerDView);
  RView::HostMirror tempSpinerHViewE
    = Kokkos::create_mirror_view(tempSpinerDViewE);
  DataBox tempSpiner_d(tempSpinerDView.data(),nFineRho,nFineT);
  DataBox tempSpiner_hm(tempSpinerHView.data(),nFineRho,nFineT);
  DataBox tempSpinerE_d(tempSpinerDViewE.data(),nFineRho,nFineT);
  DataBox tempSpinerE_hm(tempSpinerHViewE.data(),nFineRho,nFineT);
  RView rhos_v("rhos",nFineRho);
  RView Ts_v("Ts",nFineT);
  RView sies_v("sies",nFineT);
  DataBox rhos(rhos_v.data(),nFineRho);
  DataBox Ts(Ts_v.data(),nFineT);
  DataBox sies(sies_v.data(),nFineT);
  #else
  DataBox tempSpiner_d(nFineRho,nFineT);
  DataBox tempSpiner_hm = tempSpiner_d.slice(2,0,nFineRho);
  DataBox tempSpinerE_d(nFineRho,nFineT);
  DataBox tempSpinerE_hm = tempSpiner_h.slice(2,0,nFineRho);
  DataBox rhos(nFineRho);
  DataBox Ts(nFineT);
  DataBox sies(nFineT);
  #endif
  
  std::cout << "Comparing to EOSPAC for materials..." << std::endl;
  duration durationEospacTot    = duration::zero();
  duration durationSpinerTot    = duration::zero();
  duration durationSpinerDTot   = duration::zero();
  duration durationSpinerSieTot = duration::zero();
  duration durationEospac       = duration::zero();
  duration durationSpiner       = duration::zero();
  duration durationSpinerSie    = duration::zero();
  duration durationSpinerDev    = duration::zero();
  duration durationSpinerSieDev = duration::zero();
  auto start = std::chrono::high_resolution_clock::now();
  auto stop  = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < matids.size(); i++) {
    durationEospac = duration::zero();
    durationSpiner = duration::zero();

    int matid = matids[i];
    std::cout << "\t..." << matid << std::endl;
    std::cout << "\t------------------------" << std::endl;

    SesameMetadata metadata;
    eosGetMetadata(matid,metadata,Verbosity::Debug);

    hid_t idGroup = H5Gcreate(file, std::to_string(matid).c_str(),
                              H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);   
    
    std::cout << "\t\tloading EOSs\n"
              << "\t\t\tdepends rho T" << std::endl;
    SpinerEOSDependsRhoT eos_host(sp5name, matid);
    std::cout << "\t\t\tdepends rho sie" << std::endl;
    SpinerEOSDependsRhoSie eosE_host(sp5name, matid);
    std::cout << "\t\tmoving to device" << std::endl;
    EOS eos_device = eos_host.GetOnDevice();
    EOS eosE_device = eosE_host.GetOnDevice();

    Real* lambda_dp =
      (Real*)PORTABLE_MALLOC(sizeof(Real)*nFineRho*nFineT*eos_host.nlambda());
    Real* lambda_hp =
      (Real*)malloc(sizeof(Real)*nFineRho*nFineT*eos_host.nlambda());
    DataBox lambda_h(lambda_hp,nFineRho,nFineT,eos_host.nlambda());
    DataBox lambda_d(lambda_dp,nFineRho,nFineT,eos_host.nlambda());

    Real rhoMin = eosE_host.rhoMin();
    Real rhoMax = eosE_host.rhoMax();
    Real TMin = eosE_host.TMin();
    Real TMax = eosE_host.TMax();
    Real sieMin = eosE_host.sieMin();
    Real sieMax = eosE_host.sieMax();
    Real shrinklRhoBounds = 0.0;
    Real shrinklTBounds = 0.0;
    Real shrinkleBounds = 0.0;

    Bounds lRhoBounds(rhoMin,rhoMax,nFineRho,true,shrinklRhoBounds);
    Bounds lTBounds(TMin,TMax,nFineT,true,shrinklTBounds);
    Bounds leBounds(sieMin,sieMax,nFineT,true,shrinkleBounds);

    std::cout << "\t\trho bounds = [" << rhoMin << ", " << rhoMax << "]" << std::endl;
    std::cout << "\t\tT   bounds = [" << TMin << ", " << TMax << "]" << std::endl;
    std::cout << "\t\tsie bounds = [" << sieMin << ", " << sieMax << "]" << std::endl;
    
    EOS_INTEGER tableHandle[NT];
    EOS_INTEGER eospacPofRT, eospacTofRE, eospacEofRT;
    EOS_INTEGER tableType[NT] = {EOS_Pt_DT, EOS_T_DUt, EOS_Ut_DT};
    eosSafeLoad(NT, matid, tableType, tableHandle, Verbosity::Debug);
    eospacPofRT = tableHandle[0];
    eospacTofRE = tableHandle[1];
    eospacEofRT = tableHandle[2];

    std::cout << "\t\tGenerating interpolation points for rho-T tables"
              << std::endl;
    start = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < nFineRho; j++) {
      Real rho = lRhoBounds.i2lin(j);
      for (int i = 0; i < nFineT; i++) {
        Real T = lTBounds.i2lin(i);
        xVals(j,i) = rho;
        yVals(j,i) = T;
      }
    }
    stop = std::chrono::high_resolution_clock::now();
    durationEospac += std::chrono::duration_cast<duration>(stop-start);

    portableFor("densities",0,nFineRho,PORTABLE_LAMBDA(const int& j) {
        rhos(j) = lRhoBounds.i2lin(j);
      });
    portableFor("temperatures",0,nFineT,PORTABLE_LAMBDA(const int& i) {
        Ts(i) = lTBounds.i2lin(i);
        sies(i) = leBounds.i2lin(i);
      });

    std::cout << "\t\tInterpolate pressure of rho T" << std::endl;
    {
      pressSpiner_h.setRange(0,lTBounds.grid);
      pressSpiner_h.setRange(1,lRhoBounds.grid);
      pressSpiner_d.setRange(0,lTBounds.grid);
      pressSpiner_d.setRange(1,lRhoBounds.grid);
      pressSpiner_hm.setRange(0,lTBounds.grid);
      pressSpiner_hm.setRange(1,lRhoBounds.grid);
      pressEOSPAC.copyMetadata(pressSpiner_h);
      pressDiff_h.copyMetadata(pressSpiner_h);
      pressDiff_d.copyMetadata(pressSpiner_h);

      std::cout << "\t\t...eospac..." << std::endl;
      start = std::chrono::high_resolution_clock::now();
      for (int n = 0; n < NTIMES; n++) {
        for (int i = 0; i < nXYPairs; i++) vars(i) = dx(i) = dy(i) = 0;
        eosSafeInterpolate(&eospacPofRT, nXYPairs,
                           xVals.data(), yVals.data(), vars.data(),
                           dx.data(), dy.data(),
                           "PofRT", Verbosity::Quiet);
      }
      stop = std::chrono::high_resolution_clock::now();
      durationEospac += std::chrono::duration_cast<duration>(stop-start);
      
      std::cout << "\t\t...spiner on host..." << std::endl;
      start = std::chrono::high_resolution_clock::now();
      for (int n = 0; n < NTIMES; n++) {
        #pragma omp simd
        for (int j = 0; j < nFineRho; j++) {
          const Real rho = xVals(j,0);
          for (int i = 0; i < nFineT; i++) {
            const Real T = yVals(0,i);
            pressSpiner_h(j,i)
              = eos_host.PressureFromDensityTemperature(rho,T);
          }
        }
      }
      stop = std::chrono::high_resolution_clock::now();
      durationSpiner += std::chrono::duration_cast<duration>(stop-start);
      
      std::cout << "\t\t...spiner on device..." << std::endl;
      #ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
      #endif
      start = std::chrono::high_resolution_clock::now();
      portableFor("pressure from density and temperature",
                  0,NTIMES,0,nFineRho,0,nFineT,
                  PORTABLE_LAMBDA(const int& n, const int& j, const int& i) {
                    Real rho = rhos(j);
                    Real T = Ts(i);
                    pressSpiner_d(j,i)
                      = eos_device.PressureFromDensityTemperature(rho,T);
                  });
      #ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
      #endif
      stop = std::chrono::high_resolution_clock::now();
      durationSpinerDev += std::chrono::duration_cast<duration>(stop-start);

      Real diffL2;
      std::cout << "\t\t...comparing on host..." << std::endl;
      diffL2 = 0;
      for (int j = 0; j < nFineRho; j++) {
        for (int i = 0; i < nFineT; i++) {
          pressEOSPAC(j,i) = pressureFromSesame(vars(j,i));
          pressDiff_h(j,i) = pressSpiner_h(j,i) - pressEOSPAC(j,i);
          Real denom = std::max(std::abs(pressEOSPAC(j,i)),
                                std::abs(pressSpiner_h(j,i)));
          diffL2 += pressDiff_h(j,i)*pressDiff_h(j,i)/(denom*denom + 1e-5);
        }
      }
      diffL2 /= (nFineRho*nFineT);
      std::cout << "\t\t...L2 difference = " << diffL2 << std::endl;

      hid_t pressGroup = H5Gcreate(idGroup, "pressuresRhoT",
                                   H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status += pressSpiner_h.saveHDF(pressGroup,"spiner_h");
      status += pressEOSPAC.saveHDF(pressGroup,"EOSPAC");
      status += pressDiff_h.saveHDF(pressGroup,"diff_h");
      
      #ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::deep_copy(pressSpinerHView,pressSpinerDView);
      #endif
      std::cout << "\t\t...comparing device results..." << std::endl;
      diffL2 = 0;
      for (int j = 0; j < nFineRho; j++) {
        for (int i = 0; i < nFineT; i++) {
          // pressSpiner_d has been copied into pressSpiner_h
          pressDiff_d(j,i) = pressSpiner_hm(j,i) - pressEOSPAC(j,i);
          Real denom = std::max(std::abs(pressSpiner_hm(j,i)),
                                std::abs(pressEOSPAC(j,i)));
          diffL2 += pressDiff_d(j,i)*pressDiff_d(j,i)/(denom*denom + 1e-5);
        }
      }
      diffL2 /= (nFineRho*nFineT);
      std::cout << "\t\t...L2 difference = " << diffL2 << std::endl;

      status += pressSpiner_hm.saveHDF(pressGroup,"spiner_d");
      status += pressDiff_d.saveHDF(pressGroup,"diff_d");
      status += H5Gclose(pressGroup);      
      
      std::cout << "\t\t...EOSPAC time/point (microseconds) = "
                << durationEospac.count() / static_cast<Real>((nFineRho)*(nFineT)*NTIMES)
                << std::endl;
      std::cout << "\t\t...spiner host time/point (microseconds) = "
                << durationSpiner.count() / static_cast<Real>((nFineRho)*(nFineT)*NTIMES)
                << std::endl;
      std::cout << "\t\t...spiner device time/point (microseconds) = "
                << durationSpinerDev.count() / static_cast<Real>((nFineRho)*(nFineT)*NTIMES)
                << std::endl;
      std::cout << "\t\t...saving to file..." << std::endl;
    }
    durationEospacTot += durationEospac;
    durationSpinerTot += durationSpiner;
    durationSpinerDTot += durationSpinerDev;
    durationEospac = duration::zero();
    durationSpiner = duration::zero();
    durationSpinerDev = duration::zero();

    std::cout << "\t\tGenerating interpolation points for rho-sie tables"
              << std::endl;
    start = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < nFineRho; j++) {
      Real rho = lRhoBounds.i2lin(j);
      for (int i = 0; i < nFineT; i++) {
        Real sie = leBounds.i2lin(i);
        xVals(j,i) = rho;
        yVals(j,i) = sieToSesame(sie);
      }
    }
    stop = std::chrono::high_resolution_clock::now();
    durationEospac += std::chrono::duration_cast<duration>(stop-start);
    std::cout << "\t\tRoot find temperature of rho sie" << std::endl;
    {      
      tempSpiner_h.setRange(0,leBounds.grid);
      tempSpiner_h.setRange(1,lRhoBounds.grid);
      tempSpiner_d.setRange(0,leBounds.grid);
      tempSpiner_d.setRange(1,lRhoBounds.grid);
      tempSpiner_hm.setRange(0,leBounds.grid);
      tempSpiner_hm.setRange(1,lRhoBounds.grid);
      tempSpinerE_h.setRange(0,leBounds.grid);
      tempSpinerE_h.setRange(1,lRhoBounds.grid);
      tempSpinerE_d.setRange(0,leBounds.grid);
      tempSpinerE_d.setRange(1,lRhoBounds.grid);
      tempSpinerE_hm.setRange(0,leBounds.grid);
      tempSpinerE_hm.setRange(1,lRhoBounds.grid);
      tempEOSPAC.copyMetadata(tempSpiner_h);
      tempDiff_h.copyMetadata(tempSpiner_h);
      tempDiff_d.copyMetadata(tempSpiner_h);
      tempDiffE_h.copyMetadata(tempSpiner_h);
      tempDiffE_d.copyMetadata(tempSpiner_d);

      std::cout << "\t\t...eospac..." << std::endl;
      start = std::chrono::high_resolution_clock::now();
      for (int n = 0; n < NTIMES; n++) {
        eosSafeInterpolate(&eospacTofRE, nXYPairs,
                           xVals.data(), yVals.data(), vars.data(),
                           dx.data(), dy.data(),
                           "TofRE", Verbosity::Quiet);
      }
      stop = std::chrono::high_resolution_clock::now();
      durationEospac += std::chrono::duration_cast<duration>(stop-start);

      std::cout << "\t\t...spiner(rho,T) on host..." << std::endl;
      start = std::chrono::high_resolution_clock::now();
      for (int n = 0; n < NTIMES; n++) {
        for (int j = 0; j < nFineRho; j++) {
          Real rho = lRhoBounds.i2lin(j);
          DataBox lambda_hj = lambda_h.slice(j);
          for (int i = 0; i < nFineT; i++) {
            Real sie = leBounds.i2lin(i);
            DataBox lambda = lambda_hj.slice(i);
            tempSpiner_h(j,i)
              = eos_host.TemperatureFromDensityInternalEnergy(rho,sie,
                                                              lambda.data());
          }
        }
      }
      stop = std::chrono::high_resolution_clock::now();
      durationSpiner += std::chrono::duration_cast<duration>(stop-start);

      std::cout << "\t\t...spiner(rho,T) on device..." << std::endl;
      #ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
      #endif
      start = std::chrono::high_resolution_clock::now();
      portableFor("T(rho,e)",
                  0,NTIMES,0,nFineRho,0,nFineT,
                  PORTABLE_LAMBDA(const int& n, const int& j, const int& i) {
                    Real rho = rhos(j);
                    Real sie = sies(i);
                    DataBox lambda = lambda_d.slice(j).slice(i);
                    tempSpiner_d(j,i)
                      = eos_device.TemperatureFromDensityInternalEnergy(rho,sie,
                                                                        lambda.data());
                  });
      #ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
      #endif
      stop = std::chrono::high_resolution_clock::now();
      durationSpinerDev += std::chrono::duration_cast<duration>(stop-start);

      std::cout << "\t\t...spiner(rho,sie) on host..." << std::endl;
      start = std::chrono::high_resolution_clock::now();
      for (int n = 0; n < NTIMES; n++) {
        for (int j = 0; j < nFineRho; j++) {
          Real rho = lRhoBounds.i2lin(j);
          for (int i = 0; i < nFineT; i++) {
            Real sie = leBounds.i2lin(i);
            tempSpinerE_h(j,i)
              = eosE_host.TemperatureFromDensityInternalEnergy(rho,sie);
          }
        }
      }
      stop = std::chrono::high_resolution_clock::now();
      durationSpinerSie += std::chrono::duration_cast<duration>(stop-start);

      std::cout << "\t\t...spiner(rho,sie) on device..." << std::endl;
      #ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
      #endif
      start = std::chrono::high_resolution_clock::now();
        portableFor("T(rho,e)",
                    0,NTIMES,0,nFineRho,0,nFineT,
                    PORTABLE_LAMBDA(const int& n, const int& j, const int& i) {
                      Real rho = rhos(j);
                      Real sie = sies(i);
                      tempSpinerE_d(j,i)
                        = eosE_device.TemperatureFromDensityInternalEnergy(rho,
                                                                           sie);
                    });
      #ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
      #endif
      stop = std::chrono::high_resolution_clock::now();
      durationSpinerSieDev += std::chrono::duration_cast<duration>(stop-start);

      std::cout << "\t\t...comparing host data..." << std::endl;
      Real diffL1, diffL1E, diffL1_d, diffL1E_d;
      Real diffL2, diffL2E, diffL2_d, diffL2E_d;
      diffL2 = 0;
      for (int j = 0; j < nFineRho; j++) {
        for (int i = 0; i < nFineT; i++) {
          tempEOSPAC(j,i) = temperatureFromSesame(vars(j,i));
          tempDiff_h(j,i) = tempSpiner_h(j,i) - tempEOSPAC(j,i);
          Real denom = std::max(std::abs(tempSpiner_h(j,i)),
                                std::abs(tempEOSPAC(j,i)));
          Real diff = tempDiff_h(j,i)*tempDiff_h(j,i)/(denom*denom + 1e-5);
          if (isnan(diff)) {
            std:: cout << "NAN! " << j << ", " << i << ", "
                       << tempSpiner_h(j,i) << ", "
                       << tempEOSPAC(j,i) << ", "
                       << denom << ", "
                       << diff
                       << std::endl;
            exit(1);
          }
          #ifdef singularity_normalize_cold_point
          if (tempSpiner_h(j,i) <= 2.*TMin && tempEOSPAC(j,i) < TMin) {
            tempSpiner_h(j,i) = 0;
            tempEOSPAC(j,i) = 0;
            diff = 0;
          }
          #endif
          diffL2 += diff;
        }
      }
      diffL2 /= (nFineRho*nFineT);
      for (int j = 0; j < nFineRho; j++) {
        for (int i = 0; i < nFineT; i++) {
          tempDiffE_h(j,i) = tempSpinerE_h(j,i) - tempEOSPAC(j,i);
          Real denom = std::max(std::abs(tempSpinerE_h(j,i)),
                                std::abs(tempEOSPAC(j,i)));
          diffL2E += tempDiffE_h(j,i)*tempDiffE_h(j,i)/(denom*denom + 1e-5);
        }
      }
      diffL2E /= (nFineRho*nFineT);
      std::cout << "\t\tRoot finding:\n"
                << "\t\tits: counts\n";
      for (int i = 0; i < eos_host.counts.nBins(); i++) {
        std::cout << "\t\t\t"
                  << i << ": "
                  << 100.*eos_host.counts[i]/eos_host.counts.total()
                  << "\n";
      }
      hid_t pressGroup = H5Gcreate(idGroup, "temperaturesRhoSie",
                                   H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status += tempSpiner_h.saveHDF(pressGroup,"spiner_h");
      status += tempEOSPAC.saveHDF(pressGroup,"EOSPAC");
      status += tempDiff_h.saveHDF(pressGroup,"diff_h");
      status += tempSpinerE_h.saveHDF(pressGroup,"spinerE_h");
      status += tempDiffE_h.saveHDF(pressGroup,"diffE_h");

      #ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::deep_copy(tempSpinerHView,tempSpinerDView);
      Kokkos::deep_copy(tempSpinerHViewE,tempSpinerDViewE);
      #endif
      // tempSpiner_hm now contains copy of device data
      diffL2_d = 0;
      for (int j = 0; j < nFineRho; j++) {
        for (int i = 0; i < nFineT; i++) {
          tempDiff_d(j,i) = tempSpiner_hm(j,i) - tempEOSPAC(j,i);
          Real denom = std::max(std::abs(tempSpiner_hm(j,i)),
                                std::abs(tempEOSPAC(j,i)));
          Real diff = tempDiff_d(j,i)*tempDiff_d(j,i)/(denom*denom+1e-5);
          #ifdef singularity_normalize_cold_point
          if (tempSpiner_hm(j,i) < 2*TMin && tempEOSPAC(j,i) < TMin) {
            tempSpiner_hm(j,i) = 0;
            tempEOSPAC(j,i) = 0;
            diff = 0;
          }
          #endif
          diffL2_d += diff;
        }
      }
      diffL2_d /= (nFineRho*nFineT);
      diffL2E_d = 0;
      for (int j = 0; j < nFineRho; j++) {
        for (int i = 0; i < nFineT; i++) {
          tempDiffE_d(j,i) = tempSpinerE_hm(j,i) - tempEOSPAC(j,i);
          Real denom = std::max(std::abs(tempSpinerE_hm(j,i)),
                                std::abs(tempEOSPAC(j,i)));
          diffL2E_d += tempDiffE_d(j,i)*tempDiffE_d(j,i)/(denom*denom + 1e-5);
        }
      }
      diffL2E_d /= (nFineRho*nFineT);
      status += tempSpiner_hm.saveHDF(pressGroup,"spiner_d");
      status += tempDiff_d.saveHDF(pressGroup,"diff_d");
      status += tempSpinerE_hm.saveHDF(pressGroup,"spinerE_d");
      status += tempDiffE_d.saveHDF(pressGroup,"diffE_d");
      status += H5Gclose(pressGroup);

      std::cout << "\t------------------------" << std::endl;
      std::cout << "\t\t...L2 difference for spiner(rho,T) host = "
                << diffL2 << std::endl;
      std::cout << "\t\t...L2 difference for spiner(rho,sie) host = "
                << diffL2E << std::endl;
      std::cout << "\t\t...L2 difference for spiner(rho,T) device = "
                << diffL2_d << std::endl;
      std::cout << "\t\t...L2 difference for spiner(rho,sie) device = "
                << diffL2E_d << std::endl;
      std::cout << "\t\t...EOSPAC time/point (microseconds) = "
        << durationEospac.count() / static_cast<Real>((nFineRho)*(nFineT)*NTIMES)
                << std::endl;
      std::cout << "\t\t...spiner(rho,T) time/point (microseconds) host = "
        << durationSpiner.count() / static_cast<Real>((nFineRho)*(nFineT)*NTIMES)
                << std::endl;
      std::cout << "\t\t...spiner(rho,sie) time/point (microseconds) host = "
        << durationSpinerSie.count() / static_cast<Real>((nFineRho)*(nFineT)*NTIMES)
                << std::endl;
      std::cout << "\t\t...spiner(rho,T) time/point (microseconds) device = "
        << durationSpinerDev.count() / static_cast<Real>((nFineRho)*(nFineT)*NTIMES)
                << std::endl;
      std::cout << "\t\t...spiner(rho,sie) time/point (microseconds) device = "
        << durationSpinerSieDev.count() / static_cast<Real>((nFineRho)*(nFineT)*NTIMES)
                << std::endl;
    }
    std::cout << "\t------------------------" << std::endl;

    status += H5Gclose(idGroup);

    PORTABLE_FREE(lambda_dp);
    free(lambda_hp);
    eos_device.Finalize();
    eosE_device.Finalize();
  }

  status += H5Fclose(file);
  if (status != H5_SUCCESS) {
    std::cerr << "WARNING: problem with HDf5" << std::endl;
  }

  }
  #ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
  #endif
  
  return 0;
}

#endif /// SINGULARITY_USE_EOSPAC
#endif // SPINER_USE_HDF
