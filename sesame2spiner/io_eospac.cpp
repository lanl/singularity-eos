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

#include <algorithm>
#include <array>
#include <iostream>
#include <regex>

#include <eospac-wrapper/eospac_wrapper.hpp>
#include <singularity-eos/eos/eos_spiner_sie_transforms.hpp>
#include "io_eospac.hpp"


TransformDataContainer::TransformDataContainer(int matid, const Bounds& lRhoBounds, const Bounds& leBounds,
    Verbosity eospacWarn)
    : lRhoOffset(lRhoBounds.offset), lEOffset(leBounds.offset) {
    using namespace EospacWrapper;
        
    //SesameMetaData metadata;
    //eosGetMetadata(matid, metadata, Verbosity::Debug);
    //Bounds lRhoBounds, lTBounds, leBounds;

//start odinary bounds since tranform is used in getmatbounds// need to ensure bounds for 
//databoxes are the same as original. must obtain same metadata. if this is implemented,
//i should be able to use this struct for ALL transforms in all files.
/*
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
    
          using namespace singularity;
    
          eosSafeInterpolate(&eospacEofRT, nXYPairs, rho, T, sie, dx, dy, "EofRT",
                             Verbosity::Quiet);
    
          const Real sieAnchor_temp = sie[0];
          const Real sieSplitPoint_temp = [1];
       
          Bounds leBounds_tranform = Bounds(Bounds::TwoGrids(), sieMin, sieMax, sieAnchor, sieSplitPoint,
              ppdSie, ppd_factor_sie, true, shrinkleBounds);

          eosSafeDestroy(NT, tableHandle, Verbosity::Quiet);
        }
    
        const Real sieAnchor = sie[0];
        const Real sieSplitPoint = sie[1]; //sie_tranformed when testing transformations
        leBounds = Bounds(Bounds::TwoGrids(), sieMin, sieMax, sieAnchor, sieSplitPoint,
                          ppdSie, ppd_factor_sie, true, shrinkleBounds);
      } else {
        leBounds = Bounds(sieMin, sieMax, numSie, true, shrinkleBounds);
      }
    
      std::cout << "lRho bounds are\n"
                << lRhoBounds << "lT bounds are\n"
                << lTBounds << "lSie bounds are \n"
                << leBounds << std::endl;
    
      return;
    }
//end bound
*/
    constexpr int NT = 1;
    EOS_INTEGER tableHandle[NT];
    EOS_INTEGER tableType[NT] = {"EOS_Uc_D"};
    EOS_INTEGER eospacSieColdCurve, eospacTofRE;

    // Create rho and T vectors
    std::vector<EOS_REAL> rhos, Ts, sies;
    makeInterpPoints(rhos, lRhoBounds);
    makeInterpPoints(sies, leBounds);
    Ts.resize(rhos.size());

    eosSafeLoad(NT, matid, tableType, tableHandle, { "EOS_Uc_D" }, eospacWarn);
    eospacSieColdCurve = tableHandle[0];
    // Interpolatable vars
    std::vector<EOS_REAL> sie_pack(rhos.size());
    std::vector<EOS_REAL> DEDR_T(rhos.size()), dy(rhos.size()); // for derivatives and dy
    DataBox sieCold;
    sieCold.resize(rhos.size());
    sieCold.setRange(0, lRhoBounds.grid);

    for (std::size_t i = 0; i < rhos.size(); ++i){
    rhos[i] = densityToSesame(rhos[i]);
    Ts[i] = temperatureToSesame(0);
    }

    EOS_INTEGER nXYPairs = rhos.size();

    eosSafeInterpolate(&eospacSieColdCurve, nXYPairs, rhos.data(), Ts.data(),
        sie_pack.data(), DEDR_T.data(), dy.data(),
        "sieCold", eospacWarn);

    for (std::size_t i = 0; i < rhos.size(); ++i)
        sieCold(i) = sieFromSesame(sie_pack[i]);  

    EOS_INTEGER tableHandleCV[NT];
    EOS_INTEGER tableTypeCV[NT] = { impl::select(split, EOS_T_DUt, EOS_T_DUe, EOS_T_DUc);}
    //Load DTDE and T databoxes
    EospacWrapper::eospacSplit apply_splitting =
        (split == TableSplit::Total) ? eospacSplit::none : eospacSplit::splitNumProp;

    std::vector<std::string> names = ["EOS_T_DUt"];
    impl::modifyNames(split, names);

    eosSafeLoad(NT, matid, tableTypeCV, tableHandleCV, names, eospacWarn, false, 0.0,
        eospacMonotonicity::none, false, apply_splitting, false);
    eospacTofRE = tableHandleCV[0];

    DataBox T, dTde, dy;
    T.resize(rhos.size(), sies.size());
    T.setRange(0, leBounds.grid);
    T.setRange(1, lRhoBounds.grid);

    dTde.resize(rhos.size(), sies.size());
    dTde.setRange(0, leBounds.grid);
    dTde.setRange(1, lRhoBounds.grid);

    dy.resize(rhos.size(), sies.size());
    dy.setRange(0, leBounds.grid);
    dy.setRange(1, lRhoBounds.grid);

    EOS_INTEGER nXYPairs = rhos.size() * sies.size();

    std::vector<EOS_REAL> T_pack(nXYPairs), DTDR_E(nXYPairs),
        DTDE_R(nXYPairs), sie_flat(nXYPairs), rho_flat(nXYPairs);;

    std::size_t iflat = 0;
    for (std::size_t j = 0; j < rhos.size(); ++j) {
        for (std::size_t i = 0; i < sies.size(); ++i) {
            rho_flat[iflat] = densityToSesame(rhos[j]);
            sie_flat[iflat] = sieToSesame(sies[i]);
            iflat++;
        }
    }

    const bool no_errors_tofre = eosSafeInterpolate(
        &eospacTofRE, nXYPairs, rho_flat.data(), sie_flat.data(), T_pack.data(),
        DTDR_E.data(), DTDE_R.data(), "TofRE", eospacWarn);

    iflat = 0;
    for (size_t j = 0; j < rhos.size(); j++) {
        Real rho = densityToSesame(rhos[j]);
        for (size_t i = 0; i < sies.size(); i++) {
            T(j, i) = temperatureFromSesame(T_pack[iflat]);
            dTde(j, i) = sieToSesame(temperatureFromSesame(DTDE_R[iflat]));
            iflat++;
        }
    }

    sieCold_ = sieCold;
    T_ = T;
    dTde_ = dTde;

    eosSafeDestroy(NT, tableHandle, eospacWarn);
}





// TODO: more error checking of bounds?
void eosDataOfRhoSie(int matid, const TableSplit split, const Bounds &lRhoBounds,
                     const Bounds &leBounds, DataBox &Ps, DataBox &Ts, DataBox &bMods,
                     DataBox &dPdRho, DataBox &dPde, DataBox &dTdRho, DataBox &dTde,
                     DataBox &dEdRho, DataBox &mask, Verbosity eospacWarn) {
  using namespace EospacWrapper;

  constexpr int NT = 3;
  EOS_INTEGER tableHandle[NT];
  EOS_INTEGER eospacPofRT, eospacTofRE, eospacEofRT;
  EOS_INTEGER tableType[NT] = {impl::select(split, EOS_Pt_DT, EOS_Pe_DT, EOS_Pic_DT),
                               impl::select(split, EOS_T_DUt, EOS_T_DUe, EOS_T_DUic),
                               impl::select(split, EOS_Ut_DT, EOS_Ue_DT, EOS_Uic_DT)};
 
  // indep vars
  std::vector<EOS_REAL> rhos, sies;
  makeInterpPoints(rhos, lRhoBounds);
  makeInterpPoints(sies, leBounds);

  EospacWrapper::eospacSplit apply_splitting =
      (split == TableSplit::Total) ? eospacSplit::none : eospacSplit::splitNumProp;

  // Load tables
  std::vector<std::string> names = {"EOS_Pt_DT", "EOS_T_DUt", "EOS_Ut_DT"};
  impl::modifyNames(split, names);
  eosSafeLoad(NT, matid, tableType, tableHandle, names, eospacWarn, false, 0.0,
              eospacMonotonicity::none, false, apply_splitting, false);
  eospacPofRT = tableHandle[0];
  eospacTofRE = tableHandle[1];
  eospacEofRT = tableHandle[2];

  // Make sie faster moving index, for easy inversion
  // TODO: is this right?
  Ps.resize(rhos.size(), sies.size());
  Ps.setRange(0, leBounds.grid);
  Ps.setRange(1, lRhoBounds.grid);
  Ts.copyMetadata(Ps);
  bMods.copyMetadata(Ps);
  dPdRho.copyMetadata(Ps);
  dPde.copyMetadata(Ps);
  dTdRho.copyMetadata(Ps);
  dTde.copyMetadata(Ps);
  dEdRho.copyMetadata(Ps);
  mask.copyMetadata(Ps);

  using namespace singularity;

  TransformDataContainer data(matid, lRhoBounds, leBounds);
  ShiftTransform<TransformDataContainer> shift(data);

  // Interpolatable vars
  EOS_INTEGER nXYPairs = rhos.size() * sies.size();
  std::vector<EOS_REAL> T_pack(nXYPairs), P_pack(nXYPairs), sie_pack(nXYPairs),
      DTDR_E(nXYPairs), DTDE_R(nXYPairs), DPDR_T(nXYPairs), DPDT_R(nXYPairs),
      DEDR_T(nXYPairs), DEDT_R(nXYPairs), rho_flat(nXYPairs), sie_flat(nXYPairs);
  std::size_t iflat = 0;
  for (std::size_t j = 0; j < rhos.size(); ++j) {
    for (std::size_t i = 0; i < sies.size(); ++i) {
      rho_flat[iflat] = densityToSesame(rhos[j]);
      sie_flat[iflat] = sieToSesame(shift.inverse(sies[i]),rhos[j]);
      iflat++;
    }
  }

  const bool no_errors_tofre = eosSafeInterpolate(
      &eospacTofRE, nXYPairs, rho_flat.data(), sie_flat.data(), T_pack.data(),
      DTDR_E.data(), DTDE_R.data(), "TofRE", eospacWarn);
  const bool no_errors_pofrt = eosSafeInterpolate(
      &eospacPofRT, nXYPairs, rho_flat.data(), T_pack.data(), P_pack.data(),
      DPDR_T.data(), DPDT_R.data(), "PofRT", eospacWarn);
  const bool no_errors_eofrt = eosSafeInterpolate(
      &eospacEofRT, nXYPairs, rho_flat.data(), T_pack.data(), sie_pack.data(),
      DEDR_T.data(), DEDT_R.data(), "EofRT", eospacWarn);
  const bool no_errors = no_errors_tofre && no_errors_pofrt && no_errors_eofrt;

    
  DataBox Ts_temp, Ps_temp;
    Ts_temp.matadata(Ps);
    Ps_temp.metadata(Ps);

  // Loop by hand to ensure ordering ordering of independent
  // variables is under our control.
  iflat = 0;
  for (size_t j = 0; j < rhos.size(); j++) {
    Real rho = densityToSesame(rhos[j]);
    for (size_t i = 0; i < sies.size(); i++) {
      Real DPDE_R = DPDT_R[iflat] / DEDT_R[iflat];
      Real bMod =
          getBulkModulus(rho, P_pack[iflat], DPDR_T[iflat], DPDE_R, DEDR_T[iflat]);
      // Fill DataBoxes
      Ts_temp(j, i) = temperatureFromSesame(T_pack[iflat]);
      Ps_temp(j, i) = pressureFromSesame(P_pack[iflat]);
      bMods(j, i) = bulkModulusFromSesame(std::max(bMod, 0.0));
      dPdRho(j, i) = pressureFromSesame(DPDR_T[iflat] + DTDR_E[iflat] * DPDT_R[iflat]);
      dPde(j, i) = sieToSesame(pressureFromSesame(DPDT_R[iflat] * DTDE_R[iflat]));
      dTdRho(j, i) = temperatureFromSesame(DTDR_E[iflat]);
      dTde(j, i) = sieToSesame(temperatureFromSesame(DTDE_R[iflat]));
      dEdRho(j, i) = densityToSesame(sieFromSesame(DEDR_T[iflat]));
      mask(j, i) = no_errors ? 1.0 : 0.0;
      iflat++;
    }
   }

  for (size_t j = 0; j < rhos.size(); j++) {
    for (size_t i = 0; i < sies.size(); i++) {
        Real lRho = spiner_common::to_log(rhos[j], lRhoBounds.offset);
        Real lE = spiner_common::to_log(shift.inverse(sies[i]), rhos[j]), leBounds.offset);
        Real ts_orig = Ts_temp.interpToReal(lRho, lE);
        Real ps_orig = Ps_temp.interpToReal(lRho, lE);
        Real letrans = spiner_common::to_log(sies[i], leBounds.offset);
        Ts(lRho, letrans) = ts_orig;
        Ps(lRho, letrans) = ps_orig;
    }
}


  eosSafeDestroy(NT, tableHandle, eospacWarn);
}

void eosDataOfRhoT(int matid, const TableSplit split, const Bounds &lRhoBounds,
                   const Bounds &lTBounds, DataBox &Ps, DataBox &sies, DataBox &bMods,
                   DataBox &dPdRho, DataBox &dPdE, DataBox &dTdRho, DataBox &dTde,
                   DataBox &dEdRho, DataBox &dEdT, DataBox &mask, Verbosity eospacWarn) {
  using namespace EospacWrapper;

  constexpr int NT = 3;
  EOS_INTEGER tableHandle[NT];
  EOS_INTEGER eospacPofRT, eospacTofRE, eospacEofRT;
  EOS_INTEGER tableType[NT] = {impl::select(split, EOS_Pt_DT, EOS_Pe_DT, EOS_Pic_DT),
                               impl::select(split, EOS_T_DUt, EOS_T_DUe, EOS_T_DUic),
                               impl::select(split, EOS_Ut_DT, EOS_Ue_DT, EOS_Uic_DT)};

  // indep vars
  std::vector<EOS_REAL> rhos, Ts;
  makeInterpPoints(rhos, lRhoBounds);
  makeInterpPoints(Ts, lTBounds);

  // Load tables
  EospacWrapper::eospacSplit apply_splitting =
      (split == TableSplit::Total) ? eospacSplit::none : eospacSplit::splitNumProp;

  std::vector<std::string> names = {"EOS_Pt_DT", "EOS_T_DUt", "EOS_Ut_DT"};
  impl::modifyNames(split, names);
  eosSafeLoad(NT, matid, tableType, tableHandle, names, eospacWarn, false, 0.0,
              eospacMonotonicity::none, false, apply_splitting, false);
  eospacPofRT = tableHandle[0];
  eospacTofRE = tableHandle[1];
  eospacEofRT = tableHandle[2];

  // Make T faster moving index, for easy inversion
  Ps.resize(rhos.size(), Ts.size());
  Ps.setRange(0, lTBounds.grid);
  Ps.setRange(1, lRhoBounds.grid);
  sies.copyMetadata(Ps);
  bMods.copyMetadata(Ps);
  dPdRho.copyMetadata(Ps);
  dPdE.copyMetadata(Ps);
  dTdRho.copyMetadata(Ps);
  dTde.copyMetadata(Ps);
  dEdRho.copyMetadata(Ps);
  dEdT.copyMetadata(Ps);
  mask.copyMetadata(Ps);

  // Interpolatable vars
  EOS_INTEGER nXYPairs = rhos.size() * Ts.size();
  std::vector<EOS_REAL> P_pack(nXYPairs), DPDR_T(nXYPairs), DPDT_R(nXYPairs),
      E_pack(nXYPairs), DEDR_T(nXYPairs), DEDT_R(nXYPairs), T_pack(nXYPairs),
      DTDR_E(nXYPairs), DTDE_R(nXYPairs), rho_flat(nXYPairs), T_flat(nXYPairs);

  // prepare flat data structures
  std::size_t iflat = 0;
  for (std::size_t j = 0; j < rhos.size(); ++j) {
    for (std::size_t i = 0; i < Ts.size(); ++i) {
      rho_flat[iflat] = densityToSesame(rhos[j]);
      T_flat[iflat] = temperatureToSesame(Ts[i]);
      iflat++;
    }
  }

  // Pressure
  const bool no_errors_prt = eosSafeInterpolate(
      &eospacPofRT, nXYPairs, rho_flat.data(), T_flat.data(), P_pack.data(),
      DPDR_T.data(), DPDT_R.data(), "PofRT", eospacWarn);
  // Energy
  const bool no_errors_ert = eosSafeInterpolate(
      &eospacEofRT, nXYPairs, rho_flat.data(), T_flat.data(), E_pack.data(),
      DEDR_T.data(), DEDT_R.data(), "EofRT", eospacWarn);
  // T derivatives
  const bool no_errors_tre = eosSafeInterpolate(
      &eospacTofRE, nXYPairs, rho_flat.data(), E_pack.data(), T_pack.data(),
      DTDR_E.data(), DTDE_R.data(), "TofRE", eospacWarn);
  const bool no_errors = no_errors_prt && no_errors_ert && no_errors_tre;

  // fill databoxes
  iflat = 0;
  for (size_t j = 0; j < rhos.size(); j++) {
    Real rho = densityToSesame(rhos[j]);
    for (size_t i = 0; i < Ts.size(); i++) {
      Real DPDE_R = DPDT_R[iflat] / DEDT_R[iflat];
      Real bMod =
          getBulkModulus(rho, P_pack[iflat], DPDR_T[iflat], DPDE_R, DEDR_T[iflat]);
      Ps(j, i) = pressureFromSesame(P_pack[iflat]);
      sies(j, i) = sieFromSesame(E_pack[iflat]);
      bMods(j, i) = bulkModulusFromSesame(std::max(bMod, 0.0));
      dPdRho(j, i) = pressureFromSesame(DPDR_T[iflat] + DTDR_E[iflat] * DPDT_R[iflat]);
      dPdE(j, i) = sieToSesame(pressureFromSesame(DPDT_R[iflat] * DTDE_R[iflat]));
      dTdRho(j, i) = temperatureFromSesame(DTDR_E[iflat]);
      dTde(j, i) = sieToSesame(temperatureFromSesame(DTDE_R[iflat]));
      dEdRho(j, i) = densityToSesame(sieFromSesame(DEDR_T[iflat]));
      dEdT(j, i) = sieFromSesame(temperatureToSesame(DEDT_R[iflat]));
      mask(j, i) = no_errors ? 1.0 : 0.0;
      iflat++;
    }
  }
  eosSafeDestroy(NT, tableHandle, eospacWarn);
}

void eosColdCurves(int matid, const Bounds &lRhoBounds, DataBox &Ps, DataBox &sies,
                   DataBox &dPdRho, DataBox &dEdRho, DataBox &bMods, DataBox &mask,
                   Verbosity eospacWarn) {
  using namespace EospacWrapper;

  constexpr int NT = 2;
  EOS_INTEGER tableHandle[NT];
  EOS_INTEGER eospacPColdCurve, eospacSieColdCurve;
  EOS_INTEGER tableType[NT] = {EOS_Pc_D, EOS_Uc_D};

  // indep vars
  std::vector<EOS_REAL> rhos, Ts;
  makeInterpPoints(rhos, lRhoBounds);
  Ts.resize(rhos.size());

  // Load tables
  eosSafeLoad(NT, matid, tableType, tableHandle, {"EOS_Pc_D", "EOS_Uc_D"}, eospacWarn);
  eospacPColdCurve = tableHandle[0];
  eospacSieColdCurve = tableHandle[1];

  // vars
  Ps.resize(rhos.size());
  Ps.setRange(0, lRhoBounds.grid);
  sies.copyMetadata(Ps);
  dPdRho.copyMetadata(Ps);
  dEdRho.copyMetadata(Ps);
  bMods.copyMetadata(Ps);
  mask.copyMetadata(Ps);

  for (std::size_t i = 0; i < rhos.size(); i++) {
    rhos[i] = densityToSesame(rhos[i]);
    Ts[i] = temperatureToSesame(0);
  }
  EOS_INTEGER nXYPairs = rhos.size();

  // Interpolatable vars
  std::vector<EOS_REAL> P_pack, DPDR_T, sie_pack, DEDR_T, dy;
  P_pack.resize(rhos.size());
  DPDR_T.resize(rhos.size());
  sie_pack.resize(rhos.size());
  DEDR_T.resize(rhos.size());
  dy.resize(rhos.size());

  // Vector EOSPAC calls
  bool no_errors = true;
  no_errors = no_errors && eosSafeInterpolate(&eospacPColdCurve, nXYPairs, rhos.data(),
                                              Ts.data(), P_pack.data(), DPDR_T.data(),
                                              dy.data(), "PCold", eospacWarn);
  no_errors = no_errors && eosSafeInterpolate(&eospacSieColdCurve, nXYPairs, rhos.data(),
                                              Ts.data(), sie_pack.data(), DEDR_T.data(),
                                              dy.data(), "sieCold", eospacWarn);

  // fill in vals
  for (std::size_t i = 0; i < rhos.size(); i++) {
    // bulk modulus cold curve
    Real bMod = getBulkModulus(rhos[i], P_pack[i], DPDR_T[i], 0, 0);
    // fill DataBoxes
    Ps(i) = pressureFromSesame(P_pack[i]);
    sies(i) = sieFromSesame(sie_pack[i]);
    bMods(i) = bulkModulusFromSesame(std::max(bMod, 0.0));
    dPdRho(i) = pressureFromSesame(DPDR_T[i]); // on cold curve DTDR_E = 0
    dEdRho(i) = sieFromSesame(DEDR_T[i]);
    mask(i) = no_errors ? 1.0 : 0.0; // TODO(JMM): Currently unused.
  }
}

void eosColdCurveMask(int matid, const Bounds &lRhoBounds, const int numSie,
                      const DataBox &sieColdCurve, DataBox &mask, Verbosity eospacWarn) {
  using namespace EospacWrapper;

  constexpr int NT = 1;
  constexpr EOS_INTEGER nXYPairs = 1;
  EOS_INTEGER tableHandle[NT];
  EOS_INTEGER eospacTofRE;
  EOS_INTEGER tableType[NT] = {EOS_T_DUt};

  // Load tables
  eosSafeLoad(NT, matid, tableType, tableHandle, {"EOS_T_DUt"}, eospacWarn);
  eospacTofRE = tableHandle[0];

  // Interpolatable vars
  EOS_REAL var[1], dx[1], dy[1];

  // rho is as expected. sie will be filled in the loop.
  std::vector<EOS_REAL> rhos;
  makeInterpPoints(rhos, lRhoBounds);

  // vars
  mask.resize(rhos.size(), numSie);
  mask.setRange(1, lRhoBounds.grid);

  // sie
  Bounds coldBounds(-1, 1, numSie, false);
  mask.setRange(0, coldBounds.grid);

  // loop and fill dummy variable just to se.e if EOSPAC errors out
  for (size_t j = 0; j < rhos.size(); j++) {
    Real rho = densityToSesame(rhos[j]);
    Real sieCold = sieColdCurve(j);
    Real sieColdAbs = std::abs(sieCold);
    Bounds sieBounds(sieCold - sieColdAbs, sieCold + sieColdAbs, numSie, false);
    for (int i = 0; i < numSie; i++) {
      Real sie = sieToSesame(sieBounds.grid.x(i));
      // T
      bool no_errors = eosSafeInterpolate(&eospacTofRE, nXYPairs, &rho, &sie, var, dx, dy,
                                          "TofRE", eospacWarn);
      mask(j, i) = no_errors ? 1.0 : 0.0;
    }
  }
  eosSafeDestroy(NT, tableHandle, eospacWarn);
}

void makeInterpPoints(std::vector<EOS_REAL> &v, const Bounds &b) {
  v.resize(b.grid.nPoints());
  for (size_t i = 0; i < v.size(); i++) {
    v[i] = b.i2lin(i);
  }
}
namespace impl {
void modifyNames(TableSplit split, std::vector<std::string> &names) {
  if (split != TableSplit::Total) {
    auto rt = std::regex("t");
    std::string newstr = (split == TableSplit::ElectronOnly) ? "e" : "ic";
    for (auto &s : names) {
      if (s != "EOS_D_PtT") {
        s = std::regex_replace(s, rt, newstr);
      }
    }
  }
}
} // namespace impl
