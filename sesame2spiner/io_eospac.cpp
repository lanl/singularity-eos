//======================================================================
// sesame2spiner tool for converting eospac to spiner
// Author: Jonah Miller (jonahm@lanl.gov)
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
//======================================================================

#include <algorithm>
#include <array>
#include <iostream>
#include <regex>

#include <eospac-wrapper/eospac_wrapper.hpp>

#include "io_eospac.hpp"

// TODO: more error checking of bounds?
void eosDataOfRhoSie(int matid, const Bounds &lRhoBounds, const Bounds &leBounds,
                     DataBox &Ps, DataBox &Ts, DataBox &bMods, DataBox &dPdRho,
                     DataBox &dPde, DataBox &dTdRho, DataBox &dTde, DataBox &dEdRho,
                     DataBox &mask, Verbosity eospacWarn) {
  using namespace EospacWrapper;

  constexpr int NT = 4;
  constexpr EOS_INTEGER nXYPairs = 1;
  EOS_INTEGER tableHandle[NT];
  EOS_INTEGER eospacPofRT, eospacTofRE, eospacEofRT;
  EOS_INTEGER tableType[NT] = {EOS_Pt_DT, EOS_T_DUt, EOS_Ut_DT};

  // Interpolatable vars
  EOS_REAL var[1], dx[1], dy[1];

  // indep vars
  std::vector<EOS_REAL> rhos, sies;
  makeInterpPoints(rhos, lRhoBounds);
  makeInterpPoints(sies, leBounds);

  // Load tables
  eosSafeLoad(NT, matid, tableType, tableHandle, {"EOS_Pt_DT", "EOS_T_DUt", "EOS_Ut_DT"},
              eospacWarn);
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

  // Loop by hand to ensure ordering ordering of independent
  // variables is under our control.
  for (size_t j = 0; j < rhos.size(); j++) {
    Real rho = densityToSesame(rhos[j]);
    for (size_t i = 0; i < sies.size(); i++) {
      Real sie = sieToSesame(sies[i]);
      // T
      bool no_errors = true;
      no_errors = no_errors && eosSafeInterpolate(&eospacTofRE, nXYPairs, &rho, &sie, var,
                                                  dx, dy, "TofRE", eospacWarn);
      Real T = var[0];
      Real DTDR_E = dx[0];
      Real DTDE_R = dy[0];
      // P
      no_errors = no_errors && eosSafeInterpolate(&eospacPofRT, nXYPairs, &rho, &T, var,
                                                  dx, dy, "PofRT", eospacWarn);
      Real P = var[0];
      Real DPDR_T = dx[0];
      Real DPDT_R = dy[0];
      // Bulk modulus
      no_errors = no_errors && eosSafeInterpolate(&eospacEofRT, nXYPairs, &rho, &T, var,
                                                  dx, dy, "EofRT", eospacWarn);
      Real DEDR_T = dx[0];
      Real DEDT_R = dy[0];
      Real DPDE_R = DPDT_R / DEDT_R;
      Real bMod = getBulkModulus(rho, P, DPDR_T, DPDE_R, DEDR_T);
      // Fill DataBoxes
      Ts(j, i) = temperatureFromSesame(T);
      Ps(j, i) = pressureFromSesame(P);
      bMods(j, i) = bulkModulusFromSesame(std::max(bMod, 0.0));
      dPdRho(j, i) = pressureFromSesame(DPDR_T + DTDR_E * DPDT_R);
      dPde(j, i) = sieToSesame(pressureFromSesame(DPDT_R * DTDE_R));
      dTdRho(j, i) = temperatureFromSesame(DTDR_E);
      dTde(j, i) = sieToSesame(temperatureFromSesame(DTDE_R));
      dEdRho(j, i) = densityToSesame(sieFromSesame(DEDR_T));
      mask(j, i) = no_errors ? 1.0 : 0.0;
    }
  }

  eosSafeDestroy(NT, tableHandle, eospacWarn);
}

void eosDataOfRhoT(int matid, const Bounds &lRhoBounds, const Bounds &lTBounds,
                   DataBox &Ps, DataBox &sies, DataBox &bMods, DataBox &dPdRho,
                   DataBox &dPdE, DataBox &dTdRho, DataBox &dTde, DataBox &dEdRho,
                   DataBox &dEdT, DataBox &mask, Verbosity eospacWarn) {
  using namespace EospacWrapper;

  constexpr int NT = 3;
  constexpr EOS_INTEGER nXYPairs = 1;
  EOS_INTEGER tableHandle[NT];
  EOS_INTEGER eospacPofRT, eospacTofRE, eospacEofRT;
  EOS_INTEGER tableType[NT] = {EOS_Pt_DT, EOS_T_DUt, EOS_Ut_DT};

  // Interpolatable vars
  EOS_REAL var[1], dx[1], dy[1];

  // indep vars
  std::vector<EOS_REAL> rhos, Ts;
  makeInterpPoints(rhos, lRhoBounds);
  makeInterpPoints(Ts, lTBounds);

  // Load tables
  eosSafeLoad(NT, matid, tableType, tableHandle, {"EOS_Pt_DT", "EOS_T_DUt", "EOS_Ut_DT"},
              eospacWarn);
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

  // Loop by hand to ensure ordering ordering of independent
  // variables is under our control.
  for (size_t j = 0; j < rhos.size(); j++) {
    Real rho = densityToSesame(rhos[j]);
    for (size_t i = 0; i < Ts.size(); i++) {
      Real T = temperatureToSesame(Ts[i]);
      bool no_errors = true;
      // Pressure
      no_errors = no_errors && eosSafeInterpolate(&eospacPofRT, nXYPairs, &rho, &T, var,
                                                  dx, dy, "PofRT", eospacWarn);
      Real P = var[0];
      Real DPDR_T = dx[0];
      Real DPDT_R = dy[0];
      // Energy
      no_errors = no_errors && eosSafeInterpolate(&eospacEofRT, nXYPairs, &rho, &T, var,
                                                  dx, dy, "EofRT", eospacWarn);
      Real E = var[0];
      Real DEDR_T = dx[0];
      Real DEDT_R = dy[0];
      // T derivatives
      no_errors = no_errors && eosSafeInterpolate(&eospacTofRE, nXYPairs, &rho, &E, var,
                                                  dx, dy, "TofRE", eospacWarn);
      Real DTDR_E = dx[0];
      Real DTDE_R = dy[0];
      Real DPDE_R = DPDT_R / DEDT_R;
      Real bMod = getBulkModulus(rho, P, DPDR_T, DPDE_R, DEDR_T);
      // Fill DataBoxes
      Ps(j, i) = pressureFromSesame(P);
      sies(j, i) = sieFromSesame(E);
      bMods(j, i) = bulkModulusFromSesame(std::max(bMod, 0.0));
      dPdRho(j, i) = pressureFromSesame(DPDR_T + DTDR_E * DPDT_R);
      dPdE(j, i) = sieToSesame(pressureFromSesame(DPDT_R * DTDE_R));
      dTdRho(j, i) = temperatureFromSesame(DTDR_E);
      dTde(j, i) = sieToSesame(temperatureFromSesame(DTDE_R));
      dEdRho(j, i) = densityToSesame(sieFromSesame(DEDR_T));
      dEdT(j, i) = sieFromSesame(temperatureToSesame(DEDT_R));
      mask(j, i) = no_errors ? 1.0 : 0.0;
    }
  }
  eosSafeDestroy(NT, tableHandle, eospacWarn);
}

void eosColdCurves(int matid, const Bounds &lRhoBounds, DataBox &Ps, DataBox &sies,
                   DataBox &dPdRho, DataBox &dEdRho, DataBox &bMods, DataBox &mask,
                   Verbosity eospacWarn) {
  using namespace EospacWrapper;

  constexpr int NT = 2;
  constexpr EOS_INTEGER nXYPairs = 1;
  EOS_INTEGER tableHandle[NT];
  EOS_INTEGER eospacPColdCurve, eospacSieColdCurve;
  EOS_INTEGER tableType[NT] = {EOS_Pc_D, EOS_Uc_D};

  // Interpolatable vars
  EOS_REAL var[1], dx[1], dy[1];

  // indep vars
  std::vector<EOS_REAL> rhos;
  makeInterpPoints(rhos, lRhoBounds);

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

  // loop by hand to ensure ordering
  // TODO(JMM): this is probably not necessary
  for (size_t i = 0; i < rhos.size(); i++) {
    Real rho = densityToSesame(rhos[i]);
    Real T = temperatureToSesame(0);
    bool no_errors = true;
    // pressure cold curve
    no_errors = no_errors && eosSafeInterpolate(&eospacPColdCurve, nXYPairs, &rho, &T,
                                                var, dx, dy, "PCold", eospacWarn);
    Real P = var[0];
    Real DPDR_T = dx[0];
    // energy cold curve
    no_errors = no_errors && eosSafeInterpolate(&eospacSieColdCurve, nXYPairs, &rho, &T,
                                                var, dx, dy, "sieCold", eospacWarn);
    Real sie = var[0];
    Real DEDR_T = dx[0];
    // bulk modulus cold curev
    Real bMod = getBulkModulus(rho, P, DPDR_T, 0, 0);
    // fill DataBoxes
    Ps(i) = pressureFromSesame(P);
    sies(i) = sieFromSesame(sie);
    bMods(i) = bulkModulusFromSesame(std::max(bMod, 0.0));
    dPdRho(i) = pressureFromSesame(DPDR_T); // on cold curve DTDR_E = 0
    dEdRho(i) = sieFromSesame(DEDR_T);
    mask(i) = no_errors ? 1.0 : 0.0;
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

  // loop and fill dummy variable just to see if EOSPAC errors out
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
