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

#include "io_eospac.hpp"

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

  // Interpolatable vars
  EOS_INTEGER nXYPairs = rhos.size() * sies.size();
  std::vector<EOS_REAL> T_pack(nXYPairs), P_pack(nXYPairs), sie_pack(nXYPairs),
      DTDR_E(nXYPairs), DTDE_R(nXYPairs), DPDR_T(nXYPairs), DPDT_R(nXYPairs),
      DEDR_T(nXYPairs), DEDT_R(nXYPairs), rho_flat(nXYPairs), sie_flat(nXYPairs);
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
  const bool no_errors_pofrt = eosSafeInterpolate(
      &eospacPofRT, nXYPairs, rho_flat.data(), T_pack.data(), P_pack.data(),
      DPDR_T.data(), DPDT_R.data(), "PofRT", eospacWarn);
  const bool no_errors_eofrt = eosSafeInterpolate(
      &eospacEofRT, nXYPairs, rho_flat.data(), T_pack.data(), sie_pack.data(),
      DEDR_T.data(), DEDT_R.data(), "EofRT", eospacWarn);
  const bool no_errors = no_errors_tofre && no_errors_pofrt && no_errors_eofrt;

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
      Ts(j, i) = temperatureFromSesame(T_pack[iflat]);
      Ps(j, i) = pressureFromSesame(P_pack[iflat]);
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

bool eosMassFraction(int matid, const Bounds &lRhoBounds, const Bounds &lTBounds,
                     Bounds &nphBounds, DataBox &Ms, DataBox &mask,
                     std::string &phase_names, Verbosity eospacWarn) {
  using namespace EospacWrapper;

  EOS_INTEGER errorCode = EOS_OK;
  constexpr int NT = 2;
  int ntables = NT;
  EOS_INTEGER tableHandle[NT];
  EOS_INTEGER tableType[NT] = {EOS_M_DT, EOS_Comment};
  std::vector<EOS_INTEGER> matid_v(NT, matid);
  std::vector<std::string> table_names = {"EOS_M_DT", "Eos_Material_Phases"};

  const int exists = eosCheckTableExistence(EOS_M_DT, matid, eospacWarn);
  if (exists > 0) {
    eosSafeLoad(NT, matid, tableType, tableHandle, {"EOS_M_DT", "Eos_Material_Phases"},
                eospacWarn);
  } else {
    phase_names = std::string("");
    return false;
  }
  // indep vars
  // Reuses the EOS log temp and log rho bounds
  std::vector<EOS_REAL> rhos, Ts, phs;
  makeInterpPoints(rhos, lRhoBounds);
  makeInterpPoints(Ts, lTBounds);

  // Get the number of phases
  EOS_REAL nphases_r = 1.0;
  EOS_INTEGER NITEMS[1] = {1};
  EOS_INTEGER request = EOS_NUM_PHASES;
  eos_GetTableInfo(&tableHandle[0], NITEMS, &request, &nphases_r, &errorCode);
  eosCheckError(errorCode, "eos_GetTableInfo", eospacWarn);
  const int nph = static_cast<int>(nphases_r);

  // And the phase names
  EOS_CHAR infoString[EOS_META_DATA_STRLEN];
  infoString[0] = '\0';
  EOS_INTEGER infoTypes[1] = {EOS_Material_Phases};

  eos_GetTableMetaData(&tableHandle[1], infoTypes, infoString, &errorCode);
  eosCheckError(errorCode, "eos_GetTableMetaData", eospacWarn);

  phase_names = std::string(infoString);

  DataBox dMdt, dMdr;
  Ms.resize(rhos.size(), Ts.size(), nph);
  Ms.setRange(1, lTBounds.grid);
  Ms.setRange(2, lRhoBounds.grid);

  dMdr.copyMetadata(Ms);
  dMdt.copyMetadata(Ms);
  mask.copyMetadata(Ms);

  // Interpolatable vars
  EOS_INTEGER nXYPairs = rhos.size() * Ts.size();
  std::vector<EOS_REAL> M_pack(nXYPairs * nph), DMDR_T(nXYPairs), DMDT_R(nXYPairs),
      rho_flat(nXYPairs), T_flat(nXYPairs);

  // prepare flat data structures
  for (std::size_t j = 0, iflat = 0; j < rhos.size(); ++j) {
    for (std::size_t i = 0; i < Ts.size(); ++i, iflat++) {
      rho_flat[iflat] = densityToSesame(rhos[j]);
      T_flat[iflat] = temperatureToSesame(Ts[i]);
    }
  }

  const bool no_errors = eosSafeInterpolate(&tableHandle[0], nXYPairs, rho_flat.data(),
                                            T_flat.data(), M_pack.data(), DMDR_T.data(),
                                            DMDT_R.data(), "EOS_M_DT", eospacWarn);

  for (size_t j = 0, iflat = 0; j < rhos.size(); j++) {
    for (size_t i = 0; i < Ts.size(); i++, iflat++) {
      for (size_t k = 0; k < nph; k++) {
        Ms(j, i, k) = M_pack[k * nXYPairs + iflat];
      }
    }
  }
  eosSafeDestroy(NT, tableHandle, eospacWarn);
  return true;
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
