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

#include "io_eospac.hpp"
#include <algorithm>
#include <array>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

void eosGetMetadata(int matid, SesameMetadata &metadata, Verbosity eospacWarn) {
  constexpr int NT = 2;
  EOS_INTEGER tableHandle[NT];
  EOS_INTEGER tableType[NT] = {EOS_Info, EOS_Ut_DT};

  EOS_INTEGER commentsHandle[1];
  EOS_INTEGER commentsType[1] = {EOS_Comment};

  constexpr int numInfoTables = 2;
  constexpr int NI[] = {5, 11};
  std::array<std::vector<EOS_INTEGER>, numInfoTables> infoItems = {
      std::vector<EOS_INTEGER>{EOS_Exchange_Coeff, EOS_Mean_Atomic_Mass,
                               EOS_Mean_Atomic_Num, EOS_Modulus, EOS_Normal_Density},
      std::vector<EOS_INTEGER>{EOS_Rmin, EOS_Rmax, EOS_Tmin, EOS_Tmax, EOS_Fmin, EOS_Fmax,
                               EOS_NR, EOS_NT, EOS_X_Convert_Factor, EOS_Y_Convert_Factor,
                               EOS_F_Convert_Factor}};
  std::vector<EOS_REAL> infoVals[numInfoTables];
  for (int i = 0; i < numInfoTables; i++) {
    infoVals[i].resize(NI[i]);
    for (int j = 0; j < NI[i]; j++) {
      infoVals[i][j] = 0;
    }
  }

  eosSafeLoad(NT, matid, tableType, tableHandle, {"EOS_Info", "EOS_Ut_DT"}, eospacWarn);

  for (int i = 0; i < numInfoTables; i++) {
    eosSafeTableInfo(&(tableHandle[i]), NI[i], infoItems[i].data(), infoVals[i].data(),
                     eospacWarn);
  }
  metadata.matid = matid;
  metadata.exchangeCoefficient = infoVals[0][0];
  metadata.meanAtomicMass = infoVals[0][1];
  metadata.meanAtomicNumber = infoVals[0][2];
  metadata.solidBulkModulus = bulkModulusFromSesame(infoVals[0][3]);
  metadata.normalDensity = densityFromSesame(infoVals[0][4]);
  metadata.rhoMin = densityFromSesame(infoVals[1][0]);
  metadata.rhoMax = densityFromSesame(infoVals[1][1]);
  metadata.TMin = temperatureFromSesame(infoVals[1][2]);
  metadata.TMax = temperatureFromSesame(infoVals[1][3]);
  metadata.sieMin = sieFromSesame(infoVals[1][4]);
  metadata.sieMax = sieFromSesame(infoVals[1][5]);
  metadata.numRho = static_cast<int>(infoVals[1][6]);
  metadata.numT = static_cast<int>(infoVals[1][7]);
  metadata.rhoConversionFactor = infoVals[1][8];
  metadata.TConversionFactor = infoVals[1][9];
  metadata.sieConversionFactor = infoVals[1][10];

  eosSafeDestroy(NT, tableHandle, eospacWarn);

  EOS_INTEGER errorCode =
      eosSafeLoad(1, matid, commentsType, commentsHandle, {"EOS_Comments"}, eospacWarn);
  EOS_INTEGER eospacComments = commentsHandle[0];

  if (errorCode == EOS_OK) {
    std::vector<EOS_CHAR> comments;
    EOS_REAL commentLen;
    EOS_INTEGER commentItem = EOS_Cmnt_Len;
    eosSafeTableInfo(commentsHandle, 1, &commentItem, &commentLen, eospacWarn);

    comments.resize(static_cast<int>(commentLen));
    metadata.comments.resize(comments.size());
    eosSafeTableCmnts(&eospacComments, comments.data(), eospacWarn);
    for (size_t i = 0; i < comments.size(); i++) {
      metadata.comments[i] = comments[i];
    }
    metadata.name = getName(metadata.comments);

    eosSafeDestroy(1, commentsHandle, eospacWarn);
  } else {
    std::string matid_str = std::to_string(matid);
    if (eospacWarn != Verbosity::Quiet) {
      std::cerr << "eos_GetMetadata: failed to get comments table. "
                << "Using default comments and name fields." << std::endl;
    }
    metadata.name = "No name for matid " + matid_str;
    metadata.comments = "Comment unavailable for matid " + matid_str;
  }
}

// TODO: more error checking of bounds?
void eosDataOfRhoSie(int matid, const Bounds &lRhoBounds, const Bounds &leBounds,
                     DataBox &Ps, DataBox &Ts, DataBox &bMods, DataBox &dPdRho,
                     DataBox &dPde, DataBox &dTdRho, DataBox &dTde, DataBox &dEdRho,
                     DataBox &mask, Verbosity eospacWarn) {

  constexpr int NT = 3;
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

EOS_INTEGER eosSafeLoad(int ntables, int matid, EOS_INTEGER tableType[],
                        EOS_INTEGER tableHandle[], Verbosity eospacWarn) {
  std::vector<std::string> empty;
  return eosSafeLoad(ntables, matid, tableType, tableHandle, empty, eospacWarn);
}

EOS_INTEGER eosSafeLoad(int ntables, int matid, EOS_INTEGER tableType[],
                        EOS_INTEGER tableHandle[],
                        const std::vector<std::string> &table_names,
                        Verbosity eospacWarn) {

  EOS_INTEGER NTABLES[] = {ntables};
  std::vector<EOS_INTEGER> MATID(ntables, matid);

  EOS_INTEGER errorCode = EOS_OK;
  EOS_INTEGER tableHandleErrorCode = EOS_OK;
  EOS_CHAR errorMessage[EOS_MaxErrMsgLen];

  eos_CreateTables(NTABLES, tableType, MATID.data(), tableHandle, &errorCode);

#ifdef SINGULARITY_INVERT_AT_SETUP
  EOS_INTEGER options[] = {EOS_INVERT_AT_SETUP, EOS_INSERT_DATA};
  EOS_REAL values[] = {1., 4.};
  for (int i = 0; i < ntables; i++) {
    eos_SetOption(&(tableHandle[i]), &(options[0]), &(values[0]), &errorCode);
    eos_SetOption(&(tableHandle[i]), &(options[1]), &(values[1]), &errorCode);
  }
#endif

#ifdef SINGULARITY_EOS_SKIP_EXTRAP
  for (int i = 0; i < ntables; i++) {
    eos_SetOption(&(tableHandle[i]), &EOS_SKIP_EXTRAP_CHECK, NULL, &errorCode);
  }
#endif

  eos_LoadTables(&ntables, tableHandle, &errorCode);
  if (errorCode != EOS_OK && eospacWarn != Verbosity::Quiet) {
    for (int i = 0; i < ntables; i++) {
      eos_GetErrorCode(&tableHandle[i], &tableHandleErrorCode);
      eos_GetErrorMessage(&tableHandleErrorCode, errorMessage);
      std::cerr << "eos_CreateTables ERROR " << tableHandleErrorCode;
      if (table_names.size() > 0) {
        std::cerr << " for table names\n\t{";
        for (auto &name : table_names) {
          std::cerr << name << ", ";
        }
        std::cerr << "}";
      }
      std::cerr << ":\n\t" << errorMessage << std::endl;
    }
  }
  return errorCode;
}

bool eosSafeInterpolate(EOS_INTEGER *table, EOS_INTEGER nxypairs, EOS_REAL xVals[],
                        EOS_REAL yVals[], EOS_REAL var[], EOS_REAL dx[], EOS_REAL dy[],
                        const char tablename[], Verbosity eospacWarn) {
  EOS_INTEGER errorCode = EOS_OK;
  EOS_CHAR errorMessage[EOS_MaxErrMsgLen];

  eos_Interpolate(table, &nxypairs, xVals, yVals, var, dx, dy, &errorCode);
#ifndef SINGULARITY_EOS_SKIP_EXTRAP
  if (errorCode != EOS_OK && eospacWarn == Verbosity::Debug) {
    eos_GetErrorMessage(&errorCode, errorMessage);
    std::cerr << "Table " << tablename << ":" << std::endl;
    std::cerr << "eos_Interpolate ERROR " << errorCode << ": " << errorMessage
              << std::endl;
    std::vector<EOS_INTEGER> xyBounds(nxypairs);
    eos_CheckExtrap(table, &nxypairs, xVals, yVals, xyBounds.data(), &errorCode);
    for (size_t i = 0; i < xyBounds.size(); i++) {
      std::string status = eosErrorString(xyBounds[i]);
      std::cerr << "var " << i << ": " << status << std::endl;
      std::cerr << "x = " << xVals[i] << std::endl;
      std::cerr << "y = " << yVals[i] << std::endl;
    }
  }
#endif
  return (errorCode == EOS_OK); // 1 for no erros. 0 for errors.
}

void eosSafeTableInfo(EOS_INTEGER *table, EOS_INTEGER numInfoItems,
                      EOS_INTEGER infoItems[], EOS_REAL infoVals[],
                      Verbosity eospacWarn) {
  EOS_INTEGER errorCode = EOS_OK;
  // EOS_CHAR errorMessage[EOS_MaxErrMsgLen];
  EOS_INTEGER NITEMS[] = {numInfoItems};
  eos_GetTableInfo(table, NITEMS, infoItems, infoVals, &errorCode);
  eosCheckError(errorCode, "eos_GetTableInfo", eospacWarn);
}

void eosSafeTableCmnts(EOS_INTEGER *table, EOS_CHAR *comments, Verbosity eospacWarn) {
  EOS_INTEGER errorCode = EOS_OK;
  eos_GetTableCmnts(table, comments, &errorCode);
  eosCheckError(errorCode, "eos_GetTableCmnts", eospacWarn);
}

void eosCheckError(EOS_INTEGER errorCode, const std::string &name, Verbosity eospacWarn) {
  EOS_CHAR errorMessage[EOS_MaxErrMsgLen];
  if (errorCode != EOS_OK && eospacWarn != Verbosity::Quiet) {
    eos_GetErrorMessage(&errorCode, errorMessage);
    std::cerr << name << " ERROR " << errorCode << ":\n\t" << errorMessage << std::endl;
  }
}

std::string eosErrorString(EOS_INTEGER errorCode) {
  // Algorithmicallly generated by parsing the EOSPAC docs
  // I'm sorry. It's gross. ~JMM
  switch (errorCode) {
  case EOS_OK:
    return "EOS_OK";
  case EOS_BAD_DATA_TYPE:
    return "EOS_BAD_DATA_TYPE";
  case EOS_BAD_DERIVATIVE_FLAG:
    return "EOS_BAD_DERIVATIVE_FLAG";
  case EOS_BAD_INTERPOLATION_FLAG:
    return "EOS_BAD_INTERPOLATION_FLAG";
  case EOS_BAD_MATERIAL_ID:
    return "EOS_BAD_MATERIAL_ID";
  case EOS_CANT_INVERT_DATA:
    return "EOS_CANT_INVERT_DATA";
  case EOS_CANT_MAKE_MONOTONIC:
    return "EOS_CANT_MAKE_MONOTONIC";
  case EOS_CONVERGENCE_FAILED:
    return "EOS_CONVERGENCE_FAILED";
  case EOS_DATA_TYPE_NOT_FOUND:
    return "EOS_DATA_TYPE_NOT_FOUND";
  case EOS_DATA_TYPE_NO_MATCH:
    return "EOS_DATA_TYPE_NO_MATCH";
  case EOS_FAILED:
    return "EOS_FAILED";
  case EOS_INTEGRATION_FAILED:
    return "EOS_INTEGRATION_FAILED";
  case EOS_INTERP_EXTRAPOLATED:
    return "EOS_INTERP_EXTRAPOLATED";
  case EOS_INTERP_EXTRAP_PBAL:
    return "EOS_INTERP_EXTRAP_PBAL";
  case EOS_INTERP_EXTRAP_TBAL:
    return "EOS_INTERP_EXTRAP_TBAL";
  case EOS_INVALID_CONC_SUM:
    return "EOS_INVALID_CONC_SUM";
  case EOS_INVALID_DATA_TYPE:
    return "EOS_INVALID_DATA_TYPE";
  case EOS_INVALID_INFO_FLAG:
    return "EOS_INVALID_INFO_FLAG";
  case EOS_INVALID_OPTION_FLAG:
    return "EOS_INVALID_OPTION_FLAG";
  case EOS_INVALID_SUBTABLE_INDEX:
    return "EOS_INVALID_SUBTABLE_INDEX";
  case EOS_INVALID_TABLE_HANDLE:
    return "EOS_INVALID_TABLE_HANDLE";
  case EOS_MATERIAL_NOT_FOUND:
    return "EOS_MATERIAL_NOT_FOUND";
  case EOS_MEM_ALLOCATION_FAILED:
    return "EOS_MEM_ALLOCATION_FAILED";
  case EOS_NOT_ALLOCATED:
    return "EOS_NOT_ALLOCATED";
  case EOS_NOT_INITIALIZED:
    return "EOS_NOT_INITIALIZED";
  case EOS_NO_COMMENTS:
    return "EOS_NO_COMMENTS";
  case EOS_NO_DATA_TABLE:
    return "EOS_NO_DATA_TABLE";
  case EOS_NO_SESAME_FILES:
    return "EOS_NO_SESAME_FILES";
  case EOS_OPEN_SESAME_FILE_FAILED:
    return "EOS_OPEN_SESAME_FILE_FAILED";
  case EOS_READ_DATA_FAILED:
    return "EOS_READ_DATA_FAILED";
  case EOS_READ_FILE_VERSION_FAILED:
    return "EOS_READ_FILE_VERSION_FAILED";
  case EOS_READ_MASTER_DIR_FAILED:
    return "EOS_READ_MASTER_DIR_FAILED";
  case EOS_READ_MATERIAL_DIR_FAILED:
    return "EOS_READ_MATERIAL_DIR_FAILED";
  case EOS_READ_TOTAL_MATERIALS_FAILED:
    return "EOS_READ_TOTAL_MATERIALS_FAILED";
  case EOS_SPLIT_FAILED:
    return "EOS_SPLIT_FAILED";
  case EOS_xHi_yHi:
    return "EOS_xHi_yHi";
  case EOS_xHi_yOk:
    return "EOS_xHI_yOk";
  case EOS_xHi_yLo:
    return "EOS_xHi_yLo";
  case EOS_xOk_yLo:
    return "EOS_xOk_yLo";
  case EOS_xLo_yLo:
    return "EOS_xLo_yLo";
  case EOS_xLo_yOk:
    return "EOS_xLo_yOk";
  case EOS_xLo_yHi:
    return "EOS_xLo_yHi";
  case EOS_xOk_yHi:
    return "EOS_xOk_yHi";
  case EOS_WARNING:
    return "EOS_WARNING";
  case EOS_UNDEFINED:
    return "EOS_UNDEFINED";
  default:
    return "UNKNOWN ERROR: " + std::to_string(errorCode);
  }
}

void eosSafeDestroy(int ntables, EOS_INTEGER tableHandles[], Verbosity eospacWarn) {
  EOS_INTEGER errorCode = EOS_OK;
  // EOS_CHAR errorMessage[EOS_MaxErrMsgLen];
  EOS_INTEGER NTABLES[] = {ntables};
  eos_DestroyTables(NTABLES, tableHandles, &errorCode);
  eosCheckError(errorCode, "eos_DestroyTables", eospacWarn);
}

void makeInterpPoints(std::vector<EOS_REAL> &v, const Bounds &b) {
  v.resize(b.grid.nPoints());
  for (size_t i = 0; i < v.size(); i++) {
    v[i] = b.i2lin(i);
  }
}

std::string getName(std::string comment) {
  // required because lookbehinds not supported
  // constexpr int startPos=15;
  std::regex r("101: material. (.*)(?=\\s+\\(z)");
  std::smatch match;
  if (std::regex_search(comment, match, r)) {
    return std::string(match[0]).substr(15);
  }
  std::string badstring("-1");
  return badstring;
}
