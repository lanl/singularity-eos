//------------------------------------------------------------------------------
// © 2022. Triad National Security, LLC. All rights reserved.  This
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

// This is a simple utility program to print out a full thermodynamic
// state using both eospac and SpinerEOS for a material of a given
// matid. For simplicity, the utility is CPU-only.

// C++ headers
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

// Needed to import the eos models
#include <singularity-eos/eos/eos.hpp>

// For simplicity. Most things live in the singularity namespace
using namespace singularity;

constexpr int iPC = 0;
constexpr int iRT = 1;
constexpr int iRE = 2;

inline void PrintStateOfRhoT(const double rho, const double T, EOS eos[3]) {
  Real press[3];
  Real dp[2];
  Real sie[3];
  Real dsie[2];
  Real cv[3];
  Real dcv[2];
  Real bmod[3];
  Real dbmod[2];
  Real gam[3];
  Real dgam[2];
  for (int i = 0; i < 3; ++i) {
    press[i] = eos[i].PressureFromDensityTemperature(rho, T);
    sie[i] = eos[i].InternalEnergyFromDensityTemperature(rho, T);
    cv[i] = eos[i].SpecificHeatFromDensityTemperature(rho, T);
    bmod[i] = eos[i].BulkModulusFromDensityTemperature(rho, T);
    gam[i] = eos[i].GruneisenParamFromDensityTemperature(rho, T);
  }
  for (int i = 0; i < 2; ++i) {
    dp[i] = press[i + 1] - press[0];
    dsie[i] = sie[i + 1] - sie[0];
    dcv[i] = cv[i + 1] - cv[0];
    dbmod[i] = bmod[i + 1] - bmod[0];
    dgam[i] = gam[i + 1] - gam[0];
  }
  printf("# =================================================================== #\n"
         "# State of rho, T                                                     #\n"
         "# =================================================================== #\n"
         "# rho = %.8e \n"
         "# T   = %.8e \n"
         "# =================================================================== #\n"
         "#  P          | sie         | cv          | bmod        | gamma       #\n"
         "# =================================================================== #\n"
         "# eospac                                                              #\n"
         "# ------------------------------------------------------------------- #\n"
         "| %-11.4e | %-11.4e | %-11.4e | %-11.4e | %-11.4e |\n"
         "# ------------------------------------------------------------------- #\n"
         "# Spiner(rho, T)                                                      #\n"
         "# ------------------------------------------------------------------- #\n"
         "| %-11.4e | %-11.4e | %-11.4e | %-11.4e | %-11.4e |\n"
         "# ------------------------------------------------------------------- #\n"
         "# difference                                                          #\n"
         "# ------------------------------------------------------------------- #\n"
         "| %-11.4e | %-11.4e | %-11.4e | %-11.4e | %-11.4e |\n"
         "# ------------------------------------------------------------------- #\n"
         "# Spiner(rho, sie)                                                    #\n"
         "# ------------------------------------------------------------------- #\n"
         "| %-11.4e | %-11.4e | %-11.4e | %-11.4e | %-11.4e |\n"
         "# ------------------------------------------------------------------- #\n"
         "# difference                                                          #\n"
         "# ------------------------------------------------------------------- #\n"
         "| %-11.4e | %-11.4e | %-11.4e | %-11.4e | %-11.4e |\n"
         "# ------------------------------------------------------------------- #\n",
         rho, T, press[iPC], sie[iPC], cv[iPC], bmod[iPC], gam[iPC], press[iRT], sie[iRT],
         cv[iRT], bmod[iRT], gam[iRT], dp[0], dsie[0], dcv[0], dbmod[0], dgam[0],
         press[iRE], sie[iRE], cv[iRE], bmod[iRE], gam[iRE], dp[1], dsie[1], dcv[1],
         dbmod[1], dgam[1]);
}

inline void PrintStateOfRhoSie(const double rho, const double sie, EOS eos[3]) {
  Real press[3];
  Real dp[2];
  Real T[3];
  Real dT[2];
  Real cv[3];
  Real dcv[2];
  Real bmod[3];
  Real dbmod[2];
  Real gam[3];
  Real dgam[2];
  for (int i = 0; i < 3; ++i) {
    press[i] = eos[i].PressureFromDensityInternalEnergy(rho, sie);
    T[i] = eos[i].TemperatureFromDensityInternalEnergy(rho, sie);
    cv[i] = eos[i].SpecificHeatFromDensityInternalEnergy(rho, sie);
    bmod[i] = eos[i].BulkModulusFromDensityInternalEnergy(rho, sie);
    gam[i] = eos[i].GruneisenParamFromDensityInternalEnergy(rho, sie);
  }
  for (int i = 0; i < 2; ++i) {
    dp[i] = press[i + 1] - press[0];
    dT[i] = T[i + 1] - T[0];
    dcv[i] = cv[i + 1] - cv[0];
    dbmod[i] = bmod[i + 1] - bmod[0];
    dgam[i] = gam[i + 1] - gam[0];
  }
  printf("# =================================================================== #\n"
         "# State of rho, sie                                                   #\n"
         "# =================================================================== #\n"
         "# rho = %.8e \n"
         "# sie = %.8e \n"
         "# =================================================================== #\n"
         "#  P          | T           | cv          | bmod        | gamma       #\n"
         "# =================================================================== #\n"
         "# eospac                                                              #\n"
         "# ------------------------------------------------------------------- #\n"
         "| %-11.4e | %-11.4e | %-11.4e | %-11.4e | %-11.4e |\n"
         "# ------------------------------------------------------------------- #\n"
         "# Spiner(rho, T)                                                      #\n"
         "# ------------------------------------------------------------------- #\n"
         "| %-11.4e | %-11.4e | %-11.4e | %-11.4e | %-11.4e |\n"
         "# ------------------------------------------------------------------- #\n"
         "# difference                                                          #\n"
         "# ------------------------------------------------------------------- #\n"
         "| %-11.4e | %-11.4e | %-11.4e | %-11.4e | %-11.4e |\n"
         "# ------------------------------------------------------------------- #\n"
         "# Spiner(rho, sie)                                                    #\n"
         "# ------------------------------------------------------------------- #\n"
         "| %-11.4e | %-11.4e | %-11.4e | %-11.4e | %-11.4e |\n"
         "# ------------------------------------------------------------------- #\n"
         "# difference                                                          #\n"
         "# ------------------------------------------------------------------- #\n"
         "| %-11.4e | %-11.4e | %-11.4e | %-11.4e | %-11.4e |\n"
         "# ------------------------------------------------------------------- #\n",
         rho, sie, press[iPC], T[iPC], cv[iPC], bmod[iPC], gam[iPC], press[iRT], T[iRT],
         cv[iRT], bmod[iRT], gam[iRT], dp[0], dT[0], dcv[0], dbmod[0], dgam[0],
         press[iRE], T[iRE], cv[iRE], bmod[iRE], gam[iRE], dp[1], dT[1], dcv[1], dbmod[1],
         dgam[1]);
}

int main(int argc, char *argv[]) {
  // In case we're doing a Kokkos backend
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
#endif
  { // begin kokkos scope
    if (argc < 6) {
      std::cerr << "Usage: " << argv[0] << " matid sp5_filename rho T sie" << std::endl;
      // std::cerr << "Arguments given: " << argc << std::endl;
      std::exit(1);
    }

    const int matid = std::atoi(argv[1]);
    const std::string sp5_filename = argv[2];
    const Real rho = std::atof(argv[3]);
    const Real T = std::atof(argv[4]);
    const Real sie = std::atof(argv[5]);

    EOS eos[3];
    eos[iPC] = EOSPAC(matid);
    eos[iRT] = SpinerEOSDependsRhoT(sp5_filename, matid);
    eos[iRE] = SpinerEOSDependsRhoSie(sp5_filename, matid);
    PrintStateOfRhoT(rho, T, eos);
    std::cout << std::endl;
    PrintStateOfRhoSie(rho, sie, eos);
  }
}
