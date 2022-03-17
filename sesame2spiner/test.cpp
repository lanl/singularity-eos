//======================================================================
// sesame2spiner tool for converting eospac to spiner
// Author: Jonah Miller (jonahm@lanl.gov)
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
//======================================================================

#include <cmath>
#include <iostream>
#include <string>

#include <hdf5.h>
#include <hdf5_hl.h>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#ifndef SPINER_USE_HDF
#error "HDF5 must be enabled"
#endif

#include <sp5/singularity_eos_sp5.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/ports-of-call/portability.hpp>
#include <spiner/sp5.hpp>
#include <spiner/spiner_types.hpp>

#include "generate_files.hpp"
#include "io_eospac.hpp"
#include "parser.hpp"

using Spiner::DataBox;
using Spiner::RegularGrid1D;

constexpr int matid = 5030;  // dry air
constexpr int matid2 = 4272; // stainless steel
constexpr int NRHO = 64;
constexpr int NUMT = 32;
constexpr int NUME = NUMT;
constexpr Real EPS = 1e-5;

inline bool isClose(Real a, Real b) { return fabs((b - a) / (a + b + 1e-20)) <= EPS; }

SCENARIO("Reading metadata for air through eospac", "[eospac],[eosGetMetadata]") {
  SesameMetadata m;

  GIVEN("Metadata is read in via eospac") {
    eosGetMetadata(matid, m, Verbosity::Verbose);
    THEN("Known quantities are as expected") {
      REQUIRE(isClose(m.exchangeCoefficient, 0.0));
      REQUIRE(isClose(m.meanAtomicMass, 1.480303959999998e1));
      REQUIRE(isClose(m.meanAtomicNumber, 7.372955340000004));
      REQUIRE(isClose(m.normalDensity, 1.293000000000002e-3));
      REQUIRE(isClose(m.rhoMin, 1e-7));
      REQUIRE(isClose(m.rhoMax, 15));
      REQUIRE(isClose(m.TMin, 175.235));
      REQUIRE(isClose(m.TMax, 3.4815e8));
      REQUIRE(isClose(m.sieMin, 8.402e8));
      REQUIRE(isClose(m.sieMax, 2.484e16));
      REQUIRE(m.numRho == 21);
      REQUIRE(m.numT == 31);
      REQUIRE(m.name == "dry air");
    }
  }
}

SCENARIO("Reading air data for indep. variables rho, sie", "[eospac],[eosDataOfRhoSie]") {

  DataBox P, T, bMods, dPdRho, dPde, dTdRho, dTde, dEdRho;
  DataBox mask;
  // Bounds assumed to be in log
  // TODO: these bounds are totally arbitrary
  //       BUT they are experimentally verified to be
  //       on the SESAME tables.
  Bounds lRhoBounds(-2, 1, NRHO);
  Bounds leBounds(12, 16, NUME);

  GIVEN("DataBoxes are filled by eospac") {
    eosDataOfRhoSie(matid, lRhoBounds, leBounds, P, T, bMods, dPdRho, dPde, dTdRho, dTde,
                    dEdRho, mask, Verbosity::Verbose);
    THEN("The shapes are correct") {
      REQUIRE(P.dim(1) == NUME);
      REQUIRE(P.dim(2) == NRHO);
      REQUIRE(T.dim(1) == NUME);
      REQUIRE(T.dim(2) == NRHO);
      REQUIRE(bMods.dim(1) == NUME);
      REQUIRE(bMods.dim(2) == NRHO);
      REQUIRE(dPdRho.dim(1) == NUME);
      REQUIRE(dPdRho.dim(2) == NRHO);
      REQUIRE(dPde.dim(1) == NUME);
      REQUIRE(dPde.dim(2) == NRHO);
      REQUIRE(dTdRho.dim(1) == NUME);
      REQUIRE(dTdRho.dim(2) == NRHO);
      REQUIRE(dTde.dim(1) == NUME);
      REQUIRE(dTde.dim(2) == NRHO);

      AND_THEN("The tables contain all finite values") {
        bool nansPresent = false;
        for (int j = 0; j < NRHO; j++) {
          for (int i = 0; i < NUME; i++) {
            nansPresent |= std::isnan(P(j, i));
            nansPresent |= std::isnan(T(j, i));
            nansPresent |= std::isnan(bMods(j, i));
            nansPresent |= std::isnan(dPdRho(j, i));
            nansPresent |= std::isnan(dPdRho(j, i));
            nansPresent |= std::isnan(dPde(j, i));
            nansPresent |= std::isnan(dTdRho(j, i));
            nansPresent |= std::isnan(dTde(j, i));
          }
        }
        REQUIRE(!nansPresent);
      }
    }
  }
}

SCENARIO("Reading air data for indep. variables rho, T", "[eospac],[eosDataOfRhoT]") {

  DataBox P, sie, bMod, dPdRho, dPdE, dTdRho, dTdE, dEdRho, dEdT, mask;
  Bounds lRhoBounds(-4, 1, NRHO);
  Bounds lTBounds(2.4, 4, NUMT);

  GIVEN("DataBoxes are filled by eospac") {
    eosDataOfRhoT(matid, lRhoBounds, lTBounds, P, sie, bMod, dPdRho, dPdE, dTdRho, dTdE,
                  dEdRho, dEdT, mask, Verbosity::Verbose);
    THEN("The shapes are correct") {
      REQUIRE(P.dim(1) == NUMT);
      REQUIRE(P.dim(2) == NRHO);
      AND_THEN("The tables contain all finite values") {
        bool nansPresent = false;
        for (int j = 0; j < NRHO; j++) {
          for (int i = 0; i < NUMT; i++) {
            nansPresent |= std::isnan(P(j, i));
            nansPresent |= std::isnan(sie(j, i));
            nansPresent |= std::isnan(bMod(j, i));
            nansPresent |= std::isnan(dPdRho(j, i));
            nansPresent |= std::isnan(dPdRho(j, i));
            nansPresent |= std::isnan(dPdE(j, i));
            nansPresent |= std::isnan(dTdRho(j, i));
            nansPresent |= std::isnan(dTdE(j, i));
          }
        }
        REQUIRE(!nansPresent);
      }
    }
  }
}

SCENARIO("Reading air data for all indep. variables one after the other",
         "[eospac],[eosDataOfRhoT],[eosDataOfRhoSie]") {
  DataBox P, T, sie, bMod, dPdRho, dPdE, dTdRho, dTdE, dEdRho, dEdT, mask;
  Bounds lRhoBounds(-2, 1, NRHO);
  Bounds lTBounds(2.4, 4, NUMT);
  Bounds leBounds(12, 16, NUME);
  GIVEN("Databoxes filled for log rho and log T") {
    eosDataOfRhoT(matid, lRhoBounds, lTBounds, P, sie, bMod, dPdRho, dPdE, dTdRho, dTdE,
                  dEdRho, dEdT, mask, Verbosity::Verbose);
    THEN("The shapes are correct") {
      REQUIRE(P.dim(1) == NUMT);
      REQUIRE(P.dim(2) == NRHO);
      REQUIRE(sie.dim(1) == NUMT);
      REQUIRE(sie.dim(2) == NRHO);
      REQUIRE(bMod.dim(1) == NUMT);
      REQUIRE(bMod.dim(2) == NRHO);
      REQUIRE(dPdRho.dim(1) == NUMT);
      REQUIRE(dPdRho.dim(2) == NRHO);
      REQUIRE(dPdE.dim(1) == NUMT);
      REQUIRE(dPdE.dim(2) == NRHO);
      REQUIRE(dTdRho.dim(1) == NUMT);
      REQUIRE(dTdRho.dim(2) == NRHO);
      REQUIRE(dTdE.dim(1) == NUMT);
      REQUIRE(dTdE.dim(2) == NRHO);

      AND_THEN("The tables contain all finite values") {
        bool nansPresent = false;
        for (int j = 0; j < NRHO; j++) {
          for (int i = 0; i < NUMT; i++) {
            nansPresent |= std::isnan(P(j, i));
            nansPresent |= std::isnan(sie(j, i));
            nansPresent |= std::isnan(bMod(j, i));
            nansPresent |= std::isnan(dPdRho(j, i));
            nansPresent |= std::isnan(dPdRho(j, i));
            nansPresent |= std::isnan(dPdE(j, i));
            nansPresent |= std::isnan(dTdRho(j, i));
            nansPresent |= std::isnan(dTdE(j, i));
          }
        }
        REQUIRE(!nansPresent);

        GIVEN("Databoxes filled for log rho and log sie") {
          eosDataOfRhoSie(matid, lRhoBounds, leBounds, P, T, bMod, dPdRho, dPdE, dTdRho,
                          dTdE, dEdRho, mask, Verbosity::Verbose);
          THEN("The shapes are correct") {
            REQUIRE(P.dim(1) == NUME);
            REQUIRE(P.dim(2) == NRHO);
            REQUIRE(T.dim(1) == NUME);
            REQUIRE(T.dim(2) == NRHO);
            REQUIRE(bMod.dim(1) == NUME);
            REQUIRE(bMod.dim(2) == NRHO);
            REQUIRE(dPdRho.dim(1) == NUME);
            REQUIRE(dPdRho.dim(2) == NRHO);
            REQUIRE(dPdE.dim(1) == NUME);
            REQUIRE(dPdE.dim(2) == NRHO);
            REQUIRE(dTdRho.dim(1) == NUME);
            REQUIRE(dTdRho.dim(2) == NRHO);
            REQUIRE(dTdE.dim(1) == NUME);
            REQUIRE(dTdE.dim(2) == NRHO);

            AND_THEN("The tables contain all finite values") {
              bool nansPresent = false;
              for (int j = 0; j < NRHO; j++) {
                for (int i = 0; i < NUME; i++) {
                  nansPresent |= std::isnan(P(j, i));
                  nansPresent |= std::isnan(T(j, i));
                  nansPresent |= std::isnan(bMod(j, i));
                  nansPresent |= std::isnan(dPdRho(j, i));
                  nansPresent |= std::isnan(dPdRho(j, i));
                  nansPresent |= std::isnan(dPdE(j, i));
                  nansPresent |= std::isnan(dTdRho(j, i));
                  nansPresent |= std::isnan(dTdE(j, i));
                }
              }
              REQUIRE(!nansPresent);
            }
          }
        }
      }
    }
  }
}

SCENARIO("Converting air to file using saveMaterial", "[saveMaterial]") {
  SesameMetadata m;
  eosGetMetadata(matid, m, Verbosity::Verbose);

  std::string name = m.name;
  Bounds lRhoBounds(m.rhoMin, m.rhoMax, m.numRho, true, 0.15);
  Bounds lTBounds(m.TMin, m.TMax, m.numT, true, 0.15);
  // Extremely conservative because extrapolation likely to go bad
  // Don't actually do this
  Bounds leBounds(m.sieMin, m.sieMax, m.numT, true, 0.5);

  GIVEN("HDF5 file can be written to") {
    hid_t file =
        H5Fcreate(SP5::defaultSesFileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    herr_t status = saveMaterial(file, m, lRhoBounds, lTBounds, leBounds, name);
    status += H5Fclose(file);
    THEN("HDF5 returns no errors") { REQUIRE(status == H5_SUCCESS); }
  }
}
