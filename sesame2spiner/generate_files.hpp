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

#ifndef _SESAME2SPINER_GENERATE_FILES_HPP_
#define _SESAME2SPINER_GENERATE_FILES_HPP_

#include <string>
#include <vector>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <eospac-wrapper/eospac_wrapper.hpp>

#include "io_eospac.hpp"
#include "parser.hpp"

using namespace EospacWrapper;

constexpr int PPD_DEFAULT_RHO = 350;
constexpr int PPD_DEFAULT_T = 100;
constexpr Real STRICTLY_POS_MIN_RHO = 1e-8;
constexpr Real STRICTLY_POS_MIN_T = 1e-2;
constexpr Real COARSE_FACTOR_DEFAULT_RHO_LO = 3;
constexpr Real COARSE_FACTOR_DEFAULT_RHO_HI = 5;
constexpr Real COARSE_FACTOR_DEFAULT_T = 1.5;
constexpr Real RHO_FINE_DIAMETER_DEFAULT = 1.5;
constexpr Real T_SPLIT_POINT_DEFAULT = 1e4;

herr_t saveMaterial(hid_t loc, const SesameMetadata &metadata, const Bounds &lRhoBounds,
                    const Bounds &lTBounds, const Bounds &leBounds,
                    const std::string &name, Verbosity eospacWarn = Verbosity::Quiet);

herr_t saveAllMaterials(const std::string &savename,
                        const std::vector<std::string> &filenames, bool printMetadata,
                        Verbosity eospacWarn);

void getMatBounds(int i, int matid, const SesameMetadata &metadata, const Params &params,
                  Bounds &lRhoBounds, Bounds &lTBounds, Bounds &leBounds);

bool checkValInMatBounds(int matid, const std::string &name, Real val, Real vmin,
                         Real vmax);

int getNumPointsFromPPD(Real min, Real max, int ppd);

#endif // _SESAME2SPINER_GENERATE_FILES_HPP_
