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

#ifndef _SESAME2SPINER_IO_EOSPAC_HPP_
#define _SESAME2SPINER_IO_EOSPAC_HPP_

#include <string>
#include <vector>

#include <eos_Interface.h> // eospac API

#ifndef SINGULARITY_USE_SPINER_WITH_HDF5
#error "HDF5 must be enabled"
#endif

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/fast-math/logs.hpp>
#include <singularity-eos/base/spiner_table_utils.hpp>
#include <spiner/databox.hpp>

#include <eospac-wrapper/eospac_wrapper.hpp>

using EospacWrapper::Verbosity;
constexpr int NGRIDS = 3;
using singularity::TableSplit;
using Bounds = singularity::table_utils::Bounds<NGRIDS>;
using Grid_t = Spiner::PiecewiseGrid1D<Real, NGRIDS>;
using DataBox = Spiner::DataBox<Real, Grid_t>;

void eosDataOfRhoSie(int matid, const TableSplit split, const Bounds &lRhoBounds,
                     const Bounds &leBounds, DataBox &P, DataBox &T, DataBox &bMods,
                     DataBox &dPdRho, DataBox &dPdE, DataBox &dTdRho, DataBox &dTdE,
                     DataBox &dEdRho, DataBox &mask,
                     Verbosity eospacWarn = Verbosity::Quiet);
inline void eosDataOfRhoSie(int matid, const Bounds &lRhoBounds, const Bounds &leBounds,
                            DataBox &P, DataBox &T, DataBox &bMods, DataBox &dPdRho,
                            DataBox &dPdE, DataBox &dTdRho, DataBox &dTdE,
                            DataBox &dEdRho, DataBox &mask,
                            Verbosity eospacWarn = Verbosity::Quiet) {
  eosDataOfRhoSie(matid, TableSplit::Total, lRhoBounds, leBounds, P, T, bMods, dPdRho,
                  dPdE, dTdRho, dTdE, dEdRho, mask, eospacWarn);
}

void eosDataOfRhoT(int matid, const TableSplit split, const Bounds &lRhoBounds,
                   const Bounds &lTBounds, DataBox &Ps, DataBox &sies, DataBox &bMods,
                   DataBox &dPdRho, DataBox &dPdE, DataBox &dTdRho, DataBox &dTdE,
                   DataBox &dEdRho, DataBox &dEdT, DataBox &mask,
                   Verbosity eospacWarn = Verbosity::Quiet);
inline void eosDataOfRhoT(int matid, const Bounds &lRhoBounds, const Bounds &lTBounds,
                          DataBox &Ps, DataBox &sies, DataBox &bMods, DataBox &dPdRho,
                          DataBox &dPdE, DataBox &dTdRho, DataBox &dTdE, DataBox &dEdRho,
                          DataBox &dEdT, DataBox &mask,
                          Verbosity eospacWarn = Verbosity::Quiet) {
  eosDataOfRhoT(matid, TableSplit::Total, lRhoBounds, lTBounds, Ps, sies, bMods, dPdRho,
                dPdE, dTdRho, dTdE, dEdRho, dEdT, mask, eospacWarn);
}

void eosColdCurves(int matid, const Bounds &lRhoBounds, DataBox &Ps, DataBox &sies,
                   DataBox &dPdRho, DataBox &dEdRho, DataBox &bMod, DataBox &mask,
                   Verbosity eospacWarn = Verbosity::Quiet);

void eosColdCurveMask(int matid, const Bounds &lRhoBounds, const int numSie,
                      const DataBox &sieColdCurve, DataBox &mask,
                      Verbosity eospacWarn = Verbosity::Quiet);

void makeInterpPoints(std::vector<EOS_REAL> &v, const Bounds &b);

namespace impl {
template <typename T>
T select(TableSplit split, T a, T b, T c) {
  if (split == TableSplit::Total) return a;
  if (split == TableSplit::ElectronOnly)
    return b;
  else // if (split == TableSplit::IonCold)
    return c;
}
void modifyNames(TableSplit split, std::vector<std::string> &names);
} // namespace impl

#endif // _SESAME2SPINER_IO_EOSPAC_HPP_
