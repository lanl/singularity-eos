//------------------------------------------------------------------------------
// © 2026. Triad National Security, LLC. All rights reserved.  This
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

#ifndef SINGULARITY_UTILS_CONSTANTS_HPP_
#define SINGULARITY_UTILS_CONSTANTS_HPP_

namespace singularity {

// Enum for table splitting modes
// Used when converting EOSPAC tables that can be split into
// electron-only and ion-cold contributions
enum class TableSplit { Total = 0, ElectronOnly = 1, IonCold = 2 };

} // namespace singularity

#endif // SINGULARITY_UTILS_CONSTANTS_HPP_
