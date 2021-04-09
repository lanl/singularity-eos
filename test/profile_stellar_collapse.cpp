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
  Tool for profiling the Stellar Collapse Table Reader
  Authors: Jonah Miller
 */

#ifdef SPINER_USE_HDF

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <chrono>

#include "hdf5.h"
#include "hdf5_hl.h"

#include "../utils/ports-of-call/portability.hpp"
#include "../utils/spiner/spiner_types.hpp"
#include "../utils/spiner/sp5.hpp"
#include "../utils/spiner/databox.hpp"
#include "../utils/spiner/interpolation.hpp"
#include "../eos/eos.hpp"

using namespace singularity;

using duration = std::chrono::microseconds;
constexpr char DIFFS_NAME[] = "diffs.sp5";
constexpr Real MP = 1.67262171e-24; // proton mass

int main(int argc, char* argv[]) {

  herr_t status = H5_SUCCESS;

  #ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
  #endif
  {
    if (argc < 4) {
      std::cerr << "Usage: " << argv[0] << "filename gamma nfine" << std::endl;
      std::exit(1);
    }
    std::string filename = argv[1];
    Real gamma = std::atof(argv[2]);
    int nfine = std::atoi(argv[3]);

    if (nfine < 1) {
      std::cerr << "We need at least one interpolation point" << std::endl;
      std::exit(1);
    }
    if (gamma <= 1) {
      std::cerr <"gamma - 1 must be a positive number" << std::endl;
      std::exit(1);
    }

    std::cout << "Profiling a stellar collapse table" << std::endl;
    std::cout << "\t...Creating an ideal gas equation of state with gamma = "
              << gamma << " to compare to."
              << std::endl;
    
  }
  return 0;
}

#endif // SPINER_USE_HDF
