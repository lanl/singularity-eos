//======================================================================
// sesame2spiner tool for converting eospac to spiner
// Author: Jonah Miller (jonahm@lanl.gov)
// © 2021-2022. Triad National Security, LLC. All rights reserved.  This
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
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <hdf5.h>
#include <hdf5_hl.h>

#ifndef SPINER_USE_HDF
#error "HDF5 must be enabled"
#endif

#include <singularity-eos/base/sp5/singularity_eos_sp5.hpp>
#include <ports-of-call/portability.hpp>
#include <spiner/sp5.hpp>

#include "generate_files.hpp"
#include "io_eospac.hpp"
#include "parse_cli.hpp"

int main(int argc, char *argv[]) {
  std::vector<std::string> filenames;
  std::string savename, helpMessage;
  Verbosity eospacWarn = Verbosity::Quiet;
  bool printMetadata = false;
  herr_t status = H5_SUCCESS;

  parseCLI(argc, argv, savename, filenames, printMetadata, eospacWarn, helpMessage);

  std::cout << "sesame2spiner                            \n"
            << "-----------------------------------------\n"
            << "Author: Jonah Miller (jonahm@lanl.gov)   \n"
            << "-----------------------------------------\n"
            << std::endl;

  status = saveAllMaterials(savename, filenames, printMetadata, eospacWarn);

  std::cout << "Done." << std::endl;

  return (status == H5_SUCCESS) ? 0 : 1;
}
