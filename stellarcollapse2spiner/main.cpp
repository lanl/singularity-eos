//======================================================================
// stellarcollapse2sapiner tool for converting stellar collapse tables
// to spiner
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

#include <iostream>
#include <string>

#include <singularity-eos/eos/eos.hpp>

using singularity::StellarCollapse;

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cout << "Converts a Stellar Collapse EOS file into an SP5 file.\n"
              << "Performs relevant data cleanup for use in fluid codes.\n"
              << "Usage:\n"
              << argv[0] << " input_filename ouutput_filename" << std::endl;
    return 1;
  }

  const std::string input_name = argv[1];
  const std::string output_name = argv[2];
  StellarCollapse eos(input_name, false, true);
  eos.Save(output_name);

  return 0;
}
