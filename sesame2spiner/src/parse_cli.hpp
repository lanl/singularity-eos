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

#ifndef _SESAME2SPINER_PARSER_HPP_
#define _SESAME2SPINER_PARSER_HPP_

#include <sesame2spiner/io_eospac.hpp>
#include <string>
#include <vector>

const std::string DEFAULT_SAVENAME = "materials.sp5";

using EospacWrapper::Verbosity;
void parseCLI(int argc, char *argv[], std::string &savename,
              std::vector<std::string> &filenames, bool &printMetadata,
              Verbosity &eospacWarn, std::string &helpMessage);

#endif // _SESAME2SPINER_PARSER_HPP_
