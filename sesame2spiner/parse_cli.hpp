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

#ifndef _SESAME2SPINER_PARSER_HPP_
#define _SESAME2SPINER_PARSER_HPP_

#include "io_eospac.hpp"
#include <string>
#include <vector>

const std::string DEFAULT_SAVENAME = "materials.sp5";
const std::string EXAMPLESTRING = R"(
# air.dat
# These are comments. 
# The "#" character must be at the beginning of a line.
# only matid is required. All others override defaults.
matid = 5030
name = air
# rho is in g/cm^3
rhomin = 1e-2
rhomax = 10
numrho = 64
# T is in Kelvin
Tmin = 252
Tmax = 1e4
numT = 32
# sie is in erg/g
siemin = 1e12
siemax = 1e16
numsie = 32


# titanium.dat
matid = 2961
name = titanium
# These set the number of grid poitns per decade
# for each variable. The default is 50 points.
numrho/decade = 30
numT/decade = 25
numSie/decade = 15


# steel.dat
matid=4272
rhomin = 1e-2
Tmin = 1
# These shrink lograithm of bounds
# by a fraction of the total interval <= 1
shrinklRhoBounds = 0.15
shrinklTBounds = 0.15
shrinkleBounds = 0.5
)";

void parseCLI(int argc, char *argv[], std::string &savename,
              std::vector<std::string> &filenames, bool &printMetadata,
              Verbosity &eospacWarn, std::string &helpMessage);

#endif // _SESAME2SPINER_PARSER_HPP_
