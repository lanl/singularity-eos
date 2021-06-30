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

#ifndef _PARSER_HPP_
#define _PARSER_HPP_

#include <string>
#include <nlohmann/json.hpp>
#include "io_eospac.hpp"

const std::string EXAMPLESTRING = R"(
{
    "savename"  : "materials.sp5",
    "materials" : [
  { // only matid is required. All others override defaults.
      "matid"  : 5030,
      "name"   : "air",
      "rhomin" : 1e-2, // g/cm^3
      "rhomax" : 10,
      "numrho" : 64,
      "Tmin"   : 252, // kelvin
      "Tmax"   : 1e4,
      "numT"   : 32,
      "siemin" : 1e12, // erg/g
      "siemax" : 1e16, // erg/g
      "numsie" : 32
  },
  {
      "matid" : 2961,
      "name"  : "titanium",
      /* These set the number of grid points per decade
         for each variable. The default is 20 points
         per decade.
      */
      "numrho/decade" : 30,
      "numT/decade"   : 25,
      "numSie/decade" : 15
  },
  {
      "matid" : 4272, // steel
      "rhomin" : 1e-2, // g/cm^3
      "Tmin"   : 1,    // kelvin
      /* These shrink logarithm of bounds
         by a fraction of the total interval <= 1
      */
      "shrinklRhoBounds" : 0.15,
      "shrinklTBounds"   : 0.15,
      "shrinkleBounds"   : 0.5
  }
    ]
}
)";

std::string removeComments(const std::string& in);
nlohmann::json stringToJson(const std::string& in);
nlohmann::json fileToJson(const std::string& filename);
void parseCLI(int argc, char* argv[],
	      std::string& filename, bool& printMetadata,
	      Verbosity& eospacWarn,
	      std::string& helpMessage);

#endif // _PARSER_HPP_
