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

#include <regex>
#include <string>
#include <sstream>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <fstream>

#include "nlohmann/json.hpp"
#include "parser.hpp"
#include "io_eospac.hpp"

std::string removeComments(const std::string& in) {
  using std::regex;
  using std::regex_replace;

  // first part is /* comment */. Second part is // comment
  // Source:
  // https://blog.ostermiller.org/find-comment
  regex comments("(/\\*([^*]|[\r\n]|(\\*+([^*/]|[\r\n])))*\\*+/)|(//(.*))");

  return regex_replace(in,comments,"");
}

nlohmann::json stringToJson(const std::string& in) {
  return nlohmann::json::parse(removeComments(in));
}

nlohmann::json fileToJson(const std::string& filename) {
  std::ifstream file(filename.c_str());
  std::stringstream buffer;
  if (!file.good()) {
    std::cerr << "Cannot read file " << filename << std::endl;
    std::exit(1);
  }
  buffer << file.rdbuf();
  return stringToJson(buffer.str());
}

void parseCLI(int argc, char* argv[],
	      std::string& filename, bool& printMetadata,
	      Verbosity& eospacWarn,
	      std::string& helpMessage) {

  int argCount = 2;

  std::stringstream helpStream;
  helpStream << "Usage: " << argv[0] << "[-p] [-w] [-h] <parameter file>\n\n"
	     << "\t <parameter file>: input file in json format\n"
	     << "\t-p:  print metadata associated with materials "
	     << "in parameter file\n"
       << "\t-v:  print eospac warnings\n"
       << "\t-vv: print debug information\n"
	     << "\t-w:  same as -v\n"
       << "\t-d:  same as -vv\n"
	     << "\t-h:  print this message\n"
	     << "\n"
	     << "Example JSON file:\n"
	     << EXAMPLESTRING
	     << "\n"
	     << std::endl;
  helpMessage = helpStream.str();

  if (argc < 2) {
    std::cerr << helpMessage << std::endl;
    std::exit(1);
  }
  if ( argc == 2 && std::strcmp(argv[1],"-h") == 0 ) {
    std::cout << helpMessage << std::endl;
    std::exit(0);
  }
  for (int i = 1; i < argc; i++) {
    if ( std::strcmp(argv[i],"-h") == 0 ) {
      std::cout << helpMessage << std::endl;
      std::exit(0);
    }
    else if ( std::strcmp(argv[i],"-p") == 0 ) {
      printMetadata = true;
      argCount++;
    } else if ( (std::strcmp(argv[i],"-w") == 0 ||
                 std::strcmp(argv[i],"-v") == 0 )
                && eospacWarn == Verbosity::Quiet) {
      eospacWarn = Verbosity::Verbose;
      argCount++;
    } else if ( (std::strcmp(argv[i],"-d") == 0 ||
                 std::strcmp(argv[i],"-vv") == 0 )
                && eospacWarn != Verbosity::Debug) {
      eospacWarn = Verbosity::Debug;
      argCount++;
    } else {
      filename = argv[i];
    }
  }
  if (argc != argCount) {
    std::cerr << helpMessage << std::endl;
    std::exit(1);
  }
}
