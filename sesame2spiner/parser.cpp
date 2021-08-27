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

#include <algorithm>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "parser.hpp"

Params::Params(const std::string &input_file) {
  std::ifstream config_file(input_file);
  if (config_file.is_open()) {
    Parse(config_file);
  } else {
    throw std::runtime_error("Couldn't open config file " + input_file + "\n");
  }
}

void Params::Parse(std::istream &s) {
  std::string line;
  while (getline(s, line)) {
    line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
    if (line[0] == '#' || line.empty())
      continue;
    auto delimiter_pos = line.find("=");
    auto name = line.substr(0, delimiter_pos);
    auto value = line.substr(delimiter_pos + 1);
    params_[name] = value;
  }
}

template <> std::string Params::Get(const std::string &key) const {
  return params_.at(key);
}

template <> int Params::Get(const std::string &key) const {
  return std::stoi(params_.at(key));
}

template <> double Params::Get(const std::string &key) const {
  return std::stod(params_.at(key));
}

template <> float Params::Get(const std::string &key) const {
  return std::stof(params_.at(key));
}

template <> bool Params::Get(const std::string &key) const {
  std::string val = params_.at(key);
  if ((val == "true") || (val == "1")) {
    return true;
  } else if ((val == "false") || (val == "0")) {
    return false;
  } else {
    throw std::runtime_error("The value of " + key + " is not a boolean.\n");
  }
}
