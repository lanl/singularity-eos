//------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------

#include <algorithm>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "parse_cli.hpp"
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
    if (line[0] == '#' || line.empty()) continue;
    auto delimiter_pos = line.find("=");
    auto name = line.substr(0, delimiter_pos);
    auto value = line.substr(delimiter_pos + 1);
    // make sure there's no trailing comment
    auto comment_pos = value.find("#");
    if (comment_pos != std::string::npos) {
      value.erase(comment_pos);
    }
    params_[name] = value;
  }
}

void Params::Print(std::ostream &s) const {
  s << "\nParams for " << params_.at("name") << "\n";
  for (const auto &pair : params_) {
    s << pair.first << " = " << pair.second << "\n";
  }
}

template <>
std::string Params::Get(const std::string &key) const {
  return params_.at(key);
}

template <>
int Params::Get(const std::string &key) const {
  return std::stoi(params_.at(key));
}

template <>
double Params::Get(const std::string &key) const {
  return std::stod(params_.at(key));
}

template <>
float Params::Get(const std::string &key) const {
  return std::stof(params_.at(key));
}

template <>
bool Params::Get(const std::string &key) const {
  std::string val = params_.at(key);
  if ((val == "true") || (val == "1")) {
    return true;
  } else if ((val == "false") || (val == "0")) {
    return false;
  } else {
    throw std::runtime_error("The value of " + key + " is not a boolean.\n");
  }
}

void add_param(const Params &p, std::vector<Params> &params, std::vector<int> &matids,
               const std::string &filename) {

  if (!p.Contains("matid")) {
    if (p.Contains("name")) {
      const auto &name = p.Get<std::string>("name");
      std::cerr << "Material " << name << " in file " << filename << "is missing matid.\n"
                << "Example input files:\n"
                << EXAMPLESTRING << std::endl;
      std::exit(1);
    } else {
      std::cerr << "A Material in file " << filename << "has no name.\n"
                << "Example input files:\n"
                << EXAMPLESTRING << std::endl;
      std::exit(1);
    }
  }
  matids.push_back(p.Get<int>("matid"));
  params.push_back(p);
}

void parse_file(std::vector<Params> &params, std::vector<int> &matids, std::istream &s,
                const std::string &filename) {
  std::string line;
  Params p;
  size_t line_num = 1;
  while (getline(s, line)) {
    line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
    if (line[0] == '#' || line.empty()) {
      line_num++;
      continue;
    }
    if (line[0] == '<') {
      auto term = line.find('>');
      if (term == std::string::npos) {
        std::cerr << "Missing closing > on line " << line_num << " of file " << filename
                  << "\n";
        std::exit(1);
      }

      auto name = line.substr(1, term - 1);
      if (!p.Empty()) add_param(p, params, matids, filename);

      p.Clear();
      p.Set("name", name);
    } else {
      auto delimiter_pos = line.find("=");
      auto name = line.substr(0, delimiter_pos);
      auto value = line.substr(delimiter_pos + 1);
      // make sure there's no trailing comment
      auto comment_pos = value.find("#");
      if (comment_pos != std::string::npos) {
        value.erase(comment_pos);
      }
      p.Set(name, value);
    }
    line_num++;
  }
  if (!p.Empty()) add_param(p, params, matids, filename);
}

void AddMaterials(std::vector<Params> &params, std::vector<int> &matids,
                  const std::string &input_file) {
  std::ifstream config_file(input_file);
  if (config_file.is_open()) {
    parse_file(params, matids, config_file, input_file);
  } else {
    throw std::runtime_error("Couldn't open config file " + input_file + "\n");
  }
}