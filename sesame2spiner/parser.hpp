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

#ifndef SESAME2SPINER_PARSER_HPP_
#define SESAME2SPINER_PARSER_HPP_

#include <istream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// Parse a simple parameter file with
// "#" denoting comments.
class Params {
 public:
  Params() = default;
  Params(const std::string &input_file);
  Params(std::stringstream &input) { Parse(input); }

  bool Contains(const std::string &key) const { return params_.count(key); }
  template <typename T>
  T Get(const std::string &key) const;
  template <typename T>
  T Get(const std::string &key, T default_value) const {
    return Contains(key) ? Get<T>(key) : default_value;
  }
  void Print(std::ostream &s) const;
  void Clear() { params_.clear(); }
  bool Empty() const { return params_.size() == 0; }
  void Set(const std::string &key, const std::string &val) { params_[key] = val; }

 private:
  void Parse(std::istream &s);
  std::unordered_map<std::string, std::string> params_;
};

void AddMaterials(std::vector<Params> &params, std::vector<int> &matids,
                  const std::string &input_file);

#endif // SESAME2SPINER_PARSER_HPP_
