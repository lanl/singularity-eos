//------------------------------------------------------------------------------
// © 2021-2026. Triad National Security, LLC. All rights reserved.  This
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

namespace sesame2spiner {
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

} // namespace sesame2spiner
#endif // SESAME2SPINER_PARSER_HPP_
