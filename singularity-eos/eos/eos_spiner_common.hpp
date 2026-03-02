//------------------------------------------------------------------------------
// © 2025-2026. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_EOS_EOS_SPINER_COMMON_HPP_
#define _SINGULARITY_EOS_EOS_EOS_SPINER_COMMON_HPP_

#ifdef SINGULARITY_USE_SPINER_WITH_HDF5

#include <filesystem>
#include <limits>
#include <type_traits>

#include <hdf5.h>
#include <hdf5_hl.h>

// ports-of-call
#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>

// spiner
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/sp5.hpp>
#include <spiner/spiner_types.hpp>

// singularity-eos
#include <singularity-eos/base/eos_error.hpp>
#include <singularity-eos/base/fast-math/logs.hpp>
#include <singularity-eos/base/robust_utils.hpp>

#define SPINER_EOS_VERBOSE (0)
#define SP_ROOT_FINDER (RootFinding1D::regula_falsi)

namespace singularity {
namespace spiner_common {

static constexpr int NGRIDS = 3;
using Grid_t = Spiner::PiecewiseGrid1D<Real, NGRIDS>;
using DataBox = Spiner::DataBox<Real, Grid_t>;

PORTABLE_FORCEINLINE_FUNCTION Real to_log(const Real x, const Real offset) {
  return FastMath::log10(std::abs(std::max(x, -offset) + offset) + robust::EPS());
}
PORTABLE_FORCEINLINE_FUNCTION Real from_log(const Real lx, const Real offset) {
  return FastMath::pow10(lx) - offset;
}

inline Real SetRhoPMin(DataBox &P, DataBox &rho_at_pmin, const bool pmin_vapor_dome,
                       const Real VAPOR_DPDR_THRESH, const Real lRhoOffset) {
  Real PMin = std::numeric_limits<Real>::max();
  const auto lTs = P.range(0);
  const auto lRs = P.range(1);
  const Real NT = lTs.nPoints();
  const Real NR = lRs.nPoints();
  rho_at_pmin.resize(NT);
  rho_at_pmin.setRange(0, lTs);
  for (int i = 0; i < NT; ++i) {
    Real PMin_at_T = std::numeric_limits<Real>::max();
    int jmax = 0;
    for (int j = 0; j < NR; ++j) {
      if (P(j, i) < PMin_at_T) {
        PMin_at_T = P(j, i);
        jmax = j;
      }
      // check gradient if excluding vapor dome
      if ((j > 1) && pmin_vapor_dome) {
        // JMM: Use width of 2 points here because width 1 could
        // accidentally catch hierarchical grid overlap points
        Real dP = P(j, i) - P(j - 2, i);
        Real dr = from_log(lRs.x(j), lRhoOffset) - from_log(lRs.x(j - 2), lRhoOffset);
        Real dpdr = robust::ratio(dP, dr);
        if (dpdr < VAPOR_DPDR_THRESH) {
          jmax = j;
        }
      }
    }
    if ((PMin_at_T > 0) && !pmin_vapor_dome) {
      rho_at_pmin(i) = 0;
      PMin_at_T = 0;
    } else {
      rho_at_pmin(i) = from_log(lRs.x(jmax), lRhoOffset);
    }
    PMin = std::min(PMin_at_T, PMin);
  }

  // IF we are looking for the edge of the vapor dome, try to capture
  // multiple Van der Waals loops by enforcing monotonicity of
  // rho_at_pmin with respect to T.
  if (pmin_vapor_dome) {
    for (int i = NT - 2; i >= 0; i--) {
      if (rho_at_pmin(i) < rho_at_pmin(i + 1)) {
        rho_at_pmin(i) = rho_at_pmin(i + 1);
      }
    }
  }

  return PMin;
}

PORTABLE_INLINE_FUNCTION void PrintRhoPMin(const DataBox &rho_at_pmin,
                                           const Real lTOffset) {
  const auto &range = rho_at_pmin.range(0);
  for (std::size_t i = 0; i < range.nPoints(); ++i) {
    const Real lT = range.x(i);
    const Real T = from_log(lT, lTOffset);
    const Real rho = rho_at_pmin(i);
    printf("%ld %.14e %.14e\n", i, T, rho);
  }
}

inline herr_t aborting_error_handler(hid_t stack, void *client_data) {
  H5Eprint2(stack, stderr);
  PORTABLE_ALWAYS_THROW_OR_ABORT("HDF5 error detected! Erroring out!");
  return -1;
}

inline hid_t h5_safe_fopen(const char *filename, hid_t mode, hid_t property) {
  if (!std::filesystem::exists(filename)) {
    std::string msg = "Failed to open file " + std::string(filename);
    EOS_ERROR(msg);
  }
  return H5Fopen(filename, mode, property);
}
inline void h5_safe_fclose(hid_t &file) {
  if (H5Fclose(file) != H5_SUCCESS) {
    EOS_ERROR("Error closing hdf5 file");
  }
}
inline void h5_safe_gclose(hid_t &grp) {
  if (H5Gclose(grp) != H5_SUCCESS) {
    EOS_ERROR("Error closing hdf5 group");
  }
}

inline hid_t h5_safe_gopen(hid_t &file, const char *name, hid_t property) {
  if (!H5Lexists(file, name, property)) {
    std::string msg = "Failed to open group " + std::string(name);
    EOS_ERROR(msg);
  }
  return H5Gopen(file, name, property);
}

inline char *h5_safe_read_attr_string(hid_t &file, const char *grp, const char *name,
                                      std::size_t &len) {
  // NOTE THIS FUNCTION ALLOCATES THE BUFFER STRING!!!
  hsize_t dims[1] = {0};
  H5T_class_t attr_type;
  if (H5LTget_attribute_info(file, grp, name, dims, &attr_type, &len) != H5_SUCCESS) {
    std::string msg = "Could not read the string " + std::string(name);
    EOS_ERROR(msg);
    return nullptr;
  }
  if (len <= 0) {
    std::string msg = "Length of string" + std::string(name) + " is 0";
    EOS_ERROR(msg);
    return nullptr;
  }
  if (attr_type != H5T_STRING) {
    std::string msg = "Attribute type is NOT a string for " + std::string(name);
    EOS_ERROR(msg);
    return nullptr;
  }
  len += 1;
  char *buffer = (char *)malloc(len);
  if (H5LTget_attribute_string(file, grp, name, buffer) != H5_SUCCESS) {
    std::string msg = "Error reading string attribute for " + std::string(name);
    free(buffer);
    EOS_ERROR(msg);
    return nullptr;
  }
  buffer[len - 1] = '\0';
  return buffer;
}

template <typename T>
inline void h5_safe_get_attribute(hid_t &file, const char *grp, const char *name,
                                  T *output, const bool optional = false) {
  herr_t status = H5_SUCCESS;
  if constexpr (std::is_same_v<T, int>) {
    status = H5LTget_attribute_int(file, grp, name, output);
  } else if constexpr (std::is_same_v<T, double>) {
    status = H5LTget_attribute_double(file, grp, name, output);
  } else if constexpr (std::is_same_v<T, char>) {
    status = H5LTget_attribute_char(file, grp, name, output);
  } else {
    EOS_ERROR("Unsupported attribute type in h5_safe_get_attribute");
  }

  if ((!optional) && (status != H5_SUCCESS)) {
    std::string msg =
        "Attribute " + std::string(name) + " in group " + std::string(grp) + " not found";
    EOS_ERROR(msg);
  }
  return;
}

} // namespace spiner_common

// TODO: we're using log-linear interpolation, not log-log
// this may be suboptimal. We may want a way to do some variables
// in log-log. In particular, it might be good to do:
// pressure, energy, and bulk modulus in log-log.
// ~JMM

// replace lambdas with callable

namespace callable_interp {

using namespace spiner_common;

class l_interp {
 private:
  const DataBox &field;
  const Real fixed;

 public:
  PORTABLE_INLINE_FUNCTION
  l_interp(const DataBox &field_, const Real fixed_) : field{field_}, fixed{fixed_} {}

  PORTABLE_INLINE_FUNCTION Real operator()(const Real x) const {
    return field.interpToReal(x, fixed);
  }
};

class r_interp {
 private:
  const DataBox &field;
  const Real fixed;

 public:
  PORTABLE_INLINE_FUNCTION
  r_interp(const DataBox &field_, const Real fixed_) : field{field_}, fixed{fixed_} {}

  PORTABLE_INLINE_FUNCTION Real operator()(const Real x) const {
    return field.interpToReal(fixed, x);
  }
};

class prod_interp_1d {
 private:
  const DataBox &field1, field2;
  const Real r;

 public:
  PORTABLE_INLINE_FUNCTION
  prod_interp_1d(const DataBox &field1_, const DataBox &field2_, const Real r_)
      : field1{field1_}, field2{field2_}, r{r_} {}
  PORTABLE_INLINE_FUNCTION Real operator()(const Real x) const {
    return field1.interpToReal(x) * field2.interpToReal(x) * r * x;
  }
};

class interp {
 private:
  const DataBox &field;

 public:
  PORTABLE_INLINE_FUNCTION
  interp(const DataBox &field_) : field(field_) {}
  PORTABLE_INLINE_FUNCTION Real operator()(const Real x) const {
    return field.interpToReal(x);
  }
};
} // namespace callable_interp
} // namespace singularity

#endif // SINGULARITY_USE_SPINER_WITH_HDF5
#endif // _SINGULARITY_EOS_EOS_EOS_SPINER_COMMON_HPP_
