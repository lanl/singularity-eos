//------------------------------------------------------------------------------
// © 2025. Triad National Security, LLC. All rights reserved.  This
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

#include <limits>

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

PORTABLE_INLINE_FUNCTION Real SetRhoPMin(DataBox &P, DataBox &rho_at_pmin,
                                         const bool pmin_vapor_dome,
                                         const Real VAPOR_DPDR_THRESH,
                                         const Real lRhoOffset) {
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
      if ((j > 0) && pmin_vapor_dome) {
        Real dP = P(j, i) - P(j - 1, i);
        Real dr = from_log(lRs.x(j), lRhoOffset) - from_log(lRs.x(j - 1), lRhoOffset);
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
