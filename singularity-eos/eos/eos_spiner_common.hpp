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

// ports-of-call
#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>

// spiner
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/sp5.hpp>
#include <spiner/spiner_types.hpp>

// implementation details below
// ======================================================================

// TODO: we're using log-linear interpolation, not log-log
// this may be suboptimal. We may want a way to do some variables
// in log-log. In particular, it might be good to do:
// pressure, energy, and bulk modulus in log-log.
// ~JMM

// replace lambdas with callable

namespace callable_interp {

static constexpr int NGRIDS = 3;
using Grid_t = Spiner::PiecewiseGrid1D<Real, NGRIDS>;
using DataBox = Spiner::DataBox<Real, Grid_t>;

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

#endif // SINGULARITY_USE_SPINER_WITH_HDF5
#endif // _SINGULARITY_EOS_EOS_EOS_SPINER_COMMON_HPP_
