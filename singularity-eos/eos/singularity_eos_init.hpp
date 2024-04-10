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

#ifndef _SINGULARITY_EOS_EOS_SINGULARITY_EOS_INIT_
#define _SINGULARITY_EOS_EOS_SINGULARITY_EOS_INIT_

#include <cassert>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_builder.hpp>
#include <singularity-eos/eos/singularity_eos.hpp>

using namespace singularity;

static int def_en[4] = {0, 0, 0, 0};
static double def_v[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

// function that applies scale, shift, and ramp in the correct order
// and does the correct p-alpha conversion if necessary
inline void apply_modifiers(EOS& eos, bool const apply_scale,
			    bool const apply_shift, bool const apply_ramp,
			    bool const convert_pa_to_bl, Real const& scale,
			    Real const& shift, Real& r0, Real& a, Real& b,
			    Real& c) {
  if (apply_shift) {
    eos = eos.Modify<ShiftedEOS>(shift);
  }
  if (apply_scale) {
    eos = eos.Modify<ScaledEOS>(scale);
  }
  if (apply_ramp) {
    if (convert_pa_to_bl) {
      // if here the variable names passed in are p-alpha params and
      // need to be converted
      Real v1 = r0, v2 = a, v3 = b;
      singularity::pAlpha2BilinearRampParams(eos, v1, v2, v3, r0, a, b, c);
    }
    eos = eos.Modify<BilinearRampEOS>(r0, a, b, c);
  }
  return;
}

// calls the apply modifiers function given the expected arrays
inline void apply_modifiers_from_arrs(EOS& eos, int const *const enabled,
			              double *const vals) {
  apply_modifiers(eos, enabled[0] == 1, enabled[1] == 1,
		  enabled[2] == 1 || enabled[3] == 1, enabled[3] == 1,
		  vals[0], vals[1], vals[2], vals[3], vals[4], vals[5]);
  return;
}

// applies modifiers given arrays and does the device transfer
// and index assertion
inline void apply_mods_and_dev_transfer(int const matidx, EOS eos, EOS* eoss,
					int const *const enabled,
					double *const vals) {
  assert(matidx >= 0);
  apply_modifiers_from_arrs(eos, enabled, vals);
  eoss[matidx] = eos.GetOnDevice();
}

#endif
