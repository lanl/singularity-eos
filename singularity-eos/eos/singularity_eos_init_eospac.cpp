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

#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/singularity_eos.hpp>
#include <singularity-eos/eos/singularity_eos_init.hpp>

using namespace singularity;

#ifdef SINGULARITY_USE_EOSPAC
int init_sg_eospac(const int matindex, EOS *eos, const int id, double *const eospac_vals,
                   int const *const enabled, double *const vals) {

  using namespace EospacWrapper;
  assert(matindex >= 0);
  bool invert_at_setup = eospac_vals[0];
  double insert_data = eospac_vals[1];
  eospacMonotonicity monotonicity = static_cast<eospacMonotonicity>(eospac_vals[2]);
  bool apply_smoothing = eospac_vals[3];
  eospacSplit apply_splitting = static_cast<eospacSplit>(eospac_vals[4]);
  bool linear_interp = eospac_vals[5];

  EOS eos_ = 
      EOSPAC(id, invert_at_setup = invert_at_setup, insert_data = insert_data,
             monotonicity = monotonicity, apply_smoothing = apply_smoothing,
             apply_splitting = apply_splitting, linear_interp = linear_interp);
  apply_mods_and_dev_transfer(matindex, eos_, eos, enabled, vals);
  return 0;
}
int init_sg_eospac(const int matindex, EOS *eos, const int id,
                   double *const eospac_vals) {
  return init_sg_eospac(matindex, eos, id, eospac_vals, def_en, def_v);
}
#endif // SINGULARITY_USE_EOSPAC
