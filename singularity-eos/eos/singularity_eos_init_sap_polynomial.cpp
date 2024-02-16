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

int init_sg_SAP_Polynomial(const int matindex, EOS *eos, const double rho0,
                           const double a0, const double a1, const double a2c,
                           const double a2e, const double a3, const double b0,
                           const double b1, const double b2c, const double b2e,
                           const double b3, int const *const enabled,
                           double *const vals) {
  EOS eos_ = SAP_Polynomial(rho0, a0, a1, a2c, a2e, a3, b0, b1, b2c, b2e, b3);
  apply_mods_and_dev_transfer(matindex, eos_, eos, enabled, vals);
  return 0;
}
