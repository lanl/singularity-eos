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

int init_sg_JWL(const int matindex, EOS *eos, const double A, const double B,
                const double R1, const double R2, const double w, const double rho0,
                const double Cv, int const *const enabled, double *const vals) {
  EOS eos_ = JWL(A, B, R1, R2, w, rho0, Cv);
  apply_mods_and_dev_transfer(matindex, eos_, eos, enabled, vals);
  return 0;
}

int init_sg_JWL(const int matindex, EOS *eos, const double A, const double B,
                const double R1, const double R2, const double w, const double rho0,
                const double Cv) {
  return init_sg_JWL(matindex, eos, A, B, R1, R2, w, rho0, Cv, def_en, def_v);
}