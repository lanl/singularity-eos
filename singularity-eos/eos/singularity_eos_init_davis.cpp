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

int init_sg_DavisProducts(const int matindex, EOS *eos, const double a, const double b,
                          const double k, const double n, const double vc,
                          const double pc, const double Cv, const double E0,
                          int const *const enabled, double *const vals) {
  EOS eos_ = DavisProducts(a, b, k, n, vc, pc, Cv, E0);
  apply_mods_and_dev_transfer(matindex, eos_, eos, enabled, vals);
  return 0;
}

int init_sg_DavisProducts(const int matindex, EOS *eos, const double a, const double b,
                          const double k, const double n, const double vc,
                          const double pc, const double Cv, const double E0) {
  return init_sg_DavisProducts(matindex, eos, a, b, k, n, vc, pc, Cv, E0, def_en, def_v);
}

int init_sg_DavisReactants(const int matindex, EOS *eos, const double rho0,
                           const double e0, const double P0, const double T0,
                           const double A, const double B, const double C,
                           const double G0, const double Z, const double alpha,
                           const double Cv0, int const *const enabled,
                           double *const vals) {
  EOS eos_ = DavisReactants(rho0, e0, P0, T0, A, B, C, G0, Z, alpha, Cv0);
  apply_mods_and_dev_transfer(matindex, eos_, eos, enabled, vals);
  return 0;
}

int init_sg_DavisReactants(const int matindex, EOS *eos, const double rho0,
                           const double e0, const double P0, const double T0,
                           const double A, const double B, const double C,
                           const double G0, const double Z, const double alpha,
                           const double Cv0) {
  return init_sg_DavisReactants(matindex, eos, rho0, e0, P0, T0, A, B, C, G0, Z, alpha,
                                Cv0, def_en, def_v);
}
