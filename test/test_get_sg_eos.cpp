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

#include <iostream>
#include <memory>
#include <stdlib.h>

#include <ports-of-call/portability.hpp>
#include <pte_test_utils.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/singularity_eos.hpp>

using namespace singularity;

#ifdef PORTABILITY_STRATEGY_KOKKOS
// TODO DAH: when the get_sg_eos function is moved out of sg,
// this function will have to change how the PTE solutions
// are obtained.
int run_sg_get_eos_tests() {
  int nfails = 0;
  // initialize inputs outputs
  static constexpr const double ev2k = 1.160451930280894026e4;
  EOS eoss[NMAT];
  set_eos(eoss);
  // set volume fractions to equipartition
  // all variables needed to check the pte solve
  Real mfrac[NMAT];
  Real P_true = 5.e10, T_true = 800.0;
  Real T_true_ev = T_true / ev2k;
  Real mass_tot = 0.0, sie_tot_true = 0.0, v_true = 1.0, rho_tot, spvol;
  Real pmax, bmod, dpde, cv;
  mfrac[0] = 0.67;
  mfrac[1] = 0.14;
  mfrac[2] = 0.19;
  for (int i = 0; i < NMAT; ++i) {
    mass_tot += mfrac[i];
  }
  rho_tot = mass_tot / v_true;
  spvol = 1.0 / rho_tot;
  int cell_offset = 1;
  int eos_offset[NMAT];
  for (int m = 0; m < NMAT; ++m) {
    eos_offset[m] = m + 1;
  }
  // do a rho/e(P,T) solve
  Real vfrac_true[NMAT], ie_true[NMAT];
  get_sg_eos(NMAT, 1, 1, -1, eos_offset, eoss, &cell_offset, &P_true, &pmax, &v_true,
             &spvol, &sie_tot_true, &T_true_ev, &bmod, &dpde, &cv, mfrac, vfrac_true,
             ie_true, nullptr, nullptr, nullptr);
  Real sie_tot_check = 0.0;
  for (int m = 0; m < NMAT; ++m) {
    const Real r_m = mfrac[m] / vfrac_true[m];
    const Real sie_m = ie_true[m] / mfrac[m];
    const Real p_eos = eoss[m].PressureFromDensityInternalEnergy(r_m, sie_m);
    const Real p_resid = std::abs(P_true - p_eos) / P_true;
    const Real t_eos = eoss[m].TemperatureFromDensityInternalEnergy(r_m, sie_m);
    const Real t_resid = std::abs(T_true - t_eos) / T_true;
    // check internal consitency of P-T input function
    if (t_resid > 1.e-5 || p_resid > 1.e-5) {
      printf("P-T: t_resid: %e | p_resid: %e\n", t_resid, p_resid);
      nfails += 1;
    }
    sie_tot_check += ie_true[m];
  }
  sie_tot_check /= mass_tot;
  // further check of internal consintency of P-T input function
  if (std::abs(sie_tot_check - sie_tot_true) / std::abs(sie_tot_true) > 1.e-5) {
    printf("P-T: sie_tot_true: %e | sie_tot_check: %e\n", sie_tot_true, sie_tot_check);
    nfails += 1;
  }
  // obtain converged and consistent PTE solution
  // do rho-T input solve
  Real p_check, vfrac_check[NMAT], ie_check[NMAT];
  get_sg_eos(NMAT, 1, 1, -3, eos_offset, eoss, &cell_offset, &p_check, &pmax, &v_true,
             &spvol, &sie_tot_check, &T_true_ev, &bmod, &dpde, &cv, mfrac, vfrac_check,
             ie_check, nullptr, nullptr, nullptr);
  // check output pressure and sie, indicate failure if relative err is too large
  if (std::abs(P_true - p_check) / std::abs(P_true) > 1.e-5 ||
      std::abs(sie_tot_true - sie_tot_check) / std::abs(sie_tot_true) > 1.e-5) {
    printf("r-T: p_true: %e | p_check: %e\n", P_true, p_check);
    printf("r-T: sie_tot_true: %e | sie_tot_check: %e\n", sie_tot_true, sie_tot_check);
    nfails += 1;
  }
  Real max_vfrac_resid = 0.0;
  Real max_sie_resid = 0.0;
  for (int m = 0; m < NMAT; ++m) {
    max_vfrac_resid = std::max(max_vfrac_resid, std::abs(vfrac_true[m] - vfrac_check[m]) /
                                                    std::abs(vfrac_true[m]));
    max_sie_resid = std::max(max_sie_resid,
                             std::abs(ie_true[m] - ie_check[m]) / std::abs(ie_true[m]));
  }
  if (max_vfrac_resid > 1.e-5 || max_sie_resid > 1.e-5) {
    printf("r-T: vr: %e | sr: %e\n", max_vfrac_resid, max_sie_resid);
    nfails += 1;
  }
  // do rho-P input solve
  Real t_check;
  get_sg_eos(NMAT, 1, 1, -2, eos_offset, eoss, &cell_offset, &P_true, &pmax, &v_true,
             &spvol, &sie_tot_check, &t_check, &bmod, &dpde, &cv, mfrac, vfrac_check,
             ie_check, nullptr, nullptr, nullptr);
  // check output temperature and sie, indicate failure if relative err is too large
  if (std::abs(T_true_ev - t_check) / std::abs(T_true_ev) > 1.e-5 ||
      std::abs(sie_tot_true - sie_tot_check) / std::abs(sie_tot_true) > 1.e-5) {
    printf("p-T: t_true: %e | t_check: %e\n", T_true_ev, t_check);
    printf("p-T: sie_tot_true: %e | sie_tot_check: %e\n", sie_tot_true, sie_tot_check);
    nfails += 1;
  }
  max_vfrac_resid = 0.0;
  max_sie_resid = 0.0;
  for (int m = 0; m < NMAT; ++m) {
    max_vfrac_resid = std::max(max_vfrac_resid, std::abs(vfrac_true[m] - vfrac_check[m]) /
                                                    std::abs(vfrac_true[m]));
    max_sie_resid = std::max(max_sie_resid,
                             std::abs(ie_true[m] - ie_check[m]) / std::abs(ie_true[m]));
  }
  if (max_vfrac_resid > 1.e-5 || max_sie_resid > 1.e-5) {
    printf("p-T: vr: %e | sr: %e\n", max_vfrac_resid, max_sie_resid);
    nfails += 1;
  }
  return nfails;
}
#endif

int main(int argc, char *argv[]) {

  int nfails_get_sg_eos = 1;
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
  {
    // run the various solver tests:
    // since this uses the sg_get_eos, only run
    // if kokkos enable since that function requires it
    // to run the solvers
    nfails_get_sg_eos = run_sg_get_eos_tests();
    if (nfails_get_sg_eos > 0) {
      printf("nfails of fixed T/P solvers = %i\n", nfails_get_sg_eos);
    }
  }
  Kokkos::finalize();
#endif // PORTABILITY_STRATEGY_KOKKOS
  return nfails_get_sg_eos;
}
