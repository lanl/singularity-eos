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

#ifndef _SINGULARITY_EOS_EOS_GET_SG_EOS_LAMBDAS_HPP_
#define _SINGULARITY_EOS_EOS_GET_SG_EOS_LAMBDAS_HPP_

#include <ports-of-call/portability.hpp>

#ifdef PORTABILITY_STRATEGY_KOKKOS

#define SG_GET_SG_EOS_INIT_LAMBDA_DECL \
PORTABLE_LAMBDA(const int i, const int tid, double &mass_sum, int &npte, \
                const Real t_mult, const Real s_mult, const Real p_mult) { \
  /* normalize mass fractions */ \
  /* first find the mass sum */ \
  /* also set idxs as the decrement of the eos offsets */ \
  /* to take into account 1 based indexing in fortran */ \
  for (int m = 0; m < nmat; ++m) { \
    mass_sum += frac_mass_v(i, m); \
    pte_idxs(tid, m) = eos_offsets_v(m) - 1; \
    frac_vol_v(i, m) = 0.0; \
  } \
  for (int m = 0; m < nmat; ++m) { \
    frac_mass_v(i, m) /= mass_sum; \
  } \
  /* set inputs */ \
  npte = 0; \
  for (int m = 0; m < nmat; ++m) { \
    if (frac_mass_v(i, m) > 1.e-12) { \
      pte_idxs(tid, npte) = eos_offsets_v(m) - 1; \
      pte_mats(tid, npte) = m; \
      npte += 1; \
    } \
    vfrac_pte(tid, m) = 0.0; \
    sie_pte(tid, m) = 0.0; \
    temp_pte(tid, m) = 0.0; \
    press_pte(tid, m) = 0.0; \
  } \
  for (int mp = 0; mp < npte; ++mp) { \
    const int m = pte_mats(tid, mp); \
    rho_pte(tid, mp) = npte / spvol_v(i) * frac_mass_v(i, m); \
    vfrac_pte(tid, mp) = 1.0 / npte; \
    temp_pte(tid, mp) = temp_v(i) * ev2k * t_mult; \
    press_pte(tid, mp) = press_v(i) * p_mult; \
    sie_pte(tid, mp) = sie_v(i) * frac_mass_v(i, m) * s_mult; \
  } \
}

#define SG_GET_SG_EOS_FINAL_LAMBDA_DECL \
PORTABLE_LAMBDA(const int i, const int tid, const int npte, \
                const Real mass_sum, const Real t_mult, const Real s_mult, \
		const Real p_mult, \
                singularity::mix_impl::CacheAccessor &cache) { \
  /* initialize averaged quantities to 0 */ \
  const bool do_t = t_mult == 1.0; \
  const bool do_s = s_mult == 1.0; \
  const bool do_p = p_mult == 1.0; \
  if (do_t) { \
    temp_v(i) = 0.0; \
  } \
  if (do_p) { \
    press_v(i) = 0.0; \
  } \
  if (do_s) { \
    sie_v(i) = 0.0; \
  } \
  bmod_v(i) = 0.0; \
  cv_v(i) = 0.0; \
  dpde_v(i) = 0.0; \
  /* material loop for averaging and assigning per mat quantities */ \
  for (int mp = 0; mp < npte; ++mp) { \
    const int m = pte_mats(tid, mp); \
    /* pressure contribution from material m */ \
    press_v(i) += press_pte(tid, mp) * vfrac_pte(tid, mp) * p_mult; \
    /* temperature contribution from material m */ \
    temp_v(i) += temp_pte(tid, mp) * vfrac_pte(tid, mp) * t_mult; \
    const Real ie_m = sie_pte(tid, mp) * frac_mass_v(i, m) * mass_sum; \
    /* sie contribution from material m */ \
    sie_v(i) += ie_m * s_mult; \
    /* assign per material specific internal energy */ \
    frac_ie_v(i, m) = ie_m; \
    /* assign volume fraction based on pte calculation */ \
    frac_vol_v(i, m) = vfrac_pte(tid, mp) * vol_v(i); \
    /* calculate bulk modulus for material m */ \
    const Real bmod_m = eos_v(pte_idxs(tid, mp)) \
                            .BulkModulusFromDensityTemperature( \
                                rho_pte(tid, mp), \
                                temp_pte(tid, mp), \
                                cache[mp]); \
    /* add bmod contribution from material m */ \
    bmod_v(i) += bmod_m * vfrac_pte(tid, mp); \
    /* calculate specific heat for material m */ \
    const Real cv_m = ev2k * eos_v(pte_idxs(tid, mp)) \
                                 .SpecificHeatFromDensityTemperature( \
                                     rho_pte(tid, mp), \
                                     temp_pte(tid, mp), \
                                     cache[mp]); \
    /* add mass weighted contribution specific heat for material m */ \
    cv_v(i) += cv_m * frac_mass_v(i, m); \
    /* calculate gruneisen parameter for material m */ \
    const Real dpde_m = eos_v(pte_idxs(tid, mp)) \
                            .GruneisenParamFromDensityTemperature( \
                                rho_pte(tid, mp), \
                                temp_pte(tid, mp), \
                                cache[mp]);\
    /* add gruneisen param contribution from material m */ \
    dpde_v(i) += dpde_m * vfrac_pte(tid, mp); \
    /* optionally assign per material quantities to per material arrays */ \
    if (do_frac_bmod) { \
      frac_bmod_v(i, m) = bmod_m; \
    } \
    if (do_frac_cv) { \
      frac_cv_v(i, m) = cv_m; \
    } \
    if (do_frac_dpde) { \
      frac_dpde_v(i, m) = dpde_m; \
    } \
  } \
  if (do_t) { \
    temp_v(i) /= ev2k; \
  } \
  if (do_s) { \
    sie_v(i) /= mass_sum; \
  } \
  /* reset mass fractions to original values if not normalized to 1 */ \
  for (int m = 0; m < nmat; ++m) { \
    frac_mass_v(i, m) *= mass_sum; \
  } \
}
#endif // PORTABILITY_STRATEGY_KOKKOS

#endif // _SINGULARITY_EOS_EOS_GET_SG_EOS_LAMBDAS_HPP_
