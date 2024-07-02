//------------------------------------------------------------------------------
// © 2021-2024. Triad National Security, LLC. All rights reserved.  This
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
#include <singularity-eos/base/error_utils.hpp>
#include <singularity-eos/eos/get_sg_eos.hpp>

#include <cmath>

#ifdef PORTABILITY_STRATEGY_KOKKOS

namespace singularity {

struct init_functor {
 private:
  dev_frac_v frac_mass_v;
  ScratchV<int> pte_idxs;
  indirection_v eos_offsets_v;
  dev_frac_v frac_vol_v;
  dev_frac_v frac_ie_v;
  ScratchV<int> pte_mats;
  ScratchV<double> vfrac_pte;
  ScratchV<double> sie_pte;
  ScratchV<double> temp_pte;
  ScratchV<double> press_pte;
  ScratchV<double> rho_pte;
  dev_v spvol_v;
  dev_v temp_v;
  dev_v press_v;
  dev_v sie_v;
  int nmat;
  double mass_frac_cutoff;

 public:
  init_functor() = default;
  init_functor(dev_frac_v &frac_mass_v_, ScratchV<int> &pte_idxs_,
               indirection_v &eos_offsets_v_, dev_frac_v &frac_vol_v_,
               dev_frac_v frac_ie_v_, ScratchV<int> &pte_mats_,
               ScratchV<double> &vfrac_pte_, ScratchV<double> &sie_pte_,
               ScratchV<double> &temp_pte_, ScratchV<double> &press_pte_,
               ScratchV<double> &rho_pte_, dev_v &spvol_v_, dev_v &temp_v_,
               dev_v &press_v_, dev_v &sie_v_, int &nmat, double &mass_frac_cutoff_)
      : frac_mass_v{frac_mass_v_}, pte_idxs{pte_idxs_}, eos_offsets_v{eos_offsets_v_},
        frac_vol_v{frac_vol_v_}, frac_ie_v{frac_ie_v_}, pte_mats{pte_mats_},
        vfrac_pte{vfrac_pte_}, sie_pte{sie_pte_}, temp_pte{temp_pte_},
        press_pte{press_pte_}, rho_pte{rho_pte_}, spvol_v{spvol_v_}, temp_v{temp_v_},
        press_v{press_v_}, sie_v{sie_v_}, nmat{nmat}, mass_frac_cutoff{
                                                          mass_frac_cutoff_} {}

  PORTABLE_INLINE_FUNCTION
  void operator()(const int i, const int tid, double &mass_sum, int &npte,
                  double &vfrac_sum, const Real t_mult, const Real s_mult,
                  const Real p_mult) const {
    /* normalize mass fractions */
    /* first find the mass sum */
    /* also set idxs as the decrement of the eos offsets */
    /* to take into account 1 based indexing in fortran */
    for (int m = 0; m < nmat; ++m) {
      mass_sum += frac_mass_v(i, m);
      pte_idxs(tid, m) = eos_offsets_v(m) - 1;
      // Assume non-participating materials are all at cell density. PTE volumes
      // will be overwritten
      frac_vol_v(i, m) = frac_mass_v(i, m);
      // Initialize all materials to the cell sie (if provided). PTE sie values
      // will be overwritten. This also means that the sie for the PTE
      // materials is the same as the cell PTE after removing small materials
      frac_ie_v(i, m) = sie_v(i) * frac_mass_v(i, m) * s_mult;
    }
    for (int m = 0; m < nmat; ++m) {
      frac_mass_v(i, m) /= mass_sum;
    }
    check_all_vals(i);
    // count the number of participating materials
    npte = 0;
    for (int m = 0; m < nmat; ++m) {
      if (frac_mass_v(i, m) > mass_frac_cutoff) {
        // participating materials are those with non-negligible mass fractions
        pte_idxs(tid, npte) = eos_offsets_v(m) - 1;
        pte_mats(tid, npte) = m;
        npte += 1;
      }
      // zero the PTE solver material inputs
      vfrac_pte(tid, m) = 0.0;
      sie_pte(tid, m) = 0.0;
      temp_pte(tid, m) = 0.0;
      press_pte(tid, m) = 0.0;
    }
    // Populate the inputs with consistent values.
    // NOTE: the volume fractions and densities need to be consistent with the
    // total specific volume since they are used to calculate internal
    // quantities for the PTE solver
    vfrac_sum = 0.;
    for (int mp = 0; mp < npte; ++mp) {
      const int m = pte_mats(tid, mp);
      // Need to guess volume fractions
      vfrac_pte(tid, mp) = frac_mass_v(i, m);
      vfrac_sum += vfrac_pte(tid, mp);
      // Calculate densities to be consistent with these volume fractions
      rho_pte(tid, mp) = frac_mass_v(i, m) / spvol_v(i) / vfrac_pte(tid, mp);
      temp_pte(tid, mp) = temp_v(i) * ev2k * t_mult;
      press_pte(tid, mp) = press_v(i) * p_mult;
      sie_pte(tid, mp) = sie_v(i) * frac_mass_v(i, m) * s_mult;
    }
    PORTABLE_REQUIRE(vfrac_sum <= 1.0, "Volume fraction sum is greater than 1");
    return;
  }

 private:
  // Only run these checks when built with debug flags
  PORTABLE_FORCEINLINE_FUNCTION
  void check_all_vals(int const i) const {
#ifndef NDEBUG
    bool any_bad_vals = false;
    for (auto m = 0; m < nmat; ++m) {
      any_bad_vals =
          any_bad_vals || error_utils::bad_value(frac_mass_v(i, m), "frac_mass");
    }
    any_bad_vals = any_bad_vals || error_utils::bad_value(press_v(i), "pres");
    any_bad_vals = any_bad_vals || error_utils::bad_value(sie_v(i), "sie");
    any_bad_vals = any_bad_vals || error_utils::bad_value(spvol_v(i), "spvol");

    if (any_bad_vals) {
      using PortsOfCall::printf;
      printf("### Bad Value Output state:\n");
      printf("  ~~Bulk state~~\n");
      printf("    Pressure        : % #24.15g microbar\n", press_v(i));
      printf("    Specific IE     : % #24.15g erg/g\n", sie_v(i));
      printf("    Specific Volume : % #24.15g cm^3/g\n", spvol_v(i));
      printf("  ~~Mass fractions~~\n");
      printf("   mat:");
      printf(" %24s", "Mass fraction");
      printf("\n");
      for (auto m = 0; m < nmat; ++m) {
        printf("   %3i:", m + 1);
        printf(" %24.15g", frac_mass_v(i, m));
        printf("\n");
      }
      PORTABLE_ALWAYS_ABORT(
          "Bad values INPUT to singularity-eos interface. See output for details");
    }
  }
#endif // #ifndef NDEBUG
};

struct final_functor {
 private:
  dev_v spvol_v;
  dev_v temp_v;
  dev_v press_v;
  dev_v sie_v;
  dev_v bmod_v;
  dev_v cv_v;
  dev_v dpde_v;
  ScratchV<int> pte_mats;
  ScratchV<double> press_pte;
  ScratchV<double> vfrac_pte;
  ScratchV<double> temp_pte;
  ScratchV<double> sie_pte;
  dev_frac_v frac_mass_v;
  dev_frac_v frac_ie_v;
  dev_frac_v frac_vol_v;
  dev_v vol_v;
  Kokkos::View<EOS *, Llft> eos_v;
  ScratchV<int> pte_idxs;
  ScratchV<double> rho_pte;
  dev_frac_v frac_bmod_v;
  dev_frac_v frac_cv_v;
  dev_frac_v frac_dpde_v;
  int nmat;
  bool do_frac_bmod;
  bool do_frac_cv;
  bool do_frac_dpde;

 public:
  final_functor(dev_v &spvol_v_, dev_v &temp_v_, dev_v &press_v_, dev_v &sie_v_,
                dev_v &bmod_v_, dev_v &cv_v_, dev_v &dpde_v_, ScratchV<int> &pte_mats_,
                ScratchV<double> &press_pte_, ScratchV<double> &vfrac_pte_,
                ScratchV<double> &temp_pte_, ScratchV<double> &sie_pte_,
                dev_frac_v &frac_mass_v_, dev_frac_v &frac_ie_v_, dev_frac_v &frac_vol_v_,
                dev_v &vol_v_, Kokkos::View<EOS *, Llft> &eos_v_,
                ScratchV<int> &pte_idxs_, ScratchV<double> &rho_pte_,
                dev_frac_v &frac_bmod_v_, dev_frac_v &frac_cv_v_,
                dev_frac_v &frac_dpde_v_, int &nmat_, bool do_frac_bmod_,
                bool do_frac_cv_, bool do_frac_dpde_)
      : spvol_v{spvol_v_}, temp_v{temp_v_}, press_v{press_v_}, sie_v{sie_v_},
        bmod_v{bmod_v_}, cv_v{cv_v_}, dpde_v{dpde_v_}, pte_mats{pte_mats_},
        press_pte{press_pte_}, vfrac_pte{vfrac_pte_}, temp_pte{temp_pte_},
        sie_pte{sie_pte_}, frac_mass_v{frac_mass_v_}, frac_ie_v{frac_ie_v_},
        frac_vol_v{frac_vol_v_}, vol_v{vol_v_}, eos_v{eos_v_}, pte_idxs{pte_idxs_},
        rho_pte{rho_pte_}, frac_bmod_v{frac_bmod_v_}, frac_cv_v{frac_cv_v_},
        frac_dpde_v{frac_dpde_v_}, nmat{nmat_}, do_frac_bmod{do_frac_bmod_},
        do_frac_cv{do_frac_cv_}, do_frac_dpde{do_frac_dpde_} {}

 public:
  PORTABLE_INLINE_FUNCTION
  void operator()(const int i, const int tid, const int npte, const Real mass_sum,
                  const Real t_mult, const Real s_mult, const Real p_mult,
                  const bool pte_converged,
                  singularity::mix_impl::CacheAccessor &cache) const {
    /* initialize averaged quantities to 0 */
    const bool do_t = t_mult == 1.0;
    const bool do_s = s_mult == 1.0;
    const bool do_p = p_mult == 1.0;
    if (do_t) {
      temp_v(i) = 0.0;
    }
    if (do_p) {
      press_v(i) = 0.0;
    }
    if (do_s) {
      sie_v(i) = 0.0;
    }
    bmod_v(i) = 0.0;
    cv_v(i) = 0.0;
    dpde_v(i) = 0.0;
    /* material loop for averaging and assigning per mat quantities */
    for (int mp = 0; mp < npte; ++mp) {
      const int m = pte_mats(tid, mp);
      /* pressure contribution from material m */
      press_v(i) += press_pte(tid, mp) * vfrac_pte(tid, mp) * p_mult;
      /* temperature contribution from material m */
      temp_v(i) += temp_pte(tid, mp) * vfrac_pte(tid, mp) * t_mult;
      const Real ie_m = sie_pte(tid, mp) * frac_mass_v(i, m) * mass_sum;
      /* sie contribution from material m */
      sie_v(i) += ie_m * s_mult;
      /* assign per material specific internal energy */
      frac_ie_v(i, m) = ie_m;
      /* assign volume fraction based on pte calculation */
      frac_vol_v(i, m) = vfrac_pte(tid, mp) * vol_v(i);
      /* calculate bulk modulus for material m */
      const Real bmod_m = eos_v(pte_idxs(tid, mp))
                              .BulkModulusFromDensityTemperature(
                                  rho_pte(tid, mp), temp_pte(tid, mp), cache[mp]);
      /* add bmod contribution from material m */
      bmod_v(i) += bmod_m * vfrac_pte(tid, mp);
      /* calculate specific heat for material m */
      const Real cv_m = ev2k * eos_v(pte_idxs(tid, mp))
                                   .SpecificHeatFromDensityTemperature(
                                       rho_pte(tid, mp), temp_pte(tid, mp), cache[mp]);
      /* add mass weighted contribution specific heat for material m */
      cv_v(i) += cv_m * frac_mass_v(i, m);
      /* calculate gruneisen parameter for material m */
      const Real dpde_m = eos_v(pte_idxs(tid, mp))
                              .GruneisenParamFromDensityTemperature(
                                  rho_pte(tid, mp), temp_pte(tid, mp), cache[mp]);
      /* add gruneisen param contribution from material m */
      dpde_v(i) += dpde_m * vfrac_pte(tid, mp);
      /* optionally assign per material quantities to per material arrays */
      if (do_frac_bmod) {
        frac_bmod_v(i, m) = bmod_m;
      }
      if (do_frac_cv) {
        frac_cv_v(i, m) = cv_m;
      }
      if (do_frac_dpde) {
        frac_dpde_v(i, m) = dpde_m;
      }
    }
    if (do_t) {
      temp_v(i) /= ev2k;
    }
    if (do_s) {
      sie_v(i) /= mass_sum;
    }
    check_all_vals(i, npte, tid, pte_converged);
    /* reset mass fractions to original values if not normalized to 1 */
    for (int m = 0; m < nmat; ++m) {
      frac_mass_v(i, m) *= mass_sum;
    }
    return;
  }

 private:
  // Only run these checks when built with debug flags
  PORTABLE_FORCEINLINE_FUNCTION
  void check_all_vals(int const i, int const npte, int const tid,
                      bool const pte_converged) const {
#ifndef NDEBUG
    bool any_bad_vals = false;
    for (auto mp = 0; mp < npte; ++mp) {
      const auto m = pte_mats(tid, mp);
      any_bad_vals =
          any_bad_vals || error_utils::bad_value(frac_mass_v(i, m), "frac_mass");
      any_bad_vals = any_bad_vals || error_utils::bad_value(frac_ie_v(i, m), "frac_ie");
      any_bad_vals = any_bad_vals || error_utils::bad_value(frac_vol_v(i, m), "frac_vol");
      any_bad_vals =
          any_bad_vals || error_utils::negative_value(frac_mass_v(i, m), "frac_mass");
      any_bad_vals =
          any_bad_vals || error_utils::non_positive_value(frac_vol_v(i, m), "frac_vol");
      if (do_frac_bmod) {
        any_bad_vals =
            any_bad_vals || error_utils::bad_value(frac_bmod_v(i, m), "frac_bmod");
        any_bad_vals =
            any_bad_vals || error_utils::negative_value(frac_bmod_v(i, m), "frac_bmod");
      }
      if (do_frac_cv) {
        any_bad_vals = any_bad_vals || error_utils::bad_value(frac_cv_v(i, m), "frac_cv");
        any_bad_vals =
            any_bad_vals || error_utils::negative_value(frac_cv_v(i, m), "frac_cv");
      }
      if (do_frac_dpde) {
        any_bad_vals =
            any_bad_vals || error_utils::bad_value(frac_dpde_v(i, m), "frac_dpde");
        any_bad_vals =
            any_bad_vals || error_utils::negative_value(frac_dpde_v(i, m), "frac_dpde");
      }
    }
    any_bad_vals = any_bad_vals || error_utils::bad_value(press_v(i), "pres");
    any_bad_vals = any_bad_vals || error_utils::bad_value(temp_v(i), "temp");
    any_bad_vals = any_bad_vals || error_utils::negative_value(temp_v(i), "temp");
    any_bad_vals = any_bad_vals || error_utils::bad_value(sie_v(i), "sie");
    any_bad_vals = any_bad_vals || error_utils::bad_value(bmod_v(i), "bmod");
    any_bad_vals = any_bad_vals || error_utils::non_positive_value(bmod_v(i), "bmod");
    any_bad_vals = any_bad_vals || error_utils::bad_value(cv_v(i), "cv");
    any_bad_vals = any_bad_vals || error_utils::non_positive_value(cv_v(i), "cv");
    any_bad_vals = any_bad_vals || error_utils::bad_value(dpde_v(i), "dpde");
    any_bad_vals = any_bad_vals || error_utils::non_positive_value(dpde_v(i), "dpde");

    if (any_bad_vals) {
      using PortsOfCall::printf;
      printf("### Bad Value Output state:\n");
      printf("  ~~Bulk state~~\n");
      printf("    Pressure        : % #24.15e microbar\n", press_v(i));
      printf("    Specific Volume : % #24.15e cm^3/g\n", spvol_v(i));
      printf("    Temperature     : % #24.15g eV\n", temp_v(i));
      printf("    Specific IE     : % #24.15e erg/g\n", sie_v(i));
      printf("    Bulk Modulus    : % #24.15e microbar\n", bmod_v(i));
      printf("    Heat Capacity   : % #24.15e erg/g/eV\n", cv_v(i));
      printf("    dPdE            : % #24.15g g/cm^3\n", dpde_v(i));
      printf("  ~~Material States~~\n");
      // Maybe we can loop, but this is pretty straight-forward
      printf("    mat:");
      printf(" %24s", "Mass fraction (--)");
      printf(" %24s", "Material IE (erg)");
      printf(" %24s", "Material Volume (cm^3)");
      if (do_frac_bmod) {
        printf(" %24s", "Material BMod (ubar)");
      }
      if (do_frac_cv) {
        printf(" %24s", "Material Cv (erg/g/eV)");
      }
      if (do_frac_dpde) {
        printf(" %24s", "Material dPdE (g/cm^3)");
      }
      printf("\n");
      for (auto mp = 0; mp < npte; ++mp) {
        const auto m = pte_mats(tid, mp);
        printf("%7i:", m + 1);
        printf(" %24.15g", frac_mass_v(i, m));
        printf(" %24.15e", frac_ie_v(i, m));
        printf(" %24.15g", frac_vol_v(i, m));
        if (do_frac_bmod) {
          printf(" %24.15e", frac_bmod_v(i, m));
        }
        if (do_frac_cv) {
          printf(" %24.15e", frac_cv_v(i, m));
        }
        if (do_frac_dpde) {
          printf(" %24.15g", frac_dpde_v(i, m));
        }
        printf("\n");
      }
      if (npte > 1) {
        if (pte_converged) {
          printf("  ~~PTE iteration converged~~\n");
        } else {
          printf("  ~~PTE iteration *DID NOT* converge~~\n");
        }
      }
      PORTABLE_ALWAYS_ABORT(
          "Bad value RETURNED from singularity-eos. See output for details");
    }
  }
#endif // #ifndef NDEBUG
};

} // namespace singularity

#endif // PORTABILITY_STRATEGY_KOKKOS

#endif // _SINGULARITY_EOS_EOS_GET_SG_EOS_LAMBDAS_HPP_
