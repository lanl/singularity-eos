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

 public:
  init_functor() = default;
  init_functor(dev_frac_v &frac_mass_v_, ScratchV<int> &pte_idxs_,
               indirection_v &eos_offsets_v_, dev_frac_v &frac_vol_v_,
               dev_frac_v frac_ie_v_, ScratchV<int> &pte_mats_,
               ScratchV<double> &vfrac_pte_, ScratchV<double> &sie_pte_,
               ScratchV<double> &temp_pte_, ScratchV<double> &press_pte_,
               ScratchV<double> &rho_pte_, dev_v &spvol_v_, dev_v &temp_v_,
               dev_v &press_v_, dev_v &sie_v_, int &nmat)
      : frac_mass_v{frac_mass_v_}, pte_idxs{pte_idxs_}, eos_offsets_v{eos_offsets_v_},
        frac_vol_v{frac_vol_v_}, frac_ie_v{frac_ie_v_}, pte_mats{pte_mats_},
        vfrac_pte{vfrac_pte_}, sie_pte{sie_pte_}, temp_pte{temp_pte_},
        press_pte{press_pte_}, rho_pte{rho_pte_}, spvol_v{spvol_v_}, temp_v{temp_v_},
        press_v{press_v_}, sie_v{sie_v_}, nmat{nmat} {}

  PORTABLE_INLINE_FUNCTION
  void operator()(const int i, const int tid, double &mass_sum, int &npte,
                  const Real t_mult, const Real s_mult, const Real p_mult) const {
    // Debug checks for inputs
    check_val(spvol_v(i));
    check_val(sie_v(i));
    check_val(press_v(i));
    // check_val(temp_v(i));
    /* normalize mass fractions */
    /* first find the mass sum */
    /* also set idxs as the decrement of the eos offsets */
    /* to take into account 1 based indexing in fortran */
    for (int m = 0; m < nmat; ++m) {
      mass_sum += frac_mass_v(i, m);
      pte_idxs(tid, m) = eos_offsets_v(m) - 1;
      frac_vol_v(i, m) = 0.0;
    }
    for (int m = 0; m < nmat; ++m) {
      frac_mass_v(i, m) /= mass_sum;
      check_val(frac_mass_v(i, m));
    }
    // count the number of participating materials and zero the inputs
    npte = 0;
    for (int m = 0; m < nmat; ++m) {
      if (frac_mass_v(i, m) > 1.e-12) {
        // participating materials are those with non-negligible mass fractions
        pte_idxs(tid, npte) = eos_offsets_v(m) - 1;
        pte_mats(tid, npte) = m;
        npte += 1;
      } else {
        frac_ie_v(i, m) = 0.0;
      }
      // zero the inputs
      vfrac_pte(tid, m) = 0.0;
      sie_pte(tid, m) = 0.0;
      temp_pte(tid, m) = 0.0;
      press_pte(tid, m) = 0.0;
    }
    // Populate the inputs with consistent values.
    // NOTE: the volume fractions and densities need to be consistent with the
    // total specific volume since they are used to calculate internal
    // quantities for the PTE solver
    for (int mp = 0; mp < npte; ++mp) {
      const int m = pte_mats(tid, mp);
      // Need to guess volume fractions
      vfrac_pte(tid, mp) = frac_mass_v(i, m);
      // Calculate densities to be consistent with these volume fractions
      rho_pte(tid, mp) = frac_mass_v(i, m) / spvol_v(i) / vfrac_pte(tid, mp);
      temp_pte(tid, mp) = temp_v(i) * ev2k * t_mult;
      press_pte(tid, mp) = press_v(i) * p_mult;
      sie_pte(tid, mp) = sie_v(i) * frac_mass_v(i, m) * s_mult;
    }
    return;
  }
 private:
  // Debug check to make sure a value is normal or zero
  template<typename valT>
  PORTABLE_FORCEINLINE_FUNCTION
  void check_val(valT value) const {
    PORTABLE_ALWAYS_REQUIRE(value == valT{0} || (std::isnormal(1e50 * value) && std::isnormal(value)),
                     "Bad value input to singularity-eos interface");
  }
};

struct final_functor {
 private:
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
  final_functor(dev_v &temp_v_, dev_v &press_v_, dev_v &sie_v_, dev_v &bmod_v_,
                dev_v &cv_v_, dev_v &dpde_v_, ScratchV<int> &pte_mats_,
                ScratchV<double> &press_pte_, ScratchV<double> &vfrac_pte_,
                ScratchV<double> &temp_pte_, ScratchV<double> &sie_pte_,
                dev_frac_v &frac_mass_v_, dev_frac_v &frac_ie_v_, dev_frac_v &frac_vol_v_,
                dev_v &vol_v_, Kokkos::View<EOS *, Llft> &eos_v_,
                ScratchV<int> &pte_idxs_, ScratchV<double> &rho_pte_,
                dev_frac_v &frac_bmod_v_, dev_frac_v &frac_cv_v_,
                dev_frac_v &frac_dpde_v_, int &nmat_, bool do_frac_bmod_,
                bool do_frac_cv_, bool do_frac_dpde_)
      : temp_v{temp_v_}, press_v{press_v_}, sie_v{sie_v_}, bmod_v{bmod_v_}, cv_v{cv_v_},
        dpde_v{dpde_v_}, pte_mats{pte_mats_}, press_pte{press_pte_},
        vfrac_pte{vfrac_pte_}, temp_pte{temp_pte_}, sie_pte{sie_pte_},
        frac_mass_v{frac_mass_v_}, frac_ie_v{frac_ie_v_}, frac_vol_v{frac_vol_v_},
        vol_v{vol_v_}, eos_v{eos_v_}, pte_idxs{pte_idxs_}, rho_pte{rho_pte_},
        frac_bmod_v{frac_bmod_v_}, frac_cv_v{frac_cv_v_},
        frac_dpde_v{frac_dpde_v_}, nmat{nmat_}, do_frac_bmod{do_frac_bmod_},
        do_frac_cv{do_frac_cv_}, do_frac_dpde{do_frac_dpde_} {}

 public:
  PORTABLE_INLINE_FUNCTION
  void operator()(const int i, const int tid, const int npte, const Real mass_sum,
                  const Real t_mult, const Real s_mult, const Real p_mult,
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
      check_val(frac_ie_v(i, m));
      /* assign volume fraction based on pte calculation */
      frac_vol_v(i, m) = vfrac_pte(tid, mp) * vol_v(i);
      check_val(frac_vol_v(i, m));
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
        check_val(frac_bmod_v(i, m));
      }
      if (do_frac_cv) {
        frac_cv_v(i, m) = cv_m;
        check_val(frac_cv_v(i, m));
      }
      if (do_frac_dpde) {
        frac_dpde_v(i, m) = dpde_m;
        check_val(frac_cv_v(i, m));
      }
      check_val(frac_mass_v(i, m));
    }
    if (do_t) {
      temp_v(i) /= ev2k;
    }
    if (do_s) {
      sie_v(i) /= mass_sum;
    }
    /* reset mass fractions to original values if not normalized to 1 */
    for (int m = 0; m < nmat; ++m) {
      frac_mass_v(i, m) *= mass_sum;
    }
    check_val(press_v(i));
    check_val(temp_v(i));
    check_val(sie_v(i));
    check_val(bmod_v(i));
    check_val(cv_v(i));
    check_val(dpde_v(i));
    return;
  }
 private:
  // Debug check to make sure returned values are normal or zero. Extra NDEBUG
  // check probably not required, but it's there just in case the compiler tries
  // to evaluate the boolean
  template<typename valT>
  PORTABLE_FORCEINLINE_FUNCTION
  void check_val(valT value) const {
    PORTABLE_ALWAYS_REQUIRE(value == valT{0} || (std::isnormal(1e50 * value) && std::isnormal(value)),
                     "Bad value returned from singularity-eos interface");
  }
};

} // namespace singularity

#endif // PORTABILITY_STRATEGY_KOKKOS

#endif // _SINGULARITY_EOS_EOS_GET_SG_EOS_LAMBDAS_HPP_
