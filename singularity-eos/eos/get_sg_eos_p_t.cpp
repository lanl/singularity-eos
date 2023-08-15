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

#include <ports-of-call/portability.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/get_sg_eos.hpp>
#include <singularity-eos/eos/get_sg_eos_functors.hpp>

// any bmod vars, dpde, cv, frac vol, frac ie

namespace singularity {
void get_sg_eos_p_t(const char *name, int ncell, int nmat, indirection_v &offsets_v,
                    indirection_v &eos_offsets_v, Kokkos::View<EOS *, Llft> &eos_v,
                    dev_v &press_v, dev_v &pmax_v, dev_v &vol_v, dev_v &spvol_v,
                    dev_v &sie_v, dev_v &temp_v, dev_frac_v &frac_mass_v,
                    ScratchV<int> &pte_idxs, ScratchV<int> &pte_mats,
                    ScratchV<double> &press_pte, ScratchV<double> &vfrac_pte,
                    ScratchV<double> &rho_pte, ScratchV<double> &sie_pte,
                    ScratchV<double> &temp_pte, ScratchV<double> &solver_scratch,
                    Kokkos::Experimental::UniqueToken<DES, KGlobal> &tokens,
                    bool small_loop, final_functor &f_func) {
  portableFor(
      name, 0, ncell, PORTABLE_LAMBDA(const int &iloop) {
        // cell offset
        const int i{offsets_v(iloop) - 1};
        // get "thread-id" like thing with optimization
        // for small loops
        const int32_t token{tokens.acquire()};
        const int32_t tid{small_loop ? iloop : token};
        // caching mechanism
        singularity::mix_impl::CacheAccessor cache(&solver_scratch(tid, 0));
        double mass_sum{0.0};
        // normalize mass fractions
        // first find the mass sum
        // also set idxs as the decrement of the eos offsets
        // to take into account 1 based indexing in fortran
        for (int m = 0; m < nmat; ++m) {
          mass_sum += frac_mass_v(i, m);
          pte_idxs(tid, m) = eos_offsets_v(m) - 1;
          pte_mats(tid, m) = m;
          temp_pte(tid, m) = temp_v(i) * ev2k;
          press_pte(tid, m) = press_v(i);
        }
        for (int m = 0; m < nmat; ++m) {
          frac_mass_v(i, m) /= mass_sum;
        }
        // do r-e of pt for each mat
        singularity::EOSAccessor_ eos_inx(eos_v, &pte_idxs(tid, 0));
        Real vfrac_tot{0.0};
        Real sie_tot{0.0};
        for (int m = 0; m < nmat; ++m) {
          // obtain rho and sie from P-T
          eos_inx[m].DensityEnergyFromPressureTemperature(
              press_pte(tid, m), temp_pte(tid, m), cache[m], rho_pte(tid, m),
              sie_pte(tid, m));
          // assign volume fractions
          // this is a physical volume
          vfrac_pte(tid, m) = frac_mass_v(i, m) / rho_pte(tid, m) * mass_sum;
          vfrac_tot += vfrac_pte(tid, m);
          // add internal energy component
          sie_tot += sie_pte(tid, m) * frac_mass_v(i, m);
        }
        // assign volume, etc.
        // total sie is known
        sie_v(i) = sie_tot;
        vol_v(i) = vfrac_tot;
        spvol_v(i) = vol_v(i) / mass_sum;
        for (int m = 0; m < nmat; ++m) {
          vfrac_pte(tid, m) /= vfrac_tot;
        }
        // assign remaining outputs
        f_func(i, tid, nmat, mass_sum, 0.0, 0.0, 0.0, cache);
        // assign max pressure
        pmax_v(i) = press_v(i) > pmax_v(i) ? press_v(i) : pmax_v(i);
        // release the token used for scratch arrays
        tokens.release(token);
      });
  return;
}
} // namespace singularity
