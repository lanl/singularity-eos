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

namespace singularity {
void get_sg_eos_rho_t(const char *name, int ncell, indirection_v &offsets_v,
                      Kokkos::View<EOS *, Llft> &eos_v, dev_v &press_v, dev_v &pmax_v,
                      dev_v &sie_v, dev_frac_v &frac_mass_v, ScratchV<int> &pte_idxs,
                      ScratchV<int> &pte_mats, ScratchV<double> &press_pte,
                      ScratchV<double> &vfrac_pte, ScratchV<double> &rho_pte,
                      ScratchV<double> &sie_pte, ScratchV<double> &temp_pte,
                      ScratchV<double> &solver_scratch,
                      Kokkos::Experimental::UniqueToken<DES, KGlobal> &tokens,
                      bool small_loop, init_functor &i_func, final_functor &f_func) {
  portableFor(
      name, 0, ncell, PORTABLE_LAMBDA(const int &iloop) {
        // cell offset
        const int i{offsets_v(iloop) - 1};
        // get "thread-id" like thing with optimization
        // for small loops
        const int32_t token{tokens.acquire()};
        const int32_t tid{small_loop ? iloop : token};
        double mass_sum{0.0};
        int npte{0};
        // initialize values for solver / lookup
        i_func(i, tid, mass_sum, npte, 1.0, 0.0, 0.0);
        // calculate pte condition (lookup for 1 mat cell)
        Real sie_tot_true{0.0};
        const int neq = npte;
        singularity::mix_impl::CacheAccessor cache(&solver_scratch(tid, 0) +
                                                   neq * (neq + 4) + 2 * npte);
        if (npte > 1) {
          // create solver lambda
          // eos accessor
          singularity::EOSAccessor_ eos_inx(eos_v, &pte_idxs(tid, 0));
          PTESolverFixedT<singularity::EOSAccessor_, Real *, Real **> method(
              npte, eos_inx, 1.0, temp_pte(tid, 0), &rho_pte(tid, 0), &vfrac_pte(tid, 0),
              &sie_pte(tid, 0), &temp_pte(tid, 0), &press_pte(tid, 0), cache,
              &solver_scratch(tid, 0));
          const bool res_{PTESolver(method)};
          // calculate total internal energy
          for (int mp = 0; mp < npte; ++mp) {
            const int m = pte_mats(tid, mp);
            sie_tot_true += sie_pte(tid, mp) * frac_mass_v(i, m);
          }
        } else {
          // pure cell (nmat = 1)
          // calculate sie from single eos
          sie_pte(tid, 0) = eos_v(pte_idxs(tid, 0))
                                .InternalEnergyFromDensityTemperature(
                                    rho_pte(tid, 0), temp_pte(tid, 0), cache[0]);
          // set total sie to material 0 value
          sie_tot_true = sie_pte(tid, 0);
          // set pressure
          press_pte(tid, 0) = eos_v(pte_idxs(tid, 0))
                                  .PressureFromDensityTemperature(
                                      rho_pte(tid, 0), temp_pte(tid, 0), cache[0]);
        }
        // total sie is known
        sie_v(i) = sie_tot_true;
        // assign quantities
        f_func(i, tid, npte, mass_sum, 0.0, 0.0, 1.0, cache);
        // assign max pressure
        pmax_v(i) = press_v(i) > pmax_v(i) ? press_v(i) : pmax_v(i);
        // release the token used for scratch arrays
        tokens.release(token);
      });
  return;
}
} // namespace singularity
