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

#include <ports-of-call/portability.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/get_sg_eos.hpp>
#include <singularity-eos/eos/get_sg_eos_functors.hpp>

namespace singularity {
void get_sg_eos_rho_e(const char *name, int ncell, indirection_v &offsets_v,
                      Kokkos::View<EOS *, Llft> &eos_v, dev_v &press_v, dev_v &pmax_v,
                      dev_v &sie_v, ScratchV<int> &pte_idxs, ScratchV<double> &press_pte,
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
        double vfrac_sum{0.0};
        // initialize values for solver / lookup
        i_func(i, tid, mass_sum, npte, vfrac_sum, 0.0, 1.0, 0.0);
        // need to initialize the scratch before it's used to avoid undefined behavior
        for (std::size_t idx = 0; idx < solver_scratch.extent(1); ++idx) {
          solver_scratch(tid, idx) = 0.0;
        }
        // get cache from offsets into scratch
        // JMM: We allocate more space than needed and re-use solver
        // scratch for the cache accessor.
        const int neq = npte + 1;
        singularity::mix_impl::CacheAccessor cache(&solver_scratch(tid, 0) +
                                                   neq * (neq + 4) + 2 * npte);
        bool pte_converged = true;
        if (npte > 1) {
          // create solver lambda
          // eos accessor
          singularity::EOSAccessor_ eos_inx(eos_v, &pte_idxs(tid, 0));
          // reset inputs
          // JMM: The solver constructor is (deep under the hood)
          // capturing by reference. So to avoid out-of-scope access,
          // these must be "anchored" at caller scope.
          Real *prho_pte = &rho_pte(tid, 0);
          Real *pvfrac_pte = &vfrac_pte(tid, 0);
          Real *psie_pte = &sie_pte(tid, 0);
          Real *ptemp_pte = &temp_pte(tid, 0);
          Real *ppress_pte = &press_pte(tid, 0);
          Real *pscratch = &solver_scratch(tid, 0);
          PTESolverRhoT<singularity::EOSAccessor_, Real *, decltype(cache)> method(
              npte, eos_inx, vfrac_sum, sie_v(i), prho_pte, pvfrac_pte, psie_pte,
              ptemp_pte, ppress_pte, cache, pscratch);
          auto status = PTESolver(method);
          pte_converged = status.converged;
        } else {
          // pure cell (nmat = 1)
          temp_pte(tid, 0) = eos_v(pte_idxs(tid, 0))
                                 .TemperatureFromDensityInternalEnergy(
                                     rho_pte(tid, 0), sie_pte(tid, 0), cache[0]);
          press_pte(tid, 0) = eos_v(pte_idxs(tid, 0))
                                  .PressureFromDensityTemperature(
                                      rho_pte(tid, 0), temp_pte(tid, 0), cache[0]);
        }
        // assign outputs
        f_func(i, tid, npte, mass_sum, 1.0, 0.0, 1.0, pte_converged, cache);
        // assign max pressure
        pmax_v(i) = press_v(i) > pmax_v(i) ? press_v(i) : pmax_v(i);
        // release the token used for scratch arrays
        tokens.release(token);
      });
  return;
}
} // namespace singularity
