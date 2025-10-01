//------------------------------------------------------------------------------
// © 2024. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_CLOSURE_KINETIC_PHSETRANSITION_UTILS_
#define _SINGULARITY_EOS_CLOSURE_KINETIC_PHSETRANSITION_UTILS_

#ifdef SINGULARITY_BUILD_CLOSURE

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/fast-math/logs.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/eos/eos.hpp>

#include <cmath>

#ifdef SINGULARITY_USE_KOKKOSKERNELS
#include <KokkosBatched_ApplyQ_Decl.hpp>
#include <KokkosBatched_QR_Decl.hpp>
#include <KokkosBatched_Trsv_Decl.hpp>
#else
#include <Eigen/Dense>
#endif // SINGULARITY_USE_KOKKOSKERNELS

namespace singularity {

constexpr Real KPT_MIN_MASS_FRACTION = 1.0e-10;

template<typename ConstRealIndexer, typename IntIndexer>
PORTABLE_INLINE_FUNCTION void SortGibbs(const int num_phases, ConstRealIndexer &&gibbs,
                                        IntIndexer &&gibbsorder) {
  int itmp = 0;
  // initiate gibbsorder
  for (int j = 0; j < num_phases; j++) {
    gibbsorder[j] = j;
  }
  // sort from high Gibb phase to low
  for (int j = 0; j < num_phases - 1; j++) {
    for (int k = num_phases - 1; k > j; k--) {
      // from low Gibb phase to high, k>j always, gibbs_j > gibbs_k always
      if (gibbs[gibbsorder[k]] > gibbs[gibbsorder[j]]) {
        itmp = gibbsorder[j];
        gibbsorder[j] = gibbsorder[k];
        gibbsorder[k] = itmp;
      }
    }
  }
  return;
}

} // namespace singularity

#endif // SINGULARITY_BUILD_CLOSURE
#endif // _SINGULARITY_EOS_CLOSURE_KINETIC_PHSETRANSITION_UTILS_
