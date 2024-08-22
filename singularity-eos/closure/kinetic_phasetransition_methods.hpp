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

#ifndef _SINGULARITY_EOS_CLOSURE_KINETIC_PHSETRANSITION_METHODS_
#define _SINGULARITY_EOS_CLOSURE_KINETIC_PHSETRANSITION_METHODS_

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

PORTABLE_INLINE_FUNCTION void SmallStepMFUpdate(const Real logdt, const int num_phases,
                                                const Real *massfractions,
                                                const int *gibbsorder, const Real *logRjk,
                                                Real *dmfs, Real *newmfs) {

  Real minmassfraction = 1.e-10;

  // In logRjk: First is highest level to levels below, largest gibbs energy diff first.
  // Then follows 2nd highest to levels below. logRjk[jk], with jk = (j+1)(num_phases-1) -
  // (j-1)j/2 - k, contains rate for phase j+1 to k+1 where j=0 is highest gibbs free
  // energy state, and j=num_phases-1 is lowest gibbs free energy state. jk=0 is rate from
  // highest (j=0) to lowest (k=num_phases-1) gibbs energy states. dmfs[jk] is mass
  // transport from phase j+1 to k+1, in the same way. Fill lowest level first and stop
  // when origin level is exhausted. This will allow for jumping phases if needed. I also
  // think it is most physical.
  int jk = 0;
  Real dx = 0.;
  Real test = 0.;
  Real newtest = 0.;
  int offset = 0;
  for (int j = 0; j < num_phases - 1; j++) {
    // from high Gibb phase to low
    test = massfractions[gibbsorder[j]];
    // all phases except the higest gibbs free energy one (j=0) can get contributions from
    // phases with higher gibbs free energy.
    for (int k = 0; k < j; k++) {
      offset = (k + 1) * (num_phases - 1) - (k - 1) * k / 2;
      test = test + dmfs[offset - j];
    }
    for (int k = num_phases - 1; k > j; k--) {
      // from low Gibb phase to high, k>j always, gibbs_j > gibbs_k always
      // always mass transport from j->k
      if (test < minmassfraction) {
        dx = 0.;
        dmfs[jk++] = dx;
      } else {
        dx = std::exp(logdt + logRjk[jk] + std::log(massfractions[gibbsorder[j]]));
        newtest = test - dx;
        if (newtest < minmassfraction) {
          dmfs[jk++] = test;
        } else {
          dmfs[jk++] = dx;
        }
        test = newtest;
      }
    }
    if (test < minmassfraction) {
      newmfs[gibbsorder[j]] = 0.;
    } else {
      newmfs[gibbsorder[j]] = test;
    }
  }
  // lowest gibbs free energy phase update
  // lowest gibbs free energy phase (j=num_phases-1) can only amass what is coming from
  // larger gibbs free energy levels
  test = massfractions[gibbsorder[num_phases - 1]];
  for (int k = 0; k < num_phases - 1; k++) {
    offset = (k + 1) * (num_phases - 1) - (k - 1) * k / 2;
    test = test + dmfs[offset - num_phases + 1];
  }
  newmfs[gibbsorder[num_phases - 1]] = test;
  return;
}

} // namespace singularity

#endif // _SINGULARITY_EOS_CLOSURE_KINETIC_PHSETRANSITION_METHODS_
