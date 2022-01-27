//------------------------------------------------------------------------------
// © 2021. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_EOS_VECTOR_UTILS_
#define _SINGULARITY_EOS_EOS_VECTOR_UTILS_

#include <ports-of-call/portability.hpp>

namespace singularity {
namespace vector_utils {

/*
Class where the [] operator returns the null pointer for any index
*/
class NullIndexer {
 public:
  PORTABLE_FORCEINLINE_FUNCTION
  Real *operator[](int i) {
    return nullptr;
  }
};

/*
Wrappers for scalar lookups to allow them to take vectors on the CPU but then
call GPU kernels to do the actual lookups

RealIndexer must have an operator[](int) that returns a Real. e.g., Real*
ConstRealIndexer is as RealIndexer, but assumed const type.
LambdaIndexer must have an operator[](int) that returns a Real*. e.g., Real**
*/
template<typename EOS, typename LookupFunction, typename ConstRealIndexer,
         typename RealIndexer, typename LambdaIndexer>
inline
void VectorizeScalarLookup(const EOS &eos, const LookupFunction &&f,
                           ConstRealIndexer &&xs, ConstRealIndexer &&ys,
                           RealIndexer &&result, const int num,
                           LambdaIndexer &&lambdas) {
  portableFor(
      'VectorScalarLookup', 0, num, PORTABLE_LAMBDA(const int i) {
        result[i] = (eos.*f)(xs[i], ys[i], lambdas[i]);
      }
  );
}

template<typename EOS, typename FillEosFunction, typename RealIndexer,
         typename LambdaIndexer>
inline
void VectorizeFillEos(const EOS &eos, const FillEosFunction &&f,
                      RealIndexer &&rhos, RealIndexer &&temps,
                      RealIndexer &&energies, RealIndexer &&presses,
                      RealIndexer &&cvs, RealIndexer &&bmods, const int num,
                      const unsigned long output, LambdaIndexer &&lambdas) {
  portableFor(
      'VectorFillEos', 0, num, PORTABLE_LAMBDA(const int i) {
        (eos.*f)(rhos[i], temps[i], energies[i], presses[i], cvs[i], bmods[i],
                 output, lambdas[i]);
      }
  );
}
} // vector_utils
} // singularity

#endif // _SINGULARITY_EOS_EOS_VECTOR_UTILS_
