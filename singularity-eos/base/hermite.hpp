//------------------------------------------------------------------------------
// © 2023. Triad National Security, LLC. All rights reserved.  This
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

#ifndef SINGULARITY_EOS_BASE_HERMITE_HPP_
#define SINGULARITY_EOS_BASE_HERMITE_HPP_

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/math_utils.hpp>

namespace singularity {
namespace hermite {

// psi0 and its derivatives
PORTABLE_INLINE_FUNCTION Real psi0(Real z) {
  return z * z * z * (z * (-6.0 * z + 15.0) - 10.0) + 1.0;
}
PORTABLE_INLINE_FUNCTION Real dpsi0(Real z) {
  return z * z * (z * (-30.0 * z + 60.0) - 30.0);
}
PORTABLE_INLINE_FUNCTION Real ddpsi0(Real z) {
  return z * (z * (-120.0 * z + 180.0) - 60.0);
}

// psi1 and its derivatives
PORTABLE_INLINE_FUNCTION Real psi1(Real z) {
  return z * (z * z * (z * (-3.0 * z + 8.0) - 6.0) + 1.0);
}
PORTABLE_INLINE_FUNCTION Real dpsi1(Real z) {
  return z * z * (z * (-15.0 * z + 32.0) - 18.0) + 1.0;
}
PORTABLE_INLINE_FUNCTION Real ddpsi1(Real z) {
  return z * (z * (-60.0 * z + 96.0) - 36.0);
}

// psi2 and its derivatives
PORTABLE_INLINE_FUNCTION Real psi2(Real z) {
  return 0.5 * z * z * (z * (z * (-z + 3.0) - 3.0) + 1.0);
}
PORTABLE_INLINE_FUNCTION Real dpsi2(Real z) {
  return 0.5 * z * (z * (z * (-5.0 * z + 12.0) - 9.0) + 2.0);
}
PORTABLE_INLINE_FUNCTION Real ddpsi2(Real z) {
  return 0.5 * (z * (z * (-20.0 * z + 36.0) - 18.0) + 2.0);
}

// biquintic hermite polynomial
PORTABLE_INLINE_FUNCTION Real h5(Real fi[36], Real w0t, Real w1t, Real w2t, Real w0mt,
                                 Real w1mt, Real w2mt, Real w0d, Real w1d, Real w2d,
                                 Real w0md, Real w1md, Real w2md) {
  return fi[0] * w0d * w0t + fi[1] * w0md * w0t + fi[2] * w0d * w0mt +
         fi[3] * w0md * w0mt + fi[4] * w0d * w1t + fi[5] * w0md * w1t +
         fi[6] * w0d * w1mt + fi[7] * w0md * w1mt + fi[8] * w0d * w2t +
         fi[9] * w0md * w2t + fi[10] * w0d * w2mt + fi[11] * w0md * w2mt +
         fi[12] * w1d * w0t + fi[13] * w1md * w0t + fi[14] * w1d * w0mt +
         fi[15] * w1md * w0mt + fi[16] * w2d * w0t + fi[17] * w2md * w0t +
         fi[18] * w2d * w0mt + fi[19] * w2md * w0mt + fi[20] * w1d * w1t +
         fi[21] * w1md * w1t + fi[22] * w1d * w1mt + fi[23] * w1md * w1mt +
         fi[24] * w2d * w1t + fi[25] * w2md * w1t + fi[26] * w2d * w1mt +
         fi[27] * w2md * w1mt + fi[28] * w1d * w2t + fi[29] * w1md * w2t +
         fi[30] * w1d * w2mt + fi[31] * w1md * w2mt + fi[32] * w2d * w2t +
         fi[33] * w2md * w2t + fi[34] * w2d * w2mt + fi[35] * w2md * w2mt;
}

// cubic hermite polynomial
// psi0 and its derivatives
PORTABLE_INLINE_FUNCTION Real xpsi0(Real z) { return z * z * (2.0 * z - 3.0) + 1.0; }
PORTABLE_INLINE_FUNCTION Real xdpsi0(Real z) { return z * (6.0 * z - 6.0); }

// psi1 and its derivatives
PORTABLE_INLINE_FUNCTION Real xpsi1(Real z) { return z * (z * (z - 2.0) + 1.0); }

PORTABLE_INLINE_FUNCTION Real xdpsi1(Real z) { return z * (3.0 * z - 4.0) + 1.0; }

// bicubic hermite polynomial
PORTABLE_INLINE_FUNCTION Real h3(Real fi[16], Real w0t, Real w1t, Real w0mt, Real w1mt,
                                 Real w0d, Real w1d, Real w0md, Real w1md) {
  return fi[0] * w0d * w0t + fi[1] * w0md * w0t + fi[2] * w0d * w0mt +
         fi[3] * w0md * w0mt + fi[4] * w0d * w1t + fi[5] * w0md * w1t +
         fi[6] * w0d * w1mt + fi[7] * w0md * w1mt + fi[8] * w1d * w0t +
         fi[9] * w1md * w0t + fi[10] * w1d * w0mt + fi[11] * w1md * w0mt +
         fi[12] * w1d * w1t + fi[13] * w1md * w1t + fi[14] * w1d * w1mt +
         fi[15] * w1md * w1mt;
}

} // namespace hermite
} // namespace singularity

#endif // SINGULARITY_EOS_BASE_MATH_UTILS_HPP_
