//------------------------------------------------------------------------------
// © 2022. Triad National Security, LLC. All rights reserved.  This
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

#ifndef SINGULARITY_EOS_BASE_MATH_UTILS_HPP_
#define SINGULARITY_EOS_BASE_MATH_UTILS_HPP_

#include <ports-of-call/portability.hpp>

namespace singularity {
namespace math_utils {

/*
 * Replaces "SQUARE, CUBE, etc
 * Generic is recursive, but we specialize a few of them strategically
 * A generalization would be to treat even powers specially
 */
template <size_t P>
PORTABLE_FORCEINLINE_FUNCTION constexpr auto pow(const Real x) {
  return x * pow<P - 1>(x);
}

template <>
PORTABLE_FORCEINLINE_FUNCTION constexpr auto pow<0>(const Real x) {
  return 1.0;
}

template <>
PORTABLE_FORCEINLINE_FUNCTION constexpr auto pow<1>(const Real x) {
  return x;
}

template <>
PORTABLE_FORCEINLINE_FUNCTION constexpr auto pow<2>(const Real x) {
  return x * x;
}

template <>
PORTABLE_FORCEINLINE_FUNCTION constexpr auto pow<4>(const Real x) {
  return pow<2>(x) * pow<2>(x);
}

template<class F>
double GLQuad8(F&& f, double a = -1.0, double b = 1.0) {
  static constexpr double u[4]{0.183434642495649,
                               0.525532409916329,
                               0.796666477413626,
                               0.960289856497536};
  
  static constexpr double w[4]{0.362683783378362,
                               0.313706645877887,
                               0.222381034453374,
                               0.101228536290376};
  double val = 0.0;     
  const double delta = 0.5 * (b - a);
  const double offset = 0.5 * (b + a);       
  for (int i = 0; i < 4; ++i) {
    const double xm = offset - delta * u[i];
    const double xp = offset + delta * u[i];
    val += w[i] * (f(xp) + f(-xm));  
  }
  return delta * val;
}

template<class F>
double GLQuad16(F&& f, double a = -1.0, double b = 1.0) {
  static constexpr double u[8]{ 0.0950125098376374,
                                0.2816035507792589,
                                0.4580167776572274,
                                0.6178762444026438,
                                0.7554044083550030,
                                0.8656312023878318,
                                0.9445750230732326,
                                0.9894009349916499};
  
  static constexpr double w[8]{0.1894506104550685,
                               0.1826034150449236,
                               0.1691565193950025,
                               0.1495959888165767,
                               0.1246289712555339,
                               0.0951585116824928,
                               0.0622535239386479,
                               0.0271524594117541};
  double val = 0.0;     
  const double delta = 0.5 * (b - a);
  const double offset = 0.5 * (b + a);       
  for (int i=0; i < 8; ++i) {
    const double xm = offset - delta * u[i];
    const double xp = offset + delta * u[i];
    val += w[i] * (f(xp) + f(xm));  
  }
  return delta * val;
}

template<class F>
double GLQuad32(F&& f, double a = -1.0, double b = 1.0) {
  static constexpr double u[16]{0.04830766568773834,
                                0.14447196158279654,
                                0.23928736225213719,
                                0.33186860228212779,
                                0.42135127613063534,
                                0.50689990893222944,
                                0.58771575724076238,
                                0.66304426693021528,
                                0.73218211874028970,
                                0.79448379596794240,
                                0.84936761373257008,
                                0.89632115576605218,
                                0.93490607593773976,
                                0.96476225558750646,
                                0.98561151154526849,
                                0.99726386184948169};
  
  static constexpr double w[16]{0.09654008851472785,
                                0.09563872007927495,
                                0.09384439908080466,
                                0.09117387869576396,
                                0.08765209300440385,
                                0.08331192422694675,
                                0.07819389578707037,
                                0.07234579410884857,
                                0.06582222277636189,
                                0.05868409347853559,
                                0.05099805926237628,
                                0.04283589802222678,
                                0.03427386291302149,
                                0.02539206530926219,
                                0.01627439473090571,
                                0.00701861000947011}; 
  double val = 0.0;     
  const double delta = 0.5 * (b - a);
  const double offset = 0.5 * (b + a);       
  for (int i=0; i < 16; ++i) {
    const double xm = offset - delta * u[i];
    const double xp = offset + delta * u[i];
    val += w[i] * (f(xp) + f(xm));  
  }
  return delta * val;            
}         

} // namespace math_utils
} // namespace singularity

#endif // SINGULARITY_EOS_BASE_MATH_UTILS_HPP_
