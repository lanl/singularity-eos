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

#ifndef _SINGULARITY_EOS_UTILS_NUM_INTEGRATION_HPP_
#define _SINGULARITY_EOS_UTILS_NUM_INTEGRATION_HPP_

#include <cmath>
#include <cstdio>
#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <ports-of-call/math_utils.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <tuple>


#define SINGULARITY_MY_SIGN(x) (x > 0) - (x < 0)

namespace numIntegration {

//Simpsons 3/8ths rule
template <class Integrand>
PORTABLE_INLINE_FUNCTION Real integrateSimp38(Integrand&& f, Real x1, Real x2){
    Real dx = (x2 - x1)/3.0;
    Real I = f(x1) + 3.0*f(x1+dx) + 3.0*f(x2-dx) + f(x2);
    return (x2-x1)*I/8.0;
}
template <class Integrand>
PORTABLE_INLINE_FUNCTION Real integrateSimp38(Integrand&& f, double x1, double x2, int subdivides){

    Real a = x1;
    Real I = 0.0;
    Real f0 = f(a);
    Real dx = (x2 - x1);
    Real h = dx/(3.0*static_cast<double>(subdivides));

    const Real fac = 3.0*h/8.0;    

    for (int i=0;i<subdivides;i++){
        Real f1 = f(a + h);
        Real f2 = f(a + h + h);
        Real f3 = f(a + h + h + h);

        I += fac*(f0 + 3.0*(f1 + f2) + f3);
        f0 = f3;
        a += 3.0*h;
    }

    return I;
}


} // end namespace numIntegration
#endif //_SINGULARITY_EOS_UTILS_NUM_INTEGRATION_HPP_