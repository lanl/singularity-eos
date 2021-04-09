//------------------------------------------------------------------------------
// Extra definitions to the SP5 file format required for singularity-eos
// Author: Jonah Miller (jonahm@lanl.gov)
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


#ifndef _SINGULARITY_EOS_UTILS_SP5_SINGULARITY_EOS_SP5_HPP_
#define _SINGULARITY_EOS_UTILS_SP5_SINGULARITY_EOS_SP5_HPP_

namespace SP5 {

  constexpr char defaultSesFileName[]  = "sesame_table.sp5";
  
  namespace Depends {
    constexpr char logRhoLogSie[] = "dependsLogRhoLogSie";
    constexpr char logRhoLogT[]   = "dependsLogRhoLogT";
    constexpr char coldCurve[]    = "coldCurve";
  }
  
  namespace Offsets {
    constexpr char messageName[] = "interpretation";
    constexpr char message[]
    = "All quantities are functions of log_10(X)\n"
      "for X = density rho, temperature T, or specific internal energy sie\n"
      "where conversion is X = 10^{Xlog} - Xoffset\n";
    constexpr char rho[]    = "rhoOffset";
    constexpr char T[]      = "TOffset";
    constexpr char sie[]    = "sieOffset";
  }
  
  namespace Material {
    constexpr char exchangeCoefficient[] = "exchangeCoefficient";
    constexpr char meanAtomicMass[]      = "meanAtomicMass";
    constexpr char meanAtomicNumber[]    = "meanAtomicNumber";
    constexpr char solidBulkModulus[]    = "solidBulkModulus";
    constexpr char normalDensity[]       = "normalDensity";
    constexpr char comments[]            = "comments";
    constexpr char matid[]               = "matid";
    constexpr char name[]                = "name";
  }
  
  namespace Fields {
    constexpr char P[]              = "pressure";
    constexpr char sie[]            = "specific internal energy";
    constexpr char T[]              = "temperature";
    constexpr char bMod[]           = "bulk modulus";
    constexpr char dPdRho[]         = "dPdRho";
    constexpr char dPdE[]           = "dPdE";
    constexpr char dTdRho[]         = "dTdRho";
    constexpr char dTdE[]           = "dTdE";
    constexpr char dEdRho[]         = "dEdRho";
    constexpr char dEdT[]           = "dEdT";
    constexpr char mask[]           = "mask";
    constexpr char transitionMask[] = "transition mask";
  }
  
}

#endif // _SINGULARITY_EOS_UTILS_SP5_SINGULARITY_EOS_SP5_HPP_
