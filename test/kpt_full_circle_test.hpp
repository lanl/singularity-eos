//------------------------------------------------------------------------------
// © 2025. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_TEST_KPTFULLCIRCLE_TEST_
#define _SINGULARITY_EOS_TEST_KPTFULLCIRCLE_TEST_

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

// This data is taken from the .ult file for a test case in Flag,
// MPMat_KPT_3phases.ult.std in MPMat_KPT_3phases_NW. Note that 3 phases appear by mistake
// (no windows) in sn2162 where only 2 phases should be present (with windows). Note that
// at this point in time (April 2024) the global and phase times are not alinged in .ult
// but they will be in the near future.

namespace kpt_full_circle_test {

#include <test/eos_unit_test_helpers.hpp>
#include <test/pte_test_3phaseSesameSn.hpp>

using namespace pte_test_3phaseSesameSn;

constexpr Real out_gibbs0[NTRIAL] = {1.26713070e10, 1.26869849e10, 1.27015326e10,
                                     1.27135203e10, 1.27219984e10};
constexpr Real out_gibbs1[NTRIAL] = {1.25716504e10, 1.25869850e10, 1.26012338e10,
                                     1.26130226e10, 1.26214326e10};
constexpr Real out_gibbs2[NTRIAL] = {1.26707691e10, 1.26861007e10, 1.27003466e10,
                                     1.27121330e10, 1.27205413e10};

template <typename RealIndexer>
inline void gibbsprintresults(int n, RealIndexer &&rho, RealIndexer &&temp,
                              RealIndexer &&sie, RealIndexer &&gibbsrt,
                              RealIndexer &&gibbsre) {
  std::cout << "Trial number: " << n << std::endl;
  std::cout << "Density: beta, gamma, hcp: \t\t\t\t" << rho[0] << ", " << rho[1] << ", "
            << rho[2] << std::endl;
  std::cout << "Temperature: beta, gamma, hcp: \t\t\t\t" << temp[0] << ", " << temp[1] << ", "
            << temp[2] << std::endl;
  std::cout << "Gibbs(Density,Temperature): beta, gamma, hcp: \t\t" << gibbsrt[0] << ", "
            << gibbsrt[1] << ", " << gibbsrt[2] << std::endl;
  std::cout << "Internal energy: beta, gamma, hcp: \t\t\t" << sie[0] << ", " << sie[1] << ", "
            << sie[2] << std::endl;
  std::cout << "Gibbs(Density,Internal Energy): beta, gamma, hcp: \t" << gibbsre[0]
            << ", " << gibbsre[1] << ", " << gibbsre[2] << std::endl;
  std::cout << "Within 0.1 % from each other (0=False, 1=True): beta, gamma, hcp: \t "
            << isClose(gibbsrt[0], gibbsre[0], 1e-3) << ", "
            << isClose(gibbsrt[1], gibbsre[1], 1e-3) << ", "
            << isClose(gibbsrt[2], gibbsre[2], 1e-3) << std::endl;
  std::cout << "Within 0.1 % from Gibbs from LAP: beta, gamma, hcp: \t\t\t "
            << isClose(gibbsrt[0], out_gibbs0[n], 1e-3) << ", "
            << isClose(gibbsrt[1], out_gibbs1[n], 1e-3) << ", "
            << isClose(gibbsrt[2], out_gibbs2[n], 1e-3) << std::endl
            << std::endl;
}

template <typename RealIndexer>
inline void kptprintinput(int n, RealIndexer &&rho, RealIndexer &&vfrac) {
  std::cout << "Trial number: " << n << std::endl;
  std::cout << "Total Specific Internal energy: \t\t\t" << in_sie_tot[n] << std::endl;
  std::cout << "Total density: \t\t\t\t\t" << in_rho_tot[n] << std::endl;
  std::cout << "Mass fractions: beta, gamma, hcp: \t\t\t" << in_lambda[0][n] << ", "
            << in_lambda[1][n] << ", " << in_lambda[2][n] << std::endl;
  std::cout << "Assuming volume fractions: beta, gamma, hcp: \t" << vfrac[0] << ", "
            << vfrac[1] << ", " << vfrac[2] << std::endl;
  std::cout << "gives starting phase densities: beta, gamma, hcp: \t" << rho[0] << ", "
            << rho[1] << ", " << rho[2] << std::endl
            << std::endl;
}

template <typename RealIndexer>
inline void kptprintresults(int n, RealIndexer &&rho, RealIndexer &&vfrac,
                            RealIndexer &&sie, RealIndexer &&press, RealIndexer &&temp) {
  std::cout << "Trial number: " << n << std::endl;
  std::cout << "Total Specific Internal energy: \t"
            << sie[0] * in_lambda[0][n] + sie[1] * in_lambda[1][n] +
                   sie[2] * in_lambda[2][n]
            << ", (" << in_sie_tot[n] << ")" << std::endl;
  std::cout << "Total density: \t\t\t"
            << 1.0 / (1.0 / rho[0] * in_lambda[0][n] + 1.0 / rho[1] * in_lambda[1][n] +
                      1.0 / rho[2] * in_lambda[2][n])
            << ", (" << in_rho_tot[n] << ")" << std::endl;
  std::cout << "Volume fractions: beta, gamma, hcp: " << vfrac[0] << ", " << vfrac[1]
            << ", " << vfrac[2] << std::endl;
  std::cout << "Density: beta, gamma, hcp: \t\t" << rho[0] << ", " << rho[1] << ", "
            << rho[2] << ", (" << out_rho0[n] << ", " << out_rho1[n] << ", "
            << out_rho2[n] << ")" << std::endl;
  std::cout << "Pressure: beta, gamma, hcp: \t" << press[0] << ", " << press[1] << ", "
            << press[2] << ", (" << out_press[n] << ")" << std::endl;
  std::cout << "Temperature: beta, gamma, hcp: \t" << temp[0] << ", " << temp[1] << ", "
            << temp[2] << ", (" << out_temp[n] << ")" << std::endl;
  std::cout << "Internal energy: beta, gamma, hcp: \t" << sie[0] << ", " << sie[1] << ", "
            << sie[2] << ", (" << out_sie0[n] << ", " << out_sie1[n] << ", "
            << out_sie2[n] << ")" << std::endl
            << std::endl;
}

} // end namespace kpt_full_circle_test

#endif // _SINGULARITY_EOS_TEST_KPTFULLCIRCLE_TEST_
