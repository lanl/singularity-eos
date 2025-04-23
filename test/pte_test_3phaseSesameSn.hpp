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

#ifndef _SINGULARITY_EOS_TEST_PTE_TEST_3PHASESESAMESN_
#define _SINGULARITY_EOS_TEST_PTE_TEST_3PHASESESAMESN_

#include <stdlib.h>

// For getting the symlink in
#include <unistd.h>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>

#include <singularity-eos/eos/eos_models.hpp>
#include <singularity-eos/eos/eos_variant.hpp>

// This data is taken from the .ult file for a test case in Flag,
// MPMat_KPT_3phases.ult.std in MPMat_KPT_3phases_NW. Note that 3 phases appear by mistake
// (no windows) in sn2162 where only 2 phases should be present (with windows). Note that
// at this point in time (April 2024) the global and phase times are not alinged in .ult
// but they will be in the near future.

#include <material_5phaseSesameSn.hpp>

namespace pte_test_3phase {

constexpr int NMAT = 3;
std::vector<int> PHASES = {0, 1, 3};

constexpr int NTRIAL = 5;
constexpr int NPTS = NTRIAL * NMAT;
constexpr int HIST_SIZE = 30;

constexpr Real in_rho_tot[NTRIAL] = {8.41382588, 8.41550873, 8.41719191, 8.41887543,
                                     8.42055928};
constexpr Real in_sie_tot[NTRIAL] = {1.67602914e09, 1.67876139e09, 1.68149578e09,
                                     1.68423198e09, 1.68696932e09};
constexpr Real in_lambda[NMAT][NTRIAL] = {
    {0.663852654, 0.660579238, 0.656880511, 0.652218675, 0.646232106},
    {0.336093276, 0.339262144, 0.342703726, 0.346800126, 0.352135167},
    {0.000054070, 0.000158618, 0.000415763, 0.000982199, 0.001632727}};

constexpr Real trial_vfrac[NMAT] = {650.0 / 1000.0, 349.5 / 1000.0, 0.5 / 1000.0};

constexpr Real out_press[NTRIAL] = {1.0946046e11, 1.0959183e11, 1.0971241e11,
                                    1.0980857e11, 1.0987166e11};
constexpr Real out_temp[NTRIAL] = {438.96969, 438.98046, 438.95507, 438.84886, 438.64322};
constexpr Real out_rho0[NTRIAL] = {8.35319483, 8.35425038, 8.35523024, 8.35603912,
                                   8.35661359};
constexpr Real out_rho1[NTRIAL] = {8.53618806, 8.53734120, 8.53841145, 8.53929455,
                                   8.53992118};
constexpr Real out_rho2[NTRIAL] = {8.53854993, 8.53970495, 8.54077710, 8.54166211,
                                   8.54229063};
constexpr Real out_sie0[NTRIAL] = {1.54379257e09, 1.54520613e09, 1.54644403e09,
                                   1.54728448e09, 1.54760327e09};
constexpr Real out_sie1[NTRIAL] = {1.93716882e09, 1.93864981e09, 1.93994981e09,
                                   1.94084038e09, 1.94119298e09};
constexpr Real out_sie2[NTRIAL] = {2.01473728e09, 2.01621655e09, 2.01751522e09,
                                   2.01840528e09, 2.01875845e09};

template <typename RealIndexer>
inline void set_trial_state(int n, RealIndexer &&rho, RealIndexer &&vfrac,
                             RealIndexer &&sie) {

  Real vsum = 0.;
  for (int i = 0; i < NMAT; i++) {
    vfrac[i] = trial_vfrac[i];
    rho[i] = in_lambda[i][n] * in_rho_tot[n] / vfrac[i];
    // same sie in both phases gives sie_tot=sie
    sie[i] = in_sie_tot[n];
    vsum += vfrac[i];
  }

  for (int i = 0; i < NMAT; i++)
    vfrac[i] *= 1.0 / vsum;

  return;
}

template<typename RealIndexer>
inline void printinput(int n, RealIndexer &&rho, RealIndexer &&vfrac){
      std::cout << "Trial number: " << n << std::endl;
      std::cout << "Total Specific Internal energy: \t\t\t" << in_sie_tot[n] << std::endl;
      std::cout << "Total density: \t\t\t\t\t" << in_rho_tot[n] << std::endl;
      std::cout << "Mass fractions: beta, gamma, hcp: \t\t\t" << in_lambda[0][n] << ", "
                << in_lambda[1][n] << ", " << in_lambda[2][n] << std::endl;
      std::cout << "Assuming volume fractions: beta, gamma, hcp: \t" << vfrac[0]
                << ", " << vfrac[1] << ", " << vfrac[2] << std::endl;
      std::cout << "gives starting phase densities: beta, gamma, hcp: \t" << rho[0]
                << ", " << rho[1] << ", " << rho[2] << std::endl
                << std::endl;
}

template<typename RealIndexer>
inline void printresults(int n, RealIndexer &&rho, RealIndexer &&vfrac, RealIndexer &&sie, RealIndexer &&press, RealIndexer &&temp){
      std::cout << "Trial number: " << n << std::endl;
      std::cout << "Total Specific Internal energy: \t"
                << sie[0] * in_lambda[0][n] + sie[1] * in_lambda[1][n] +
                       sie[2] * in_lambda[2][n]
                << ", (" << in_sie_tot[n] << ")" << std::endl;
      std::cout << "Total density: \t\t\t"
                << 1.0 / (1.0 / rho[0] * in_lambda[0][n] +
                          1.0 / rho[1] * in_lambda[1][n] +
                          1.0 / rho[2] * in_lambda[2][n])
                << ", (" << in_rho_tot[n] << ")" << std::endl;
      std::cout << "Volume fractions: beta, gamma, hcp: " << vfrac[0] << ", "
                << vfrac[1] << ", " << vfrac[2] << std::endl;
      std::cout << "Density: beta, gamma, hcp: \t\t" << rho[0] << ", "
                << rho[1] << ", " << rho[2] << ", (" << out_rho0[n] << ", "
                << out_rho1[n] << ", " << out_rho2[n] << ")" << std::endl;
      std::cout << "Pressure: beta, gamma, hcp: \t" << press[0] << ", "
                << press[1] << ", " << press[2] << ", (" << out_press[n]
                << ")" << std::endl;
      std::cout << "Temperature: beta, gamma, hcp: \t" << temp[0] << ", "
                << temp[1] << ", " << temp[2] << ", (" << out_temp[n] << ")"
                << std::endl;
      std::cout << "Internal energy: beta, gamma, hcp: \t" << sie[0] << ", "
                << sie[1] << ", " << sie[2] << ", (" << out_sie0[n] << ", "
                << out_sie1[n] << ", " << out_sie2[n] << ")" << std::endl
                << std::endl;

}

} //end namespace pte_test_3phase 

#endif // _SINGULARITY_EOS_TEST_PTE_TEST_3PHASESESAMESN_
