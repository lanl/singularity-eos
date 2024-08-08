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

#ifndef _SINGULARITY_EOS_TEST_PTE_LONGTEST_2PHASEVINETSN_
#define _SINGULARITY_EOS_TEST_PTE_LONGTEST_2PHASEVINETSN_

#include <stdlib.h>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <singularity-eos/eos/eos.hpp>

// This data is taken from the .ult file for a test case in Flag, MPMat_KPT_PTsolv.ult.std
// Note that at this point in time (April 2024) the global and phase times are not alinged
// in .ult but they will be in the near future.

namespace pte_longtest_2phaseVinetSn {
constexpr int NMAT = 2;
constexpr int NTRIAL = 20;
constexpr int NPTS = NTRIAL * NMAT;
constexpr int HIST_SIZE = 10;

constexpr Real in_rho_tot[NTRIAL] = {
    7.278,      7.46221983, 7.65110201, 7.84476448, 8.04332821, 8.24691719, 8.45565857,
    8.66968272, 8.88912327, 9.11411728, 9.34480523, 9.58133117, 9.8238428,  10.0724915,
    10.3274327, 10.5888253, 10.8568327, 11.1316222, 11.4133653, 11.702238};
constexpr Real in_sie_tot[NTRIAL] = {
    8.41323509e8, 8.62513154e8,         9.28665206e8,  1.04377144e9,  1.21199344e9,
    1.43767572e9, 1.72535639e9,         2.07977629e9,  2.50588681e9,  3.00885704e9,
    3.59408011e9, 4.2671792100000005e9, 5.03401308e9,  5.90068122e9,  6.87352878e9,
    7.95915117e9, 9.16439846e9,         1.04963795e10, 1.19624656e10, 1.35702945e10};
constexpr Real in_lambda[NMAT] = {0.500480901, 0.499519099};
constexpr Real trial_vfrac[NMAT] = {47.0 / 100.0, 53.0 / 100.0};
constexpr Real out_press[NTRIAL] = {
    -3.29164604e-6,        1.284232715e10,        2.753234765e10,
    4.423652945000001e10,  6.313670939999999e10,  8.443021825e10,
    1.083303075e11,        1.3506674800000002e11, 1.64886557e11,
    1.9805482549999997e11, 2.3485562650000003e11, 2.755929975e11,
    3.20591986e11,         3.70199756e11,         4.247867485e11,
    4.84747905e11,         5.505039405000001e11,  6.225026735e11,
    7.012204135e11,        7.871634035e11};
constexpr Real out_temp[NTRIAL] = {298.,        317.657834,         338.05429449999997,
                                   359.18724,   381.0524905,        403.64396550000004,
                                   426.953795,  450.97240750000003, 475.6886215,
                                   501.0897275, 527.1615735,        553.8886454999999,
                                   581.254151,  609.2401025,        637.8273995,
                                   666.995915,  696.7245780000001,  726.9914615,
                                   757.7738645, 789.048397};
constexpr Real out_rho0[NTRIAL] = {
    7.285,      7.44383086, 7.61073863, 7.78523193, 7.9669718,  8.15572753, 8.35135008,
    8.55375288, 8.76289814, 8.97878664, 9.20145031, 9.43094663, 9.6673544,  9.91077059,
    10.1613078, 10.4190924, 10.6842631, 10.9569698, 11.2373726, 11.5256413};
constexpr Real out_rho1[NTRIAL] = {
    7.271,      7.48073556, 7.69197476, 7.90533181, 8.12131372, 8.34035068, 8.56281419,
    8.78903065, 9.01929178, 9.25386249, 9.49298694, 9.73689325, 9.98579716, 10.2399049,
    10.4994157, 10.7645232, 11.0354173, 11.3122855, 11.5953136, 11.8846866};
constexpr Real out_vfrac0[NTRIAL] = {0.5,         0.501717271, 0.50313519,  0.504308007,
                                     0.50527757,  0.506076807, 0.506731916, 0.507263967,
                                     0.507690077, 0.508024281, 0.508278194, 0.508461499,
                                     0.508582337, 0.508647597, 0.508663148, 0.50863402,
                                     0.508564548, 0.508458492, 0.508319122, 0.508149303};
constexpr Real out_vfrac1[NTRIAL] = {0.5,         0.498282729, 0.49686481,  0.495691993,
                                     0.49472243,  0.493923193, 0.493268084, 0.492736033,
                                     0.492309923, 0.491975719, 0.491721806, 0.491538501,
                                     0.491417663, 0.491352403, 0.491336852, 0.49136598,
                                     0.491435452, 0.491541508, 0.491680878, 0.491850697};

template <typename T>
inline void set_eos(T *eos) {
  constexpr Real d2to40[39] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                               0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                               0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

  singularity::EOS Snbeta =
      singularity::Vinet(7.285, 298.0, 0.529e12, 5.3345, 0.000072977, 0.2149e07, 0.658e09,
                         0.4419e07, d2to40);
  singularity::EOS Sngamma =
      singularity::Vinet(7.271, 298.0, 0.3878e12, 6.0532, 0.0001085405, 0.2161e07,
                         1.025e09, 0.5051e07, d2to40);
  eos[0] = Snbeta.GetOnDevice();
  eos[1] = Sngamma.GetOnDevice();
  return;
}

template <typename RealIndexer, typename EOSIndexer>
inline void set_trial_state(int n, RealIndexer &&rho, RealIndexer &&vfrac,
                            RealIndexer &&sie, EOSIndexer &&eos) {

  Real vsum = 0.;
  for (int i = 0; i < NMAT; i++) {
    vfrac[i] = trial_vfrac[i];
    rho[i] = in_lambda[i] * in_rho_tot[n] / vfrac[i];
    // same sie in both phases gives sie_tot=sie
    sie[i] = in_sie_tot[n];
    vsum += vfrac[i];
  }

  for (int i = 0; i < NMAT; i++)
    vfrac[i] *= 1.0 / vsum;

  return;
}

} // namespace pte_longtest_2phaseVinetSn

#endif // _SINGULARITY_EOS_TEST_PTE_LONGTEST_2PHASEVINETSN_
