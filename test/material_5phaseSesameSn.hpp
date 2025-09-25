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

#ifndef _SINGULARITY_EOS_TEST_MATERIAL_5PHASESESAMESN_
#define _SINGULARITY_EOS_TEST_MATERIAL_5PHASESESAMESN_

using EOS = singularity::Variant<singularity::EOSPAC>;

constexpr Real w[25] = {0.,  0.8510E+01, 20., 20., 20., 0.8510E+01, 0.,  20., 20.,
                        20., 20.,        20., 0.,  20., 20.,        20., 20., 20.,
                        0.,  20.,        20., 20., 20., 20.,        0.};

constexpr Real b[25] = {0.,         0.3060E-04, 0.1000E-05, 0.1000E-05, 0.1000E-05,
                        0.3060E-04, 0.,         0.1000E-05, 0.1000E-05, 0.1000E-05,
                        0.1000E-05, 0.1000E-05, 0.,         0.1000E-05, 0.1000E-05,
                        0.1000E-05, 0.1000E-05, 0.1000E-05, 0.,         0.1000E-05,
                        0.1000E-05, 0.1000E-05, 0.1000E-05, 0.1000E-05, 0.};
template <typename T>
inline void set_mfupdate_params(const int n, const std::vector<int> &nphases, T *wparam, T *bparam) {

  int numphases = 5;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      wparam[i*n + j] = w[nphases[i]*numphases+nphases[j]];
      bparam[i*n + j] = b[nphases[i]*numphases+nphases[j]];
    }
  }

  return;
}

template <typename T>
inline void set_eos(const int n, const std::vector<int> &nphases, T *eos) {

  int sl = symlink("/usr/projects/data/eos/eos-developmental/Sn2162/v01/sn2162-v01.bin",
                   "sesameu");

  int numphases = 5;

  int SnbetaID = 2102;
  int SngammaID = 2103;
  int SndeltaID = 2104;
  int SnhcpID = 2105;
  int SnliquidID = 2106;

  // bool invert_at_setup = true;

  EOS Snbeta = singularity::EOSPAC(SnbetaID);
  EOS Sngamma = singularity::EOSPAC(SngammaID);
  EOS Sndelta = singularity::EOSPAC(SndeltaID);
  EOS Snhcp = singularity::EOSPAC(SnhcpID);
  EOS Snliquid = singularity::EOSPAC(SnliquidID);

  std::vector<EOS> alleos(numphases);

  alleos[0] = Snbeta.GetOnDevice();
  alleos[1] = Sngamma.GetOnDevice();
  alleos[2] = Sndelta.GetOnDevice();
  alleos[3] = Snhcp.GetOnDevice();
  alleos[4] = Snliquid.GetOnDevice();

  for (int i = 0; i < n; i++) {
    eos[i] = alleos[nphases[i]];
  }

  return;
}

#endif // _SINGULARITY_EOS_TEST_MATERIAL_5PHASESESAMESN_
