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

template <typename T>
inline void set_eos(const int n,
                    // const int n, const std::array<std::size_t, sizeof(n)> &phases,
                    T *eos) {

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

  int phases[n] = {0, 1, 3};

  for (int i = 0; i < n; i++) {
    eos[i] = alleos[phases[i]];
  }

  return;
}

#endif // _SINGULARITY_EOS_TEST_MATERIAL_5PHASESESAMESN_
