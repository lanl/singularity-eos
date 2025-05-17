//------------------------------------------------------------------------------
// © 2021-2025. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_TEST_MATERIAL_2PHASEVINETSN_
#define _SINGULARITY_EOS_TEST_MATERIAL_2PHASEVINETSN_

using EOS = singularity::Variant<singularity::Vinet>;

template <typename T>
inline void set_eos(const int n, const std::vector<int> &nphases, T *eos) {

  int numphases = 2;

  constexpr Real d2to40[39] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                               0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                               0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

  EOS Snbeta = singularity::Vinet(7.285, 298.0, 0.529e12, 5.3345, 0.000072977, 0.2149e07,
                                  0.658e09, 0.4419e07, d2to40);
  EOS Sngamma = singularity::Vinet(7.271, 298.0, 0.3878e12, 6.0532, 0.0001085405,
                                   0.2161e07, 1.025e09, 0.5051e07, d2to40);
  std::vector<EOS> alleos(numphases);

  alleos[0] = Snbeta.GetOnDevice();
  alleos[1] = Sngamma.GetOnDevice();

  for (int i = 0; i < n; i++) {
    eos[i] = alleos[nphases[i]];
  }

  return;
}

#endif // _SINGULARITY_EOS_TEST_MATERIAL_2PHASEVINETSN_
