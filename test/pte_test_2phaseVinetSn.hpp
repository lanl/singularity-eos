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

#ifndef _SINGULARITY_EOS_TEST_PTE_TEST_2PHASEVINETSN_
#define _SINGULARITY_EOS_TEST_PTE_TEST_2PHASEVINETSN_

#include <stdlib.h>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <singularity-eos/eos/eos.hpp>

constexpr int NMAT = 2;
constexpr int NTRIAL = 5;
constexpr int NPTS = NTRIAL * NMAT;
constexpr int HIST_SIZE = 10;

constexpr Real in_rho_tot[5] = {7.27800000, 7.27981950, 7.28163945, 7.28345986, 7.28528073};
constexpr Real in_sie_tot[5] = {8.41323509e08, 8.41325565e08, 8.41331734e08, 8.41342022e08, 8.41356432e08};
constexpr Real in_lambda0 = 0.500480901;
constexpr Real in_lambda1 = 0.499519099;

constexpr Real trial_vfrac0 = 47.0/100.0;
constexpr Real trial_vfrac1 = 53.0/100.0;

constexpr Real out_press[5] = {-3.29164604e-6,1.19722694e8,2.3968114450000003e8,3.598077865e8,4.80104422e8};
constexpr Real out_temp[5] = {298.,298.192952,298.38594450000005,298.5790125,298.7721535};
constexpr Real out_rho0[5] = {7.28500000, 7.28653963, 7.28808705, 7.28963516, 7.29118416};
constexpr Real out_rho1[5] = {7.27100000, 7.27309885, 7.27519088, 7.27728316, 7.27937551};
constexpr Real out_vfrac0[5] = {0.5, 0.500019325, 0.500038138, 0.500056927, 0.500075679};
constexpr Real out_vfrac1[5] = {0.5, 0.499980675, 0.499961862, 0.499943073, 0.499924321};

template <typename T>
class LinearIndexer {
 public:
  PORTABLE_FUNCTION LinearIndexer() = default;
  LinearIndexer(const T &t) : data_(t) {}
  PORTABLE_INLINE_FUNCTION
  auto &operator[](const int i) const { return data_(i); }

 private:
  T data_;
};

template <typename T>
class Indexer2D {
 public:
  PORTABLE_FUNCTION Indexer2D() = default;
  PORTABLE_FUNCTION Indexer2D(const int j, const T &t) : j_(j), data_(t) {}
  Indexer2D(const int j, const T &&t) = delete; // prevents r-value binding
  PORTABLE_INLINE_FUNCTION
  auto &operator[](const int i) const { return data_(j_, i); }

 private:
  const int j_;
  const T &data_;
};

template <typename T>
inline void set_eos(T *eos) {
  constexpr Real d2to40[39] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

  singularity::EOS Snbeta = singularity::Vinet(7.285, 298.0, 0.529e12, 5.3345, 0.000072977, 0.2149e07,
                                               0.658e09, 0.4419e07, d2to40);
  singularity::EOS Sngamma = singularity::Vinet(7.271, 298.0, 0.3878e12, 6.0532, 0.0001085405, 0.2161e07,
                                               1.025e09, 0.5051e07, d2to40);
  eos[0] = Snbeta.GetOnDevice();
  eos[1] = Sngamma.GetOnDevice();
  return;
}

template <typename RealIndexer, typename EOSIndexer>
inline void set_state(RealIndexer &&rho, RealIndexer &&vfrac, RealIndexer &&sie,
                      RealIndexer &&temp, EOSIndexer &&eos) {

  vfrac[0] = 47.0/100.0;
  vfrac[1] = 53.0/100.0;
  rho[0] = 7.278 * 0.500480901 / vfrac[0];
  rho[1] = 7.278 * 0.499519099 / vfrac[1];
  sie[0] = 8.41323509e08;
  sie[1] = 8.41323509e08;


  return;
}

#endif // _SINGULARITY_EOS_TEST_PTE_TEST_2PHASEVINETSN_
