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

#ifndef _SINGULARITY_EOS_TEST_PTE_TEST_UTILS_
#define _SINGULARITY_EOS_TEST_PTE_TEST_UTILS_

#include <stdlib.h>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <singularity-eos/eos/eos.hpp>

constexpr int NMAT = 3;
constexpr int NTRIAL = 100;
constexpr int NPTS = NTRIAL * NMAT;
constexpr int HIST_SIZE = 10;

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
  singularity::EOS gr = singularity::Gruneisen(394000.0, 1.489, 0.0, 0.0, 2.02, 0.47,
                                               8.93, 297.0, 1.0e6, 0.383e7);
  singularity::EOS dr = singularity::DavisReactants(
      1.890, 4.115e10, 1.0e6, 297.0, 1.8e5, 4.6, 0.34, 0.56, 0.0, 0.4265, 0.001074e10);
  singularity::EOS dp = singularity::DavisProducts(0.798311, 0.58, 1.35, 2.66182, 0.75419,
                                                   3.2e10, 0.001072e10, 0.0);
  eos[0] = gr.GetOnDevice();
  eos[1] = dr.GetOnDevice();
  eos[2] = dp.GetOnDevice();
  return;
}

template <typename RealIndexer, typename EOSIndexer>
inline void set_state(RealIndexer &&rho, RealIndexer &&vfrac, RealIndexer &&sie,
                      RealIndexer &&temp, EOSIndexer &&eos) {

  rho[0] = 8.93;
  rho[1] = 1.89;
  rho[2] = 2.5;

  Real vsum = 0.;
  for (int i = 0; i < NMAT; i++) {
    temp[i] = 600.0;
    sie[i] = eos[i].InternalEnergyFromDensityTemperature(rho[i], temp[i]);
    vfrac[i] = rand() / (1.0 * RAND_MAX);
    vsum += vfrac[i];
  }

  for (int i = 0; i < NMAT; i++)
    vfrac[i] *= 1.0 / vsum;

  return;
}

#endif // _SINGULARITY_EOS_TEST_PTE_TEST_UTILS_
