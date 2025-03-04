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

#include <cassert>
#include <ports-of-call/portability.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_builder.hpp>
#include <singularity-eos/eos/singularity_eos.hpp>
#include <singularity-eos/eos/singularity_eos_init_utils.hpp>

using namespace singularity;

int get_sg_PressureFromDensityInternalEnergy(int matindex, EOS *eos, const double *rhos,
                                             const double *sies, double *pressures,
                                             const int len) {
  eos[matindex].PressureFromDensityInternalEnergy(rhos, sies, pressures, len);
  return 0;
}
int get_sg_MinInternalEnergyFromDensity(int matindex, EOS *eos, const double *rhos,
                                        double *sies, const int len) {
  eos[matindex].MinInternalEnergyFromDensity(rhos, sies, len);
  return 0;
}
int get_sg_BulkModulusFromDensityInternalEnergy(int matindex, EOS *eos,
                                                const double *rhos, const double *sies,
                                                double *bmods, const int len) {
  eos[matindex].BulkModulusFromDensityInternalEnergy(rhos, sies, bmods, len);
  return 0;
}
