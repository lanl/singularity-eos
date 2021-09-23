//------------------------------------------------------------------------------
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

#ifndef SINGULARITY_EOS_BASE_EOS_ERROR_HPP_
#define SINGULARITY_EOS_BASE_EOS_ERROR_HPP_

#ifdef SINGULARITY_ENABLE_EXCEPTIONS
#include <stdexcept>
#define EOS_ERROR(x) (throw std::runtime_error(x))
#else
#include <cassert> // got errors about untyped object in template
#define EOS_ERROR(x)                                                                     \
  printf("%s\n", x);                                                                     \
  assert(false);
#endif
#define UNDEFINED_ERROR EOS_ERROR("DEFINE ME\n")

#endif // SINGULARITY_EOS_BASE_EOS_ERROR_HPP_
