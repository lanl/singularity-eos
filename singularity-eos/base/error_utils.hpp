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

#ifndef SINGULARITY_EOS_BASE_EOS_ERROR_HPP_
#define SINGULARITY_EOS_BASE_EOS_ERROR_HPP_

#include <ports-of-call/portable_errors.hpp>
#define EOS_ERROR(x) PORTABLE_ABORT(x)
#define UNDEFINED_ERROR EOS_ERROR("DEFINE ME\n")

namespace singularity {
namespace error_utils {
constexpr std::size_t MAX_NUM_CHARS = 121;
// Cuda doesn't have strcat, so we implement it ourselves
PORTABLE_FORCEINLINE_FUNCTION
char *StrCat(char *destination, const char *source) {
  int i, j; // not in loops because they're re-used.

  // specifically avoid strlen, which isn't on GPU
  for (i = 0; destination[i] != '\0'; i++) {
  }
  // assumes destination has enough memory allocated
  for (j = 0; source[j] != '\0'; j++) {
    // MAX_NUM_CHARS-1 to leave room for null terminator
    PORTABLE_REQUIRE((i + j) < MAX_NUM_CHARS - 1,
                     "Concat string must be within allowed size");
    destination[i + j] = source[j];
  }
  // null terminate destination string
  destination[i + j] = '\0';

  // the destination is returned by standard `strcat()`
  return destination;
}
} // namespace error_util
} // namespace singularity


#endif // SINGULARITY_EOS_BASE_EOS_ERROR_HPP_
