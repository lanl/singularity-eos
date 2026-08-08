#!/bin/bash

#------------------------------------------------------------------------------
# © 2021-2025. Triad National Security, LLC. All rights reserved.  This
# program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by Triad
# National Security, LLC for the U.S.  Department of Energy/National
# Nuclear Security Administration. All rights in the program are
# reserved by Triad National Security, LLC, and the U.S. Department of
# Energy/National Nuclear Security Administration. The Government is
# granted for itself and others acting on its behalf a nonexclusive,
# paid-up, irrevocable worldwide license in this material to reproduce,
# prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.
#------------------------------------------------------------------------------

# Report singularity-eos headers that are missing install logic.
#
# Every .hpp under singularity-eos/ is meant to be installed, which means it
# must appear as an argument to a register_headers(...) call in some
# CMakeLists.txt (see cmake/plugins.cmake). A header that is never registered
# is silently left out of the install. This script flags any such header so
# it can be wired into CI as a guard (non-zero exit when any are found).

set -euo pipefail

# Move to the repo root (two levels up from utils/scripts/).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

# Directory whose headers must all be installed.
HEADER_ROOT="singularity-eos"

# Collect the text of every CMakeLists.txt so we can search registrations.
CMAKE_FILES=$(find . -name CMakeLists.txt -not -path './build/*' -not -path './.git/*')

missing=0
for f in $(find "${HEADER_ROOT}" -name '*.hpp' | sed "s|^${HEADER_ROOT}/||" | sort); do
  base=$(basename "$f")
  # A header is considered registered if either its path relative to the
  # header root, or just its file name, appears in a CMakeLists.txt. The
  # register_headers() paths are relative to the calling CMakeLists.txt, so
  # matching on the basename keeps this robust without re-implementing
  # CMake's path resolution.
  if ! grep -rqF "$f" ${CMAKE_FILES} && ! grep -rqF "$base" ${CMAKE_FILES}; then
    if [ "$missing" -eq 0 ]; then
      echo "Headers missing register_headers() install logic:"
    fi
    echo "  ${HEADER_ROOT}/$f"
    missing=$((missing + 1))
  fi
done

if [ "$missing" -eq 0 ]; then
  echo "OK: every ${HEADER_ROOT} header has register_headers() install logic."
  exit 0
fi

echo ""
echo "${missing} header(s) have no register_headers call."
echo "Add each to a register_headers(...) call in the appropriate CMakeLists.txt"
echo "(with the matching build condition if it is only needed under an option),"
echo "or exclude it deliberately."
exit 1
