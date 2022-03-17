#!/bin/bash

#------------------------------------------------------------------------------
# © 2022. Triad National Security, LLC. All rights reserved.  This
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

# get version number
A=$(grep "singularity-eos VERSION" ../../CMakeLists.txt)
for word in ${A}
do 
  B=$(echo ${word} | tr -d ')')
done

# go out 3 levels
cd ../../../
tar czvf singularity-eos-${B}.tar.gz --exclude-vcs --exclude=".git*" \
--exclude="build*" --exclude="nbproject" --exclude="data" \
--exclude="eigen/doc" --exclude="eigen/demos" --exclude="eigen/failtest" \
--exclude="eigen/test" --exclude="eigen/bench" --exclude="eigen/unsupported" \
--exclude="json/benchmarks" --exclude="json/doc" --exclude="variant/3rdparty" \
--exclude="variant/test" --exclude="utils/kokkos" --exclude="test" \
--exclude="Catch2" --exclude="singularity-eos-data" singularity-eos/
