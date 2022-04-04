!------------------------------------------------------------------------------
! © 2021. Triad National Security, LLC. All rights reserved.  This
! program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by Triad
! National Security, LLC for the U.S.  Department of Energy/National
! Nuclear Security Administration. All rights in the program are
! reserved by Triad National Security, LLC, and the U.S. Department of
! Energy/National Nuclear Security Administration. The Government is
! granted for itself and others acting on its behalf a nonexclusive,
! paid-up, irrevocable worldwide license in this material to reproduce,
! prepare derivative works, distribute copies to the public, perform
! publicly and display publicly, and to permit others to do so.
!------------------------------------------------------------------------------

program test
use singularity_eos
integer                                   :: res, nmat, m
type(sg_eos_ary_t)                        :: eos
real :: Cv, gm1
nmat=2
m=1
Cv=2
gm1=1.4-1
res = init_sg_eos_f(nmat, eos)
res_sg = init_sg_IdealGas_f(m, eos, gm1, Cv)
end program test
