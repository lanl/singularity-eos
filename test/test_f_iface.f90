!------------------------------------------------------------------------------
! © 2021-2023. Triad National Security, LLC. All rights reserved.  This
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

program test_sg_fortran_interface
! modules
use singularity_eos_types
use singularity_eos
! no implicit vars
implicit none
! variable declaration
integer                                   :: nmat, res, mat
type(sg_eos_ary_t)                        :: eos

! set test parameters
nmat = 5

! allocate and initialize eos's
res = init_sg_eos_f(nmat, eos)
mat = 1

res = init_sg_IdealGas_f(mat, eos, 1.4d0, 1.0d7)

mat = mat + 1
res = init_sg_Gruneisen_f(mat, eos, 394000.d0, 1.489d0, 0.d0, 0.d0, 2.02d0, 0.47d0,&
                       8.93d0, 297.0d0, 1.0d6, 0.383d7)

mat = mat + 1
res = init_sg_JWL_f(mat, eos, 1.36177d13, 7.199d11, 6.2d0, 2.2d0, 0.5d0, 1.895d0, 1.0d7)

mat = mat + 1
res = init_sg_DavisReactants_f(mat, eos, 1.890d0, 4.115d10, 1.0d6, 297.0d0, 1.8d0,&
                            4.6d0, 0.34d0, 0.56d0,0.d0, 0.4265d0, 0.001074d10)
mat = mat + 1
res = init_sg_DavisProducts_f(mat, eos, 0.798311d0, 0.58d0, 1.35d0, 2.66182d0,&
                           0.75419d0, 3.2d10, 0.001072d10, 0.d0)

! cleanup
res = finalize_sg_eos_f(nmat, eos)

end program test_sg_fortran_interface
