!------------------------------------------------------------------------------
! © 2021-2022. Triad National Security, LLC. All rights reserved.  This
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

program test_ftniface
use iso_c_binding
use singularity_eos
use singularity_eos_types
implicit none

integer                                   :: res, nmat, ncell, i, j, matid
type(sg_eos_ary_t)                        :: eos

integer,      allocatable, dimension(:)   :: offsets
integer,      dimension(6)                :: sg_mods_enabled
real(kind=8), dimension(4)                :: sg_mods_values
real(kind=8), allocatable, dimension(:)   :: spvol, sie_tot, vsum, press, pmax
real(kind=8), allocatable, dimension(:)   :: temp, bmod, dpde, cv
real(kind=8), allocatable, dimension(:,:) :: frac_mass, frac_vol, frac_sie
real(kind=8)                              :: vsumint
character(len=64)                         :: filename

ncell = 1
nmat = 3
!filename = "../data/materials.sp5"

sg_mods_enabled = 0
sg_mods_values = 0

! allocate arrays
allocate(spvol(ncell), sie_tot(ncell), temp(ncell), press(ncell), vsum(ncell))
allocate(bmod(ncell), dpde(ncell), cv(ncell), pmax(ncell))
allocate(offsets(ncell))
allocate(frac_mass(ncell, nmat))
allocate(frac_vol(ncell, nmat), frac_sie(ncell, nmat))
! initialize eos's
res = init_sg_eos_f(nmat, eos)

res = init_sg_Gruneisen_f(1, eos, 394000.d0, 1.489d0, 0.d0, 0.d0, 2.02d0, 0.47d0,&
                       8.93d0, 297.0d0, 1.0d6, 0.383d7, sg_mods_enabled, sg_mods_values)
res = init_sg_DavisReactants_f(2, eos, 1.890d0, 4.115d10, 1.0d6, 297.0d0, 1.8d0,&
                            4.6d0, 0.34d0, 0.56d0,0.d0, 0.4265d0, 0.001074d10, sg_mods_enabled, sg_mods_values)
res = init_sg_DavisProducts_f(3, eos, 0.798311d0, 0.58d0, 1.35d0, 2.66182d0,&
                           0.75419d0, 3.2d10, 0.001072d10, 0.d0, sg_mods_enabled, sg_mods_values)

!TODO: run some evals

res = finalize_sg_eos_f(nmat, eos)

deallocate(vsum, sie_tot, spvol, frac_mass, temp, press, frac_vol, frac_sie)
deallocate(bmod, dpde, cv, offsets, pmax)
end program test_ftniface
