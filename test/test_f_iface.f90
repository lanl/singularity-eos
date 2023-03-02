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

program test
use bedroom_door_eos
integer                                   :: res, nmat, ncell, i, j, matid
type(bd_eos_ary_t)                        :: eos
integer,      allocatable, dimension(:)   :: offsets
real(kind=8), allocatable, dimension(:)   :: spvol, sie_tot, vsum, press, pmax
real(kind=8), allocatable, dimension(:)   :: temp, bmod, dpde, cv
real(kind=8), allocatable, dimension(:,:) :: frac_mass, frac_vol, frac_sie
real(kind=8)                              :: vsumint
character(len=64)                         :: filename

ncell = 1
nmat = 3
filename = "../data/materials.sp5"

! allocate arrays
allocate(spvol(ncell), sie_tot(ncell), temp(ncell), press(ncell), vsum(ncell))
allocate(bmod(ncell), dpde(ncell), cv(ncell), pmax(ncell))
allocate(offsets(ncell))
allocate(frac_mass(ncell, nmat))
allocate(frac_vol(ncell, nmat), frac_sie(ncell, nmat))
! initialize eos's
res = init_bd_eos_f(nmat, eos)

res = init_Gruneisen_f(1, eos, 394000.d0, 1.489d0, 0.d0, 0.d0, 2.02d0, 0.47d0,&
                       8.93d0, 297.0d0, 1.0d6, 0.383d7)
res = init_DavisReactants_f(2, eos, 1.890d0, 4.115d10, 1.0d6, 297.0d0, 1.8d0,&
                            4.6d0, 0.34d0, 0.56d0,0.d0, 0.4265d0, 0.001074d10)
res = init_DavisProducts_f(3, eos, 0.798311d0, 0.58d0, 1.35d0, 2.66182d0,&
                           0.75419d0, 3.2d10, 0.001072d10, 0.d0)
!matid = 4272
!res = init_SpinerDependsRhoT_f(2, eos, filename, matid)
!matid = 5030
!res = init_SpinerDependsRhoT_f(3, eos, filename, matid)
! initialize state variables
call srand(10)
press = 0.d0
do i=1,ncell
  vsum(i) = 1.0d-2
  frac_mass(i,1) = 8.93d0*vsum(i)*0.33
  frac_mass(i,2) = 1.89d0*vsum(i)*.033
  frac_mass(i,3) = 2.5d0*vsum(i)*0.34
  sie_tot(i) = 0.d0
  spvol(i) = 0.d0
  !vsum(i) = 0.d0
  offsets(i) = i
  frac_sie(i,1) = 1.16049000000000000d+09
  frac_sie(i,2) = 4.50111706015744858d+10
  frac_sie(i,3) = 3.29477034927098885d+10
  temp(i) = 0.d0
  ! frac_vols will be normalized in eos call
  do j=1,nmat
    frac_vol(i,j) = 0.d0
  enddo
  do j=1,nmat
    spvol(i) = spvol(i) + frac_mass(i,j)
    sie_tot(i) = sie_tot(i) + frac_mass(i,j)*frac_sie(i,j)
    frac_sie(i,j) = frac_sie(i,j)*frac_mass(i,j) ! needs to be energy
  enddo
  sie_tot(i) = sie_tot(i)/spvol(i)
  spvol(i) = vsum(i)/spvol(i)
enddo

!write(*,*) 'spvol', spvol
write(*,*) frac_mass, frac_vol, frac_sie
! calculate PTE eos
res = get_bd_eos_f(nmat, ncell, ncell, 0, eos, offsets, press, pmax, vsum,&
                   spvol, sie_tot, temp, bmod, dpde, cv, frac_mass, frac_vol,&
                   frac_sie)

write(*,*) press, cv, bmod, temp
write(*,*) frac_mass, frac_vol, frac_sie
! cleanup
res = finalize_bd_eos_f(nmat, eos)

deallocate(vsum, sie_tot, spvol, frac_mass, temp, press, frac_vol, frac_sie)
deallocate(bmod, dpde, cv, offsets, pmax)
end program test
