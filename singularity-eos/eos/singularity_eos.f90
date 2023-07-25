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

module singularity_eos_types
  use iso_c_binding
  implicit none

  ! data types
  public :: sg_eos_ary_t

  type :: sg_eos_ary_t
    type(c_ptr) :: ptr = C_NULL_PTR
  end type sg_eos_ary_t
end module singularity_eos_types

module singularity_eos
  use iso_c_binding
  use singularity_eos_types
  implicit none

! fortran functions that call  the interfaces
  public :: &
    init_sg_eos_f,&
    init_sg_IdealGas_f,&
    init_sg_Gruneisen_f,&
    init_sg_JWL_f,&
    init_sg_DavisProducts_f,&
    init_sg_DavisReactants_f,&
    init_sg_eospac_f,&
    get_sg_PressureFromDensityInternalEnergy_f,&
    get_sg_MinInternalEnergyFromDensity_f,&
    get_sg_BulkModulusFromDensityInternalEnergy_f,&
    get_sg_eos_f,&
    finalize_sg_eos_f

! interface functions
  interface
    integer(kind=c_int) function &
      init_sg_eos(nmat, eos) &
      bind(C, name='init_sg_eos')
      import
      integer(c_int), value, intent(in) :: nmat
      type(c_ptr), intent(in)           :: eos
    end function init_sg_eos
  end interface

  interface
    integer(kind=c_int) function &
      init_sg_IdealGas(matindex, eos, gm1, Cv, sg_mods_enabled, &
                       sg_mods_values) &
      bind(C, name='init_sg_IdealGas')
      import
      integer(c_int), value, intent(in)      :: matindex
      type(c_ptr), value, intent(in)         :: eos
      real(kind=c_double), value, intent(in) :: gm1, Cv
      type(c_ptr), value, intent(in)         :: sg_mods_enabled, sg_mods_values
    end function init_sg_IdealGas
  end interface

  interface
    integer(kind=c_int) function &
      init_sg_Gruneisen(matindex, eos, C0, s1, s2, s3, G0, b, rho0, T0, P0,&
                        Cv, sg_mods_enabled, sg_mods_values) &
      bind(C, name='init_sg_Gruneisen')
      import
      integer(c_int), value, intent(in)      :: matindex
      type(c_ptr), value, intent(in)         :: eos
      real(kind=c_double), value, intent(in) :: C0, s1, s2, s3, G0, b, rho0
      real(kind=c_double), value, intent(in) :: T0, P0, Cv
      type(c_ptr), value, intent(in)         :: sg_mods_enabled, sg_mods_values
    end function init_sg_Gruneisen
  end interface
  
  interface
    integer(kind=c_int) function &
      init_sg_JWL(matindex, eos, A, B, R1, R2, w, rho0, Cv, sg_mods_enabled, &
                  sg_mods_values) &
      bind(C, name='init_sg_JWL')
      import
      integer(c_int), value, intent(in)      :: matindex
      type(c_ptr), value, intent(in)         :: eos
      real(kind=c_double), value, intent(in) :: A, B, R1, R2, w, rho0, Cv
      type(c_ptr), value, intent(in)         :: sg_mods_enabled, sg_mods_values
    end function init_sg_JWL
  end interface

  interface
    integer(kind=c_int) function &
      init_sg_DavisProducts(matindex, eos, a, b, k, n, vc, pc, Cv, E0, &
                            sg_mods_enabled, sg_mods_values) &
      bind(C, name='init_sg_DavisProducts')
      import
      integer(c_int), value, intent(in)      :: matindex
      type(c_ptr), value, intent(in)         :: eos
      real(kind=c_double), value, intent(in) :: a, b, k, n, vc, pc, Cv, E0
      type(c_ptr), value, intent(in)         :: sg_mods_enabled, sg_mods_values
    end function init_sg_DavisProducts
  end interface

  interface
    integer(kind=c_int) function &
      init_sg_DavisReactants(matindex, eos, rho0, e0, P0, T0, A, B, C, G0, Z,&
                             alpha, Cv0, sg_mods_enabled, sg_mods_values) &
      bind(C, name='init_sg_DavisReactants')
      import
      integer(c_int), value, intent(in)      :: matindex
      type(c_ptr), value, intent(in)         :: eos
      real(kind=c_double), value, intent(in) :: rho0, e0, P0, T0, A, B, C, G0,&
                                                Z, alpha, Cv0
      type(c_ptr), value, intent(in)         :: sg_mods_enabled, sg_mods_values
    end function init_sg_DavisReactants
  end interface

  interface
    integer(kind=c_int) function &
      init_sg_eospac(matindex, eos, id, eospac_opts_values, &
      sg_mods_enabled, sg_mods_values) &
      bind(C, name='init_sg_eospac')
      import
      integer(c_int), value, intent(in) :: matindex, id
      type(c_ptr), value, intent(in)    :: eos
      type(c_ptr), value, intent(in)    :: sg_mods_enabled, sg_mods_values
      type(c_ptr), value, intent(in)    :: eospac_opts_values
    end function init_sg_eospac
  end interface

  interface
   integer(kind=c_int) function &
       get_sg_PressureFromDensityInternalEnergy(matindex, eos, rhos, sies,&
                                               pressures, len) &
       bind(C, name='get_sg_PressureFromDensityInternalEnergy')
       import
       integer(c_int), value, intent(in) :: matindex, len
       type(c_ptr), value, intent(in) :: eos, rhos, sies
       type(c_ptr), value, intent(in) :: pressures
    end function
  end interface

  interface
   integer(kind=c_int) function &
       get_sg_MinInternalEnergyFromDensity(matindex, eos, rhos, sies,&
                                           len) &
       bind(C, name='get_sg_MinInternalEnergyFromDensity')
       import
       integer(c_int), value, intent(in) :: matindex, len
       type(c_ptr), value, intent(in) :: eos, rhos, sies
    end function
  end interface


  interface
     integer(kind=c_int) function &
       get_sg_BulkModulusFromDensityInternalEnergy(matindex, eos, rhos, sies,&
                                               bmods, len) &
       bind(C, name='get_sg_BulkModulusFromDensityInternalEnergy')
       import
       integer(c_int),value, intent(in) :: matindex, len
       type(c_ptr), value, intent(in) :: eos, rhos, sies
       type(c_ptr), value, intent(in) :: bmods
    end function
  end interface

  interface
    integer(kind=c_int) function &
      get_sg_eos(nmat, ncell, cell_dim,&
                 option,&
                 eos_offsets,&
                 eos,&
                 offsets,&
                 press, pmax, vol, spvol, sie, temp, bmod, dpde, cv,&
                 frac_mass, frac_vol, frac_sie,&
                 frac_bmod, frac_dpde, frac_cv)&
      bind(C, name='get_sg_eos')
      import
      integer(kind=c_int), value, intent(in) :: nmat
      integer(kind=c_int), value, intent(in) :: ncell
      integer(kind=c_int), value, intent(in) :: cell_dim
      integer(kind=c_int), value, intent(in) :: option
      type(c_ptr), value, intent(in) :: eos_offsets
      ! better eos ptrs
      type(c_ptr), value, intent(in) :: eos
      ! other inputs
      type(c_ptr), value, intent(in) :: offsets
      type(c_ptr), value, intent(in) :: press
      type(c_ptr), value, intent(in) :: pmax
      type(c_ptr), value, intent(in) :: vol
      type(c_ptr), value, intent(in) :: spvol
      type(c_ptr), value, intent(in) :: sie
      type(c_ptr), value, intent(in) :: temp
      type(c_ptr), value, intent(in) :: bmod
      type(c_ptr), value, intent(in) :: dpde
      type(c_ptr), value, intent(in) :: cv
      type(c_ptr), value, intent(in) :: frac_mass
      type(c_ptr), value, intent(in) :: frac_vol
      type(c_ptr), value, intent(in) :: frac_sie
      type(c_ptr), value, intent(in) :: frac_bmod
      type(c_ptr), value, intent(in) :: frac_dpde
      type(c_ptr), value, intent(in) :: frac_cv
    end function get_sg_eos
  end interface

  interface
    integer(kind=c_int) function &
      finalize_sg_eos(nmat, eos, own_kokkos) &
      bind(C, name='finalize_sg_eos')
      import
      integer(kind=c_int), value, intent(in) :: nmat
      type(c_ptr), intent(in)                :: eos
      integer(kind=c_int), value, intent(in) :: own_kokkos
    end function finalize_sg_eos
  end interface

! fortran functions
contains

  integer function get_sg_eos_f(nmat, ncell, cell_dim,&
                                option,&
                                eos_offsets,&
                                eos,&
                                offsets,&
                                press, pmax, vol, spvol, sie, temp, bmod,&
                                dpde, cv,&
                                frac_mass, frac_vol, frac_sie,&
                                frac_bmod, frac_dpde, frac_cv) &
    result(err)
    integer(kind=c_int), intent(in) :: nmat
    integer(kind=c_int), intent(in) :: ncell
    integer(kind=c_int), intent(in) :: cell_dim
    integer(kind=c_int), intent(in) :: option
    integer(kind=c_int), dimension(:), target, intent(in) :: eos_offsets
    type(sg_eos_ary_t), intent(in)  :: eos
    integer(kind=c_int), dimension(:), target, intent(in) :: offsets
    real(kind=8), dimension(:),   target, intent(in)    :: press
    real(kind=8), dimension(:),   target, intent(in)    :: pmax
    real(kind=8), dimension(:),   target, intent(in)    :: vol
    real(kind=8), dimension(:),   target, intent(in)    :: spvol
    real(kind=8), dimension(:),   target, intent(in)    :: sie
    real(kind=8), dimension(:),   target, intent(in)    :: temp
    real(kind=8), dimension(:),   target, intent(in)    :: bmod
    real(kind=8), dimension(:),   target, intent(in)    :: dpde
    real(kind=8), dimension(:),   target, intent(in)    :: cv
    real(kind=8), dimension(:,:), target, intent(in) :: frac_mass
    real(kind=8), dimension(:,:), target, intent(inout) :: frac_vol
    real(kind=8), dimension(:,:), target, intent(inout) :: frac_sie
    ! optionals
    real(kind=8), dimension(:,:), target, optional, intent(inout) :: frac_bmod
    real(kind=8), dimension(:,:), target, optional, intent(inout) :: frac_dpde
    real(kind=8), dimension(:,:), target, optional, intent(inout) :: frac_cv

    ! pointers
    type(c_ptr) :: bmod_ptr, dpde_ptr, cv_ptr

    bmod_ptr = C_NULL_PTR
    dpde_ptr = C_NULL_PTR
    cv_ptr = C_NULL_PTR
    if(present(frac_bmod)) then
      bmod_ptr = c_loc(frac_bmod)
    endif
    if(present(frac_dpde)) then
      dpde_ptr = c_loc(frac_dpde)
    endif
    if(present(frac_cv)) then
      cv_ptr = c_loc(frac_cv)
    endif

    err = get_sg_eos(nmat, ncell, cell_dim, option, c_loc(eos_offsets),&
                     eos%ptr, c_loc(offsets), c_loc(press), c_loc(pmax),&
                     c_loc(vol), c_loc(spvol), c_loc(sie), c_loc(temp),&
                     c_loc(bmod), c_loc(dpde),c_loc(cv), c_loc(frac_mass),&
                     c_loc(frac_vol),c_loc(frac_sie), bmod_ptr, dpde_ptr,&
                     cv_ptr)
  end function get_sg_eos_f

  integer function init_sg_eos_f(nmat, eos) &
    result(err)
    integer(kind=c_int), intent(in) :: nmat
    type(sg_eos_ary_t), intent(in)  :: eos
    err = init_sg_eos(nmat, eos%ptr)
  end function init_sg_eos_f

  integer function init_sg_IdealGas_f(matindex, eos, gm1, Cv, &
                                      sg_mods_enabled, sg_mods_values) &
    result(err)
    integer(c_int), value, intent(in) :: matindex
    type(sg_eos_ary_t), intent(in)    :: eos
    real(kind=8), value, intent(in)   :: gm1, Cv
    integer(kind=c_int), dimension(:), target, intent(inout) :: sg_mods_enabled
    real(kind=8), dimension(:), target, intent(inout)        :: sg_mods_values
    err = init_sg_IdealGas(matindex-1, eos%ptr, gm1, Cv, &
                           c_loc(sg_mods_enabled), c_loc(sg_mods_values))
  end function init_sg_IdealGas_f

  integer function init_sg_Gruneisen_f(matindex, eos, C0, s1, s2, s3, G0, b,&
                                       rho0, T0, P0, Cv, sg_mods_enabled, &
                                       sg_mods_values) &
    result(err)
    integer(c_int), value, intent(in) :: matindex
    type(sg_eos_ary_t), intent(in)    :: eos
    real(kind=8), value, intent(in)   :: C0, s1, s2, s3, G0, b, rho0
    real(kind=8), value, intent(in)   :: T0, P0, Cv
    integer(kind=c_int), dimension(:), target, intent(inout) :: sg_mods_enabled
    real(kind=8), dimension(:), target, intent(inout)        :: sg_mods_values
    err = init_sg_Gruneisen(matindex-1, eos%ptr, C0, s1, s2, s3, G0, b, rho0,&
                            T0, P0, Cv, c_loc(sg_mods_enabled), &
                            c_loc(sg_mods_values))
  end function init_sg_Gruneisen_f

  integer function init_sg_JWL_f(matindex, eos, A, B, R1, R2, w, rho0, Cv, &
                                 sg_mods_enabled, sg_mods_values) &
    result(err)
    integer(c_int), value, intent(in) :: matindex
    type(sg_eos_ary_t), intent(in)    :: eos
    real(kind=8), value, intent(in)   :: A, B, R1, R2, w, rho0, Cv
    integer(kind=c_int), dimension(:), target, intent(inout) :: sg_mods_enabled
    real(kind=8), dimension(:), target, intent(inout)        :: sg_mods_values
    err = init_sg_JWL(matindex-1, eos%ptr, A, B, R1, R2, w, rho0, Cv, &
                      c_loc(sg_mods_enabled), c_loc(sg_mods_values))
  end function init_sg_JWL_f
  
  integer function init_sg_DavisProducts_f(matindex, eos, a, b, k, n, vc, pc, &
                                           Cv, E0, sg_mods_enabled, &
                                           sg_mods_values) &
    result(err)
    integer(c_int), value, intent(in) :: matindex
    type(sg_eos_ary_t), intent(in)    :: eos
    real(kind=8), value, intent(in)   :: a, b, k, n, vc, pc, Cv, E0
    integer(kind=c_int), dimension(:), target, intent(inout) :: sg_mods_enabled
    real(kind=8), dimension(:), target, intent(inout)        :: sg_mods_values
    err = init_sg_DavisProducts(matindex-1, eos%ptr, a, b, k, n, vc, pc, Cv, &
                                E0, c_loc(sg_mods_enabled), &
                                c_loc(sg_mods_values))
  end function init_sg_DavisProducts_f

  integer function init_sg_DavisReactants_f(matindex, eos, rho0, e0, P0, T0, &
                                            A, B, C, G0, Z, alpha, Cv0, &
                                            sg_mods_enabled, sg_mods_values) &
    result(err)
    integer(c_int), value, intent(in) :: matindex
    type(sg_eos_ary_t), intent(in)    :: eos
    real(kind=8), value, intent(in)   :: rho0, e0, P0, T0, A, B, C, G0, Z
    real(kind=8), value, intent(in)   :: alpha, Cv0
    integer(kind=c_int), dimension(:), target, intent(inout) :: sg_mods_enabled
    real(kind=8), dimension(:), target, intent(inout)        :: sg_mods_values
    err = init_sg_DavisReactants(matindex-1, eos%ptr, rho0, e0, P0, T0, A, B, &
                                 C, G0, Z, alpha, Cv0, c_loc(sg_mods_enabled), &
                                 c_loc(sg_mods_values))
  end function init_sg_DavisReactants_f

  integer function init_sg_eospac_f(matindex, eos, id, eospac_opts_values, &
       sg_mods_enabled, sg_mods_values) &
    result(err)
    integer(c_int), value, intent(in) :: matindex, id
    type(sg_eos_ary_t), intent(in)    :: eos
    integer(kind=c_int), dimension(:), target, intent(inout) :: sg_mods_enabled
    real(kind=8), dimension(:), target, intent(inout)        :: sg_mods_values
    real(kind=8), dimension(:), target, intent(inout)        :: eospac_opts_values
    err = init_sg_eospac(matindex-1, eos%ptr, id, c_loc(eospac_opts_values), &
         c_loc(sg_mods_enabled), c_loc(sg_mods_values))
  end function init_sg_eospac_f

integer function get_sg_PressureFromDensityInternalEnergy_f(matindex, &
    eos, rhos, sies, pressures, len) &
    result(err)
    integer(c_int), intent(in) :: matindex, len
    real(kind=8), dimension(:,:,:), intent(in), target:: rhos, sies
    real(kind=8), dimension(:,:,:), intent(inout), target:: pressures
    type(sg_eos_ary_t), intent(in)    :: eos
    err = get_sg_PressureFromDensityInternalEnergy(matindex-1, &
           eos%ptr, c_loc(rhos(1,1,1)), c_loc(sies(1,1,1)), c_loc(pressures(1,1,1)), len)
  end function get_sg_PressureFromDensityInternalEnergy_f


integer function get_sg_MinInternalEnergyFromDensity_f(matindex, &
    eos, rhos, sies, len) &
    result(err)
    integer(c_int), intent(in) :: matindex, len
    real(kind=8), dimension(:,:,:), intent(in), target:: rhos
    real(kind=8), dimension(:,:,:), intent(inout), target:: sies
    type(sg_eos_ary_t), intent(in)    :: eos
    err = get_sg_MinInternalEnergyFromDensity(matindex-1, &
           eos%ptr, c_loc(rhos(1,1,1)), c_loc(sies(1,1,1)), len)
  end function get_sg_MinInternalEnergyFromDensity_f

  integer function get_sg_BulkModulusFromDensityInternalEnergy_f(matindex, &
    eos, rhos, sies, bmods, len) &
    result(err)
    integer(c_int), intent(in) :: matindex, len
    real(kind=8), dimension(:,:,:), intent(in), target:: rhos, sies
    real(kind=8), dimension(:,:,:), intent(inout), target:: bmods
    type(sg_eos_ary_t), intent(in)    :: eos
    err = get_sg_BulkModulusFromDensityInternalEnergy(matindex-1, &
       eos%ptr, c_loc(rhos(1,1,1)), c_loc(sies(1,1,1)), c_loc(bmods(1,1,1)), len)
  end function get_sg_BulkModulusFromDensityInternalEnergy_f

  integer function finalize_sg_eos_f(nmat, eos) &
    result(err)
    integer(c_int), value, intent(in) :: nmat
    type(sg_eos_ary_t), intent(in)    :: eos
    err = finalize_sg_eos(nmat, eos%ptr, 1)
  end function finalize_sg_eos_f
end module singularity_eos
