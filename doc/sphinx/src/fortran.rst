.. _fortran:

Fortran Interface
================

Singularity-EOS provides a Fortran interface which can be enabled with the CMake
``SINGULARITY_USE_FORTRAN`` option. The interface provides Fortran bindings to the C++ EOS
types and gives access to both initialization functions and EOS evaluation functions.

The Fortran interface is defined in ``singularity-eos/eos/singularity_eos.f90``.


Modules
-------

The Fortran interface consists of two modules:

1. ``singularity_eos_types`` - Defines the data types used by the interface
2. ``singularity_eos`` - Contains the functions for initializing and using EOS models

Data Types
----------

.. code-block:: fortran

   type :: sg_eos_ary_t
     type(c_ptr) :: ptr = C_NULL_PTR
   end type sg_eos_ary_t

The ``sg_eos_ary_t`` type is a wrapper around a C pointer to the underlying EOS object.

Initialization Functions
-----------------------

The following functions are available for initializing EOS models:

.. code-block:: fortran

   integer function init_sg_eos_f(nmat, eos) result(err)
     integer(kind=c_int), intent(in) :: nmat
     type(sg_eos_ary_t), intent(in)  :: eos
   end function init_sg_eos_f

Initializes the EOS array with the specified number of materials.

.. code-block:: fortran

   integer function init_sg_IdealGas_f(matindex, eos, gm1, Cv, &
                                      sg_mods_enabled, sg_mods_values) &
     result(err)
     integer(c_int), value, intent(in) :: matindex
     type(sg_eos_ary_t), intent(in)    :: eos
     real(kind=8), value, intent(in)   :: gm1, Cv
     integer(kind=c_int), dimension(:), target, optional, intent(inout) :: sg_mods_enabled
     real(kind=8), dimension(:), target, optional, intent(inout)        :: sg_mods_values
   end function init_sg_IdealGas_f

Initializes an Ideal Gas EOS model at the specified material index.

- ``matindex``: Material index (1-based in Fortran)
- ``eos``: EOS array
- ``gm1``: Gamma minus 1 (γ-1)
- ``Cv``: Specific heat at constant volume
- ``sg_mods_enabled``: Optional array for enabling modifiers
- ``sg_mods_values``: Optional array for modifier values

.. code-block:: fortran

   integer function init_sg_Gruneisen_f(matindex, eos, C0, s1, s2, s3, G0, b, &
                                       rho0, T0, P0, Cv, sg_mods_enabled, &
                                       sg_mods_values) &
     result(err)
     integer(c_int), value, intent(in) :: matindex
     type(sg_eos_ary_t), intent(in)    :: eos
     real(kind=8), value, intent(in)   :: C0, s1, s2, s3, G0, b, rho0
     real(kind=8), value, intent(in)   :: T0, P0, Cv
     integer(kind=c_int), dimension(:), target, optional, intent(inout) :: sg_mods_enabled
     real(kind=8), dimension(:), target, optional, intent(inout)        :: sg_mods_values
   end function init_sg_Gruneisen_f

Initializes a Gruneisen EOS model at the specified material index.

- ``matindex``: Material index (1-based in Fortran)
- ``eos``: EOS array
- ``C0``: Bulk sound speed
- ``s1, s2, s3``: Coefficients for the shock Hugoniot
- ``G0``: Gruneisen parameter
- ``b``: Gruneisen parameter correction
- ``rho0``: Reference density
- ``T0``: Reference temperature
- ``P0``: Reference pressure
- ``Cv``: Specific heat at constant volume
- ``sg_mods_enabled``: Optional array for enabling modifiers
- ``sg_mods_values``: Optional array for modifier values

.. code-block:: fortran

   integer function init_sg_JWL_f(matindex, eos, A, B, R1, R2, w, rho0, Cv, &
                                 sg_mods_enabled, sg_mods_values) &
     result(err)
     integer(c_int), value, intent(in) :: matindex
     type(sg_eos_ary_t), intent(in)    :: eos
     real(kind=8), value, intent(in)   :: A, B, R1, R2, w, rho0, Cv
     integer(kind=c_int), dimension(:), target, optional, intent(inout) :: sg_mods_enabled
     real(kind=8), dimension(:), target, optional, intent(inout)        :: sg_mods_values
   end function init_sg_JWL_f

Initializes a JWL (Jones-Wilkins-Lee) EOS model at the specified material index.

- ``matindex``: Material index (1-based in Fortran)
- ``eos``: EOS array
- ``A, B``: JWL parameters
- ``R1, R2``: JWL parameters
- ``w``: JWL parameter
- ``rho0``: Reference density
- ``Cv``: Specific heat at constant volume
- ``sg_mods_enabled``: Optional array for enabling modifiers
- ``sg_mods_values``: Optional array for modifier values

.. code-block:: fortran

   integer function init_sg_DavisProducts_f(matindex, eos, a, b, k, n, vc, pc, &
                                           Cv, sg_mods_enabled, &
                                           sg_mods_values) &
     result(err)
     integer(c_int), value, intent(in) :: matindex
     type(sg_eos_ary_t), intent(in)    :: eos
     real(kind=8), value, intent(in)   :: a, b, k, n, vc, pc, Cv
     integer(kind=c_int), dimension(:), target, optional, intent(inout) :: sg_mods_enabled
     real(kind=8), dimension(:), target, optional, intent(inout)        :: sg_mods_values
   end function init_sg_DavisProducts_f

Initializes a Davis Products EOS model at the specified material index.

- ``matindex``: Material index (1-based in Fortran)
- ``eos``: EOS array
- ``a, b, k, n``: Davis Products parameters
- ``vc``: Critical specific volume
- ``pc``: Critical pressure
- ``Cv``: Specific heat at constant volume
- ``sg_mods_enabled``: Optional array for enabling modifiers
- ``sg_mods_values``: Optional array for modifier values

.. code-block:: fortran

   integer function init_sg_DavisReactants_f(matindex, eos, rho0, e0, P0, T0, &
                                            A, B, C, G0, Z, alpha, Cv0, &
                                            sg_mods_enabled, sg_mods_values) &
     result(err)
     integer(c_int), value, intent(in) :: matindex
     type(sg_eos_ary_t), intent(in)    :: eos
     real(kind=8), value, intent(in)   :: rho0, e0, P0, T0, A, B, C, G0, Z
     real(kind=8), value, intent(in)   :: alpha, Cv0
     integer(kind=c_int), dimension(:), target, optional, intent(inout) :: sg_mods_enabled
     real(kind=8), dimension(:), target, optional, intent(inout)        :: sg_mods_values
   end function init_sg_DavisReactants_f

Initializes a Davis Reactants EOS model at the specified material index.

- ``matindex``: Material index (1-based in Fortran)
- ``eos``: EOS array
- ``rho0``: Reference density
- ``e0``: Reference specific internal energy
- ``P0``: Reference pressure
- ``T0``: Reference temperature
- ``A, B, C``: Davis Reactants parameters
- ``G0``: Gruneisen parameter
- ``Z``: Davis Reactants parameter
- ``alpha``: Davis Reactants parameter
- ``Cv0``: Specific heat at constant volume
- ``sg_mods_enabled``: Optional array for enabling modifiers
- ``sg_mods_values``: Optional array for modifier values

.. code-block:: fortran

   integer function init_sg_NobleAbel_f(matindex, eos, gm1, Cv, &
                                      bb, qq, &
                                      sg_mods_enabled, sg_mods_values) &
     result(err)
     integer(c_int), value, intent(in) :: matindex
     type(sg_eos_ary_t), intent(in)    :: eos
     real(kind=8), value, intent(in)   :: gm1, Cv, bb, qq
     integer(kind=c_int), dimension(:), target, optional, intent(inout) :: sg_mods_enabled
     real(kind=8), dimension(:), target, optional, intent(inout)        :: sg_mods_values
   end function init_sg_NobleAbel_f

Initializes a Noble-Abel EOS model at the specified material index.

- ``matindex``: Material index (1-based in Fortran)
- ``eos``: EOS array
- ``gm1``: Gamma minus 1 (γ-1)
- ``Cv``: Specific heat at constant volume
- ``bb``: Covolume parameter
- ``qq``: Energy shift parameter
- ``sg_mods_enabled``: Optional array for enabling modifiers
- ``sg_mods_values``: Optional array for modifier values

.. code-block:: fortran

   integer function init_sg_StiffGas_f(matindex, eos, gm1, Cv, &
                                      Pinf, qq, &
                                      sg_mods_enabled, sg_mods_values) &
     result(err)
     integer(c_int), value, intent(in) :: matindex
     type(sg_eos_ary_t), intent(in)    :: eos
     real(kind=8), value, intent(in)   :: gm1, Cv, Pinf, qq
     integer(kind=c_int), dimension(:), target, optional, intent(inout) :: sg_mods_enabled
     real(kind=8), dimension(:), target, optional, intent(inout)        :: sg_mods_values
   end function init_sg_StiffGas_f

Initializes a Stiff Gas EOS model at the specified material index.

- ``matindex``: Material index (1-based in Fortran)
- ``eos``: EOS array
- ``gm1``: Gamma minus 1 (γ-1)
- ``Cv``: Specific heat at constant volume
- ``Pinf``: Pressure shift parameter
- ``qq``: Energy shift parameter
- ``sg_mods_enabled``: Optional array for enabling modifiers
- ``sg_mods_values``: Optional array for modifier values

.. code-block:: fortran

   integer function init_sg_SAP_Polynomial_f(matindex, eos, rho0, a0, a1, a2c, &
                                            a2e, a3, b0, b1, b2c, b2e, b3, &
                                            sg_mods_enabled, sg_mods_values) &
     result(err)
     integer(c_int), value, intent(in) :: matindex
     type(sg_eos_ary_t), intent(in)    :: eos
     real(kind=8), value, intent(in)   :: rho0, a0, a1, a2c, a2e, a3,&
              b0, b1, b2c, b2e, b3
     integer(kind=c_int), dimension(:), target, optional, intent(inout) :: sg_mods_enabled
     real(kind=8), dimension(:), target, optional, intent(inout)        :: sg_mods_values
   end function init_sg_SAP_Polynomial_f

Initializes a SAP (Separate Analytic Polynomials) Polynomial EOS model at the specified material index.

- ``matindex``: Material index (1-based in Fortran)
- ``eos``: EOS array
- ``rho0``: Reference density
- ``a0, a1, a2c, a2e, a3``: Compression polynomial coefficients
- ``b0, b1, b2c, b2e, b3``: Expansion polynomial coefficients
- ``sg_mods_enabled``: Optional array for enabling modifiers
- ``sg_mods_values``: Optional array for modifier values

Additional initialization functions are available when Singularity-EOS is built with specific options:

.. code-block:: fortran

   integer function init_sg_Helmholtz_f(matindex, eos, filename, rad, gas, coul, ion, ele, verbose, sg_mods_enabled, sg_mods_values) result(err)
     integer(c_int), value, intent(in) :: matindex
     type(sg_eos_ary_t), intent(in)    :: eos
     character(kind=c_char), dimension(*), intent(in) :: filename
     logical(c_bool), value, intent(in) :: rad, gas, coul, ion, ele, verbose
     integer(kind=c_int), dimension(:), target, optional, intent(inout) :: sg_mods_enabled
     real(kind=8), dimension(:), target, optional, intent(inout)        :: sg_mods_values
   end function init_sg_Helmholtz_f

Initializes a Helmholtz equation of state model at the specified material index. Available when built with ``SINGULARITY_USE_HELMHOLTZ`` and ``SINGULARITY_USE_SPINER_WITH_HDF5``.

- ``matindex``: Material index (1-based in Fortran)
- ``eos``: EOS array
- ``filename``: Path to the Helmholtz table file (typically "helm_table.dat")
- ``rad``: Boolean to enable/disable radiation term
- ``gas``: Boolean to enable/disable gas term
- ``coul``: Boolean to enable/disable Coulomb corrections
- ``ion``: Boolean to indicate if gas is ionized
- ``ele``: Boolean to indicate if electrons are degenerate
- ``verbose``: Boolean to enable verbose output
- ``sg_mods_enabled``: Optional array for enabling modifiers
- ``sg_mods_values``: Optional array for modifier values

.. code-block:: fortran

   integer function init_sg_SpinerDependsRhoT_f(matindex, eos, filename, id, split, sg_mods_enabled, sg_mods_values) result(err)
     integer(c_int), value, intent(in) :: matindex
     type(sg_eos_ary_t), intent(in)    :: eos
     character(kind=c_char), dimension(*), intent(in) :: filename
     integer(c_int), value, intent(in) :: id
     integer(c_int), value, intent(in) :: split
     integer(kind=c_int), dimension(:), target, optional, intent(inout) :: sg_mods_enabled
     real(kind=8), dimension(:), target, optional, intent(inout)        :: sg_mods_values
   end function init_sg_SpinerDependsRhoT_f

Initializes a tabulated EOS model that depends on density (ρ) and temperature (T) using the Spiner interpolation library. Available when built with ``SINGULARITY_USE_SPINER_WITH_HDF5``.

- ``matindex``: Material index (1-based in Fortran)
- ``eos``: EOS array
- ``filename``: Path to the HDF5 file containing the tabulated data
- ``id``: Material ID in the HDF5 file
- ``split``: Table split option (Total, ElectronOnly, or IonCold)
- ``sg_mods_enabled``: Optional array for enabling modifiers
- ``sg_mods_values``: Optional array for modifier values

.. code-block:: fortran

   integer function init_sg_SpinerDependsRhoSie_f(matindex, eos, filename, id, split, sg_mods_enabled, sg_mods_values) result(err)
     integer(c_int), value, intent(in) :: matindex
     type(sg_eos_ary_t), intent(in)    :: eos
     character(kind=c_char), dimension(*), intent(in) :: filename
     integer(c_int), value, intent(in) :: id
     integer(c_int), value, intent(in) :: split
     integer(kind=c_int), dimension(:), target, optional, intent(inout) :: sg_mods_enabled
     real(kind=8), dimension(:), target, optional, intent(inout)        :: sg_mods_values
   end function init_sg_SpinerDependsRhoSie_f

Initializes a tabulated EOS model that depends on density (ρ) and specific internal energy (sie) using the Spiner interpolation library. Available when built with ``SINGULARITY_USE_SPINER_WITH_HDF5``.

- ``matindex``: Material index (1-based in Fortran)
- ``eos``: EOS array
- ``filename``: Path to the HDF5 file containing the tabulated data
- ``id``: Material ID in the HDF5 file
- ``split``: Table split option (Total, ElectronOnly, or IonCold)
- ``sg_mods_enabled``: Optional array for enabling modifiers
- ``sg_mods_values``: Optional array for modifier values

.. code-block:: fortran

   integer function init_sg_eospac_f(matindex, eos, id, split, eospac_opts_values, sg_mods_enabled, sg_mods_values) result(err)
     integer(c_int), value, intent(in) :: matindex
     type(sg_eos_ary_t), intent(in)    :: eos
     integer(c_int), value, intent(in) :: id
     integer(c_int), value, intent(in) :: split
     real(kind=8), dimension(:), target, intent(in) :: eospac_opts_values
     integer(kind=c_int), dimension(:), target, optional, intent(inout) :: sg_mods_enabled
     real(kind=8), dimension(:), target, optional, intent(inout)        :: sg_mods_values
   end function init_sg_eospac_f

Initializes an EOSPAC equation of state model at the specified material index. Available when built with ``SINGULARITY_USE_EOSPAC``.

- ``matindex``: Material index (1-based in Fortran)
- ``eos``: EOS array
- ``id``: SESAME material ID
- ``split``: Table split option (Total, ElectronOnly, or IonCold)
- ``eospac_opts_values``: Array of EOSPAC-specific options:

  - ``eospac_opts_values(1)``: invert_at_setup (0/1)
  - ``eospac_opts_values(2)``: insert_data value
  - ``eospac_opts_values(3)``: monotonicity option
  - ``eospac_opts_values(4)``: apply_smoothing (0/1)
  - ``eospac_opts_values(5)``: apply_splitting option
  - ``eospac_opts_values(6)``: linear_interp (0/1)
- ``sg_mods_enabled``: Optional array for enabling modifiers
- ``sg_mods_values``: Optional array for modifier values


EOS Evaluation Functions
-----------------------

The following functions are available for evaluating EOS models:

.. code-block:: fortran

   integer function get_sg_EntropyFromDensityInternalEnergy_f(matindex, &
     eos, rhos, sies, entropies, len, stride, lambda_data) &
     result(err)
     integer(c_int), intent(in) :: matindex, len
     real(kind=8), dimension(:,:,:), intent(in), target:: rhos, sies
     real(kind=8), dimension(:,:,:), intent(inout), target:: entropies
     type(sg_eos_ary_t), intent(in)    :: eos
     integer(c_int), intent(in), optional :: stride
     real(kind=8), dimension(:,:,:,:), intent(inout), target, optional::lambda_data
   end function get_sg_EntropyFromDensityInternalEnergy_f

Calculates entropy from density and specific internal energy.

- ``matindex``: Material index (1-based in Fortran)
- ``eos``: EOS array
- ``rhos``: Input density array
- ``sies``: Input specific internal energy array
- ``entropies``: Output entropy array
- ``len``: Number of elements to process
- ``stride``: Optional stride for processing elements
- ``lambda_data``: Optional array for additional parameters

.. code-block:: fortran

   integer function get_sg_PressureFromDensityInternalEnergy_f(matindex, &
     eos, rhos, sies, pressures, len, stride, lambda_data) &
     result(err)
     integer(c_int), intent(in) :: matindex, len
     real(kind=8), dimension(:,:,:), intent(in), target:: rhos, sies
     real(kind=8), dimension(:,:,:), intent(inout), target:: pressures
     type(sg_eos_ary_t), intent(in)    :: eos
     integer(c_int), intent(in), optional :: stride
     real(kind=8), dimension(:,:,:,:), intent(inout), target, optional::lambda_data
   end function get_sg_PressureFromDensityInternalEnergy_f

Calculates pressure from density and specific internal energy.

- ``matindex``: Material index (1-based in Fortran)
- ``eos``: EOS array
- ``rhos``: Input density array
- ``sies``: Input specific internal energy array
- ``pressures``: Output pressure array
- ``len``: Number of elements to process
- ``stride``: Optional stride for processing elements
- ``lambda_data``: Optional array for additional parameters

.. code-block:: fortran

   integer function get_sg_MinInternalEnergyFromDensity_f(matindex, &
     eos, rhos, sies, len) &
     result(err)
     integer(c_int), intent(in) :: matindex, len
     real(kind=8), dimension(:,:,:), intent(in), target:: rhos
     real(kind=8), dimension(:,:,:), intent(inout), target:: sies
     type(sg_eos_ary_t), intent(in)    :: eos
   end function get_sg_MinInternalEnergyFromDensity_f

Calculates minimum internal energy from density.

- ``matindex``: Material index (1-based in Fortran)
- ``eos``: EOS array
- ``rhos``: Input density array
- ``sies``: Output minimum specific internal energy array
- ``len``: Number of elements to process

.. code-block:: fortran

   integer function get_sg_BulkModulusFromDensityInternalEnergy_f(matindex, &
     eos, rhos, sies, bmods, len, stride, lambda_data) &
     result(err)
     integer(c_int), intent(in) :: matindex, len
     real(kind=8), dimension(:,:,:), intent(in), target:: rhos, sies
     real(kind=8), dimension(:,:,:), intent(inout), target:: bmods
     type(sg_eos_ary_t), intent(in)    :: eos
     integer(c_int), intent(in), optional :: stride
     real(kind=8), dimension(:,:,:,:), intent(inout), target, optional::lambda_data
   end function get_sg_BulkModulusFromDensityInternalEnergy_f

Calculates bulk modulus from density and specific internal energy.

- ``matindex``: Material index (1-based in Fortran)
- ``eos``: EOS array
- ``rhos``: Input density array
- ``sies``: Input specific internal energy array
- ``bmods``: Output bulk modulus array
- ``len``: Number of elements to process
- ``stride``: Optional stride for processing elements
- ``lambda_data``: Optional array for additional parameters

.. code-block:: fortran

   integer function get_sg_eos_f(nmat, ncell, cell_dim,&
                                option,&
                                eos_offsets,&
                                eos,&
                                offsets,&
                                press, pmax, vol, spvol, sie, temp, bmod,&
                                dpde, cv,&
                                frac_mass, frac_vol, frac_sie,&
                                frac_bmod, frac_dpde, frac_cv,&
                                mass_frac_cutoff) &
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
     real(kind=8),                         optional, intent(in)    :: mass_frac_cutoff
   end function get_sg_eos_f

Main EOS evaluation function for mixed cells. This function is used to calculate various thermodynamic quantities for multiple materials in mixed cells.

- ``nmat``: Number of materials
- ``ncell``: Number of cells
- ``cell_dim``: Cell dimension
- ``option``: Option flag for calculation
- ``eos_offsets``: Array of EOS offsets
- ``eos``: EOS array
- ``offsets``: Array of offsets
- ``press``: Pressure array
- ``pmax``: Maximum pressure array
- ``vol``: Volume array
- ``spvol``: Specific volume array
- ``sie``: Specific internal energy array
- ``temp``: Temperature array
- ``bmod``: Bulk modulus array
- ``dpde``: Derivative of pressure with respect to energy array
- ``cv``: Specific heat array
- ``frac_mass``: Mass fraction array
- ``frac_vol``: Volume fraction array
- ``frac_sie``: Specific internal energy fraction array
- ``frac_bmod``: Optional bulk modulus fraction array
- ``frac_dpde``: Optional derivative of pressure with respect to energy fraction array
- ``frac_cv``: Optional specific heat fraction array
- ``mass_frac_cutoff``: Optional mass fraction cutoff

Finalization Function
--------------------

.. code-block:: fortran

   integer function finalize_sg_eos_f(nmat, eos) &
     result(err)
     integer(c_int), value, intent(in) :: nmat
     type(sg_eos_ary_t), intent(in)    :: eos
   end function finalize_sg_eos_f

Finalizes the EOS array, freeing any allocated resources.

- ``nmat``: Number of materials
- ``eos``: EOS array

Examples
--------

Basic Example
~~~~~~~~~~~~

Here's a simple example of using the Fortran interface to initialize and use an Ideal Gas EOS:

.. code-block:: fortran

   program example_sg_fortran_interface
   ! modules
   use singularity_eos_types
   use singularity_eos
   ! no implicit vars
   implicit none
   ! variable declaration
   integer                                   :: nmat, res, mat
   type(sg_eos_ary_t)                        :: eos
   real(kind=8), dimension(1,1,1)            :: rhos, sies, pressures

   ! set test parameters
   nmat = 1

   ! allocate and initialize eos's
   res = init_sg_eos_f(nmat, eos)
   mat = 1

   ! Initialize an ideal gas EOS
   res = init_sg_IdealGas_f(mat, eos, 0.4d0, 1.0d7)

   ! Set input values
   rhos(1,1,1) = 1.0d0  ! density in g/cm^3
   sies(1,1,1) = 1.0d6  ! specific internal energy in erg/g

   ! Calculate pressure
   res = get_sg_PressureFromDensityInternalEnergy_f(mat, eos, rhos, sies, pressures, 1)

   ! Print result
   print *, "Pressure = ", pressures(1,1,1), " dyne/cm^2"

   ! cleanup
   res = finalize_sg_eos_f(nmat, eos)

   end program example_sg_fortran_interface

Advanced Example with Tabulated EOS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example demonstrates using tabulated EOS models with a density-temperature grid:

.. code-block:: fortran

   program tabulated_eos_example
       use singularity_eos_types
       use singularity_eos
       implicit none
   
       integer :: nmat, res, mat, ncell, cell_dim, option, split
       type(sg_eos_ary_t) :: eos
       integer :: matid
       integer, allocatable :: eos_offsets(:), offsets(:)
       real(8), allocatable :: press(:), pmax(:), vol(:), spvol(:), sie(:), temp(:)
       real(8), allocatable :: bmod(:), dpde(:), cv(:)
       real(8), allocatable :: density(:)
       
       ! Grid parameters
       integer :: ndens, ntemp, idx
       real(8) :: dens_min, dens_max, temp_min, temp_max
       real(8), dimension(6) :: eospac_opts
       integer :: i, j
   
       ! Set up parameters
       nmat = 1
       ndens = 10
       ntemp = 10
       ncell = ndens * ntemp
       
       ! Define grid ranges
       dens_min = 0.1    ! g/cm^3
       dens_max = 10.0   ! g/cm^3
       temp_min = 0.1    ! eV (note: temperature in eV, not Kelvin)
       temp_max = 100.0  ! eV
       
       cell_dim = ncell
       option = -3  ! Option flag: -3 means density and temperature as inputs
       split = 1 ! Use ion subtable
   
       ! Allocate arrays
       allocate(eos_offsets(nmat), offsets(ncell), density(ncell))
       allocate(press(ncell), pmax(ncell), vol(ncell), spvol(ncell))
       allocate(sie(ncell), temp(ncell), bmod(ncell), dpde(ncell), cv(ncell))
   
       ! Initialize arrays
       eos_offsets = 1
       offsets = [(i, i=1,ncell)]
       press = 0.0d0
       pmax = 0.0d0
       vol = 0.0d0
       spvol = 0.0d0
       sie = 0.0d0
       bmod = 0.0d0
       dpde = 0.0d0
       cv = 0.0d0
       
       ! Set up density-temperature grid
       idx = 0
       do i = 1, ndens
           do j = 1, ntemp
               idx = idx + 1
               ! Logarithmic spacing
               density(idx) = dens_min * (dens_max/dens_min)**((i-1.0d0)/(ndens-1.0d0))
               spvol(idx) = 1.0d0 / density(idx)
               temp(idx) = temp_min * (temp_max/temp_min)**((j-1.0d0)/(ntemp-1.0d0))
           end do
       end do
   
       ! Initialize EOS
       res = init_sg_eos_f(nmat, eos)
       if (res /= 0) then
           print *, "Error initializing EOS:", res
           stop
       end if
   
       ! Set EOSPAC options (all zeros for default behavior)
       eospac_opts = 0.0d0
   
       ! Initialize material with EOSPAC
       mat = 1
       matid = 9999  ! Example material ID
       res = init_sg_eospac_f(mat, eos, matid, split, eospac_opts)
       if (res /= 0) then
           print *, "Error initializing EOSPAC:", res
           stop
       end if
   
       ! Calculate EOS properties
       ! Note: Temperature is in eV for get_sg_eos_f
       res = get_sg_eos_f(nmat, ncell, cell_dim, option, &
                         eos_offsets, eos, offsets, &
                         press, pmax, vol, spvol, sie, temp, bmod, dpde, cv)
       if (res /= 0) then
           print *, "Error in get_sg_eos_f:", res
           stop
       end if
   
       ! Print some results
       print *, "Results for first grid point:"
       print *, "Density = ", density(1), " g/cm^3"
       print *, "Temperature = ", temp(1), " eV"
       print *, "Pressure = ", press(1), " dyne/cm^2"
       print *, "SIE = ", sie(1), " erg/g"
   
       ! Clean up
       deallocate(eos_offsets, offsets, press, pmax, vol, spvol)
       deallocate(sie, temp, bmod, dpde, cv, density)
   
       res = finalize_sg_eos_f(nmat, eos)
       if (res /= 0) then
           print *, "Error finalizing EOS:", res
           stop
       end if
   
   end program tabulated_eos_example

.. note::
   
   When using the ``get_sg_eos_f`` function, temperature inputs should be provided in electron volts (eV), not Kelvin. This is different from some other functions in the interface that expect temperature in Kelvin.

