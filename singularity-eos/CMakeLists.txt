#------------------------------------------------------------------------------#
# © 2021-2024. Triad National Security, LLC. All rights reserved.  This
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
#------------------------------------------------------------------------------#

# Special sauce so generated file has proper include path
configure_file(eos/eos.hpp.in
  ${CMAKE_BINARY_DIR}/generated/singularity-eos/eos/eos.hpp
  @ONLY)

register_headers(
    # tell cmake:
    # (1) that if this file is changed, we gotta rebuild
    # (2) to copy this file into the install for posterity
    eos/eos.hpp.in

    # Normal files
    base/fast-math/logs.hpp
    base/indexable_types.hpp
    base/robust_utils.hpp
    base/root-finding-1d/root_finding.hpp
    base/serialization_utils.hpp
    base/spiner_table_utils.hpp
    base/variadic_utils.hpp
    base/math_utils.hpp
    base/constants.hpp
    base/eos_error.hpp
    base/error_utils.hpp
    base/sp5/singularity_eos_sp5.hpp
    eos/default_variant.hpp
    base/hermite.hpp
    eos/eos_variant.hpp
    eos/eos_stellar_collapse.hpp
    eos/eos_ideal.hpp
    eos/eos_models.hpp
    eos/eos_mgusup.hpp
    eos/eos_powermg.hpp
    eos/eos_spiner.hpp
    eos/eos_davis.hpp
    eos/eos_gruneisen.hpp
    eos/eos_vinet.hpp
    eos/eos_builder.hpp
    eos/eos_jwl.hpp
    eos/eos_helmholtz.hpp
    eos/eos_sap_polynomial.hpp
    eos/eos_type_lists.hpp
    eos/modifiers/relativistic_eos.hpp
    eos/modifiers/scaled_eos.hpp
    eos/modifiers/ramps_eos.hpp
    eos/modifiers/shifted_eos.hpp
    eos/modifiers/eos_unitsystem.hpp
    eos/modifiers/zsplit_eos.hpp
    eos/eos_base.hpp
    eos/eos_eospac.hpp
    eos/eos_noble_abel.hpp
    eos/eos_stiff.hpp
    eos/singularity_eos_init_utils.hpp
    eos/variant_utils.hpp
    eos/eos_carnahan_starling.hpp
    eos/eos_electrons.hpp
)

if (SINGULARITY_INSTANTIATE_CLASSES)
  register_srcs(
    eos/eos_carnahan_starling.cpp
    eos/eos_davis.cpp
    eos/eos_electrons.cpp
    eos/eos_eospac.cpp
    eos/eos_gruneisen.cpp
    eos/eos_helmholtz.cpp
    eos/eos_ideal.cpp
    eos/eos_jwl.cpp
    eos/eos_mgusup.cpp
    eos/eos_noble_abel.cpp
    eos/eos_powermg.cpp
    eos/eos_sap_polynomial.cpp
    eos/eos_spiner.cpp
    eos/eos_stellar_collapse.cpp
    eos/eos_stiff.cpp
    eos/eos_vinet.cpp
    eos/modifiers/shifted_eos.cpp
    eos/modifiers/scaled_eos.cpp
    )
endif()

if (SINGULARITY_BUILD_CLOSURE)
  register_headers(
    closure/mixed_cell_models.hpp
    closure/kinetic_phasetransition_models.hpp
    closure/kinetic_phasetransition_utils.hpp
    closure/kinetic_phasetransition_methods.hpp
  )
  if (SINGULARITY_BUILD_FORTRAN_BACKEND OR SINGULARITY_BUILD_TESTS)
    # while these are C++ files, they
    # are only needed for the fortran backend or unit testing
    register_srcs(eos/singularity_eos.cpp)
    register_srcs(eos/get_sg_eos_misc.cpp)
    register_headers(eos/singularity_eos.hpp)
  endif()
  if (SINGULARITY_BUILD_FORTRAN_BACKEND)
    register_srcs(eos/get_sg_eos.cpp)
    if (SINGULARITY_USE_KOKKOS)
      register_srcs(
           eos/get_sg_eos_p_t.cpp
           eos/get_sg_eos_rho_t.cpp
           eos/get_sg_eos_rho_p.cpp
           eos/get_sg_eos_rho_e.cpp)
    endif()
    register_headers(
         eos/get_sg_eos.hpp
         eos/get_sg_eos_functors.hpp)
  endif()
endif()

if (SINGULARITY_USE_FORTRAN)
  register_srcs(
    eos/singularity_eos.f90
  )
endif()
