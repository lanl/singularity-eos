//------------------------------------------------------------------------------
// © 2021-2025. Triad National Security, LLC. All rights reserved.  This
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

#include <singularity-eos/eos/eos_carnahan_starling.hpp>
#include <singularity-eos/eos/eos_davis.hpp>
#include <singularity-eos/eos/eos_electrons.hpp>
#include <singularity-eos/eos/eos_eospac.hpp>
#include <singularity-eos/eos/eos_gruneisen.hpp>
#include <singularity-eos/eos/eos_helmholtz.hpp>
#include <singularity-eos/eos/eos_ideal.hpp>
#include <singularity-eos/eos/eos_jwl.hpp>
#include <singularity-eos/eos/eos_mgusup.hpp>
#include <singularity-eos/eos/eos_noble_abel.hpp>
#include <singularity-eos/eos/eos_powermg.hpp>
#include <singularity-eos/eos/eos_sap_polynomial.hpp>
#include <singularity-eos/eos/eos_spiner.hpp>
#include <singularity-eos/eos/eos_stellar_collapse.hpp>
#include <singularity-eos/eos/eos_stiff.hpp>
#include <singularity-eos/eos/eos_vinet.hpp>

#include <singularity-eos/eos/modifiers/shifted_eos.hpp>

namespace singularity {

template class ShiftedEOS<IdealGas>;
template class ShiftedEOS<Gruneisen>;
template class ShiftedEOS<Vinet>;
template class ShiftedEOS<MGUsup>;
template class ShiftedEOS<PowerMG>;
template class ShiftedEOS<JWL>;
template class ShiftedEOS<DavisReactants>;
template class ShiftedEOS<DavisProducts>;

} // namespace singularity
