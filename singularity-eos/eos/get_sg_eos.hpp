//------------------------------------------------------------------------------
// © 2021-2023. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_EOS_GET_SG_EOS_HPP_
#define _SINGULARITY_EOS_EOS_GET_SG_EOS_HPP_

#include <ports-of-call/portability.hpp>
#include <singularity-eos/eos/eos.hpp>

using namespace singularity;

#ifdef PORTABILITY_STRATEGY_KOKKOS
using Lrgt = Kokkos::LayoutRight;
using Llft = Kokkos::LayoutLeft;
template <typename T>
using ScratchV = Kokkos::View<T **, Lrgt>;
#endif // PORTABILITY_STRATEGY_KOKKOS

#ifdef PORTABILITY_STRATEGY_KOKKOS
namespace singularity {
struct EOSAccessor_ {
  PORTABLE_INLINE_FUNCTION
  EOSAccessor_(const Kokkos::View<EOS *, Llft> &eos_v, int *mats)
      : eos_v_(eos_v), mats_(mats) {}
  PORTABLE_INLINE_FUNCTION
  EOS &operator[](const int m) const { return eos_v_(mats_[m]); }
  Kokkos::View<EOS *, Llft> eos_v_;
  int *mats_;
};
} // namespace singularity
#endif // PORTABILITY_STRATEGY_KOKKOS

#ifdef PORTABILITY_STRATEGY_KOKKOS
// convenience aliases
using Unmgd = Kokkos::MemoryUnmanaged;
using HS = Kokkos::HostSpace;
using Kokkos::create_mirror_view_and_copy;
using host_v = Kokkos::View<double *, Llft, HS, Unmgd>;
using dev_v = Kokkos::View<double *, Llft>;
using host_frac_v = Kokkos::View<double **, Llft, HS, Unmgd>;
using dev_frac_v = Kokkos::View<double **, Llft>;
using DES = Kokkos::DefaultExecutionSpace;
using DMS = DES::memory_space;
using Kokkos::MemoryTraits;
constexpr const unsigned int ra{0 | Kokkos::RandomAccess};
// constexpr const unsigned int ra_u{Kokkos::Unmanaged | Kokkos::RandomAccess};
using indirection_v = Kokkos::View<const int *, Llft, MemoryTraits<ra>>;
using VAWI = Kokkos::ViewAllocateWithoutInitializing;
using Kokkos::deep_copy;
constexpr auto KGlobal = Kokkos::Experimental::UniqueTokenScope::Global;
static constexpr const double ev2k = 1.160451930280894026e4;
#endif // PORTABILITY_STRATEGY_KOKKOS

namespace singularity {
// rho t input
void get_sg_eos_rho_t(const char* name, int ncell, int nmat, 
                      indirection_v& offsets_v,
		      indirection_v& eos_offsets_v,
		      Kokkos::View<EOS *, Llft>& eos_v,
                      dev_v& press_v,
                      dev_v& pmax_v,
                      dev_v& vol_v,
                      dev_v& spvol_v,
                      dev_v& sie_v,
                      dev_v& temp_v,
                      dev_v& bmod_v,
                      dev_v& dpde_v,
                      dev_v& cv_v,
                      dev_frac_v& frac_mass_v,
                      dev_frac_v& frac_vol_v,
                      dev_frac_v& frac_ie_v,
                      dev_frac_v& frac_bmod_v,
                      dev_frac_v& frac_dpde_v,
                      dev_frac_v& frac_cv_v,
                      ScratchV<int>& pte_idxs,
                      ScratchV<int>& pte_mats,
                      ScratchV<double>& press_pte,
                      ScratchV<double>& vfrac_pte,
                      ScratchV<double>& rho_pte,
                      ScratchV<double>& sie_pte,
                      ScratchV<double>& temp_pte,
                      ScratchV<double>& solver_scratch,
                      Kokkos::Experimental::UniqueToken<DES, KGlobal>& tokens,
		      bool small_loop,
                      bool do_frac_bmod, bool do_frac_dpde, bool do_frac_cv);
// rho P input
void get_sg_eos_rho_p(const char* name, int ncell, int nmat, 
                      indirection_v& offsets_v,
		      indirection_v& eos_offsets_v,
		      Kokkos::View<EOS *, Llft>& eos_v,
                      dev_v& press_v,
                      dev_v& pmax_v,
                      dev_v& vol_v,
                      dev_v& spvol_v,
                      dev_v& sie_v,
                      dev_v& temp_v,
                      dev_v& bmod_v,
                      dev_v& dpde_v,
                      dev_v& cv_v,
                      dev_frac_v& frac_mass_v,
                      dev_frac_v& frac_vol_v,
                      dev_frac_v& frac_ie_v,
                      dev_frac_v& frac_bmod_v,
                      dev_frac_v& frac_dpde_v,
                      dev_frac_v& frac_cv_v,
                      ScratchV<int>& pte_idxs,
                      ScratchV<int>& pte_mats,
                      ScratchV<double>& press_pte,
                      ScratchV<double>& vfrac_pte,
                      ScratchV<double>& rho_pte,
                      ScratchV<double>& sie_pte,
                      ScratchV<double>& temp_pte,
                      ScratchV<double>& solver_scratch,
                      Kokkos::Experimental::UniqueToken<DES, KGlobal>& tokens,
		      bool small_loop,
                      bool do_frac_bmod, bool do_frac_dpde, bool do_frac_cv);
// PT input
void get_sg_eos_p_t(const char* name, int ncell, int nmat, 
                    indirection_v& offsets_v,
		    indirection_v& eos_offsets_v,
		    Kokkos::View<EOS *, Llft>& eos_v,
                    dev_v& press_v,
                    dev_v& pmax_v,
                    dev_v& vol_v,
                    dev_v& spvol_v,
                    dev_v& sie_v,
                    dev_v& temp_v,
                    dev_v& bmod_v,
                    dev_v& dpde_v,
                    dev_v& cv_v,
                    dev_frac_v& frac_mass_v,
                    dev_frac_v& frac_vol_v,
                    dev_frac_v& frac_ie_v,
                    dev_frac_v& frac_bmod_v,
                    dev_frac_v& frac_dpde_v,
                    dev_frac_v& frac_cv_v,
                    ScratchV<int>& pte_idxs,
                    ScratchV<int>& pte_mats,
                    ScratchV<double>& press_pte,
                    ScratchV<double>& vfrac_pte,
                    ScratchV<double>& rho_pte,
                    ScratchV<double>& sie_pte,
                    ScratchV<double>& temp_pte,
                    ScratchV<double>& solver_scratch,
                    Kokkos::Experimental::UniqueToken<DES, KGlobal>& tokens,
		    bool small_loop,
                    bool do_frac_bmod, bool do_frac_dpde, bool do_frac_cv);
// rho e input
void get_sg_eos_rho_e(const char* name, int ncell, int nmat, 
                      indirection_v& offsets_v,
		      indirection_v& eos_offsets_v,
		      Kokkos::View<EOS *, Llft>& eos_v,
                      dev_v& press_v,
                      dev_v& pmax_v,
                      dev_v& vol_v,
                      dev_v& spvol_v,
                      dev_v& sie_v,
                      dev_v& temp_v,
                      dev_v& bmod_v,
                      dev_v& dpde_v,
                      dev_v& cv_v,
                      dev_frac_v& frac_mass_v,
                      dev_frac_v& frac_vol_v,
                      dev_frac_v& frac_ie_v,
                      dev_frac_v& frac_bmod_v,
                      dev_frac_v& frac_dpde_v,
                      dev_frac_v& frac_cv_v,
                      ScratchV<int>& pte_idxs,
                      ScratchV<int>& pte_mats,
                      ScratchV<double>& press_pte,
                      ScratchV<double>& vfrac_pte,
                      ScratchV<double>& rho_pte,
                      ScratchV<double>& sie_pte,
                      ScratchV<double>& temp_pte,
                      ScratchV<double>& solver_scratch,
                      Kokkos::Experimental::UniqueToken<DES, KGlobal>& tokens,
		      bool small_loop,
                      bool do_frac_bmod, bool do_frac_dpde, bool do_frac_cv);
} // namespace singularity

#endif // _SINGULARITY_EOS_EOS_GET_SG_EOS_HPP_
