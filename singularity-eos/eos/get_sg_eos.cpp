//------------------------------------------------------------------------------
// © 2021-2024. Triad National Security, LLC. All rights reserved.  This
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

#include <map>
#include <ports-of-call/portability.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/get_sg_eos.hpp>
#include <singularity-eos/eos/get_sg_eos_functors.hpp>
#include <singularity-eos/eos/singularity_eos.hpp>

using namespace singularity;

// mapping from EAP integer to
static const std::map<const int, const unsigned long> EAPInputToBD = {
    {-3, thermalqs::temperature | thermalqs::density},
    {-2, thermalqs::density | thermalqs::pressure},
    {-1, thermalqs::pressure | thermalqs::temperature},
    {0, thermalqs::specific_internal_energy | thermalqs::density},
    {1, thermalqs::specific_internal_energy | thermalqs::density},
};

// EAP centric arguments and function signature
int get_sg_eos( // sizing information
    int nmat, int ncell, int cell_dim,
    // Input parameters
    int input_int,
    // eos index offsets
    int *eos_offsets,
    // equation of state array
    EOS *eos,
    // index offsets
    int *offsets,
    // per cell quantities
    double *press, double *pmax, double *vol, double *spvol, double *sie, double *temp,
    double *bmod, double *dpde, double *cv,
    // per material quantities
    double *frac_mass, double *frac_vol, double *frac_ie,
    // optional per material quantities
    double *frac_bmod, double *frac_dpde, double *frac_cv,
    // Mass fraction cutoff for PTE
    double mass_frac_cutoff) {
  // printBacktrace();
  // kernel return value will be the number of failures
  int ret{0};
  // handle optionals
  bool do_frac_bmod{true};
  if (frac_bmod == NULL || frac_bmod == nullptr) do_frac_bmod = false;
  bool do_frac_dpde{true};
  if (frac_dpde == NULL || frac_dpde == nullptr) do_frac_dpde = false;
  bool do_frac_cv{true};
  if (frac_dpde == NULL || frac_dpde == nullptr) do_frac_cv = false;
  // get inputs
  enum class input_condition {
    RHO_T_INPUT = -3,
    RHO_P_INPUT = -2,
    P_T_INPUT = -1,
    NORM_RHO_E_INPUT = 0,
    RHO_E_INPUT = 1
  };
  const auto input{EAPInputToBD.at(input_int)};
  const bool p_is_inp{static_cast<bool>(input & thermalqs::pressure)};
  const bool r_is_inp{static_cast<bool>(input & thermalqs::density)};
  const bool t_is_inp{static_cast<bool>(input & thermalqs::temperature)};
  const bool s_is_inp{static_cast<bool>(input & thermalqs::specific_internal_energy)};
#ifdef PORTABILITY_STRATEGY_KOKKOS
  // convert pointers to host side views
  Kokkos::View<int *, Llft, HS, Unmgd> eos_offsets_hv(eos_offsets, nmat);
  Kokkos::View<int *, Llft, HS, Unmgd> offsets_hv(offsets, ncell);
  host_v press_hv(press, cell_dim);
  host_v pmax_hv(pmax, cell_dim);
  host_v vol_hv(vol, cell_dim);
  host_v spvol_hv(spvol, cell_dim);
  host_v sie_hv(sie, cell_dim);
  host_v temp_hv(temp, cell_dim);
  host_v bmod_hv(bmod, cell_dim);
  host_v dpde_hv(dpde, cell_dim);
  host_v cv_hv(cv, cell_dim);
  host_frac_v frac_mass_hv(frac_mass, cell_dim, nmat);
  host_frac_v frac_vol_hv(frac_vol, cell_dim, nmat);
  host_frac_v frac_ie_hv(frac_ie, cell_dim, nmat);
  host_frac_v frac_bmod_hv, frac_dpde_hv, frac_cv_hv;
  if (do_frac_bmod) frac_bmod_hv = host_frac_v(frac_bmod, cell_dim, nmat);
  if (do_frac_dpde) frac_dpde_hv = host_frac_v(frac_dpde, cell_dim, nmat);
  if (do_frac_cv) frac_cv_hv = host_frac_v(frac_cv, cell_dim, nmat);

  // get device views if necessary
  indirection_v offsets_v{create_mirror_view_and_copy(DMS(), offsets_hv)};
  indirection_v eos_offsets_v{create_mirror_view_and_copy(DMS(), eos_offsets_hv)};
  dev_v press_v{create_mirror_view_and_copy(DMS(), press_hv)};
  dev_v pmax_v{create_mirror_view_and_copy(DMS(), pmax_hv)};
  dev_v spvol_v{create_mirror_view_and_copy(DMS(), spvol_hv)};
  dev_v vol_v{create_mirror_view_and_copy(DMS(), vol_hv)};
  dev_v sie_v{create_mirror_view_and_copy(DMS(), sie_hv)};
  dev_v temp_v{create_mirror_view_and_copy(DMS(), temp_hv)};
  dev_v bmod_v{create_mirror_view_and_copy(DMS(), bmod_hv)};
  dev_v dpde_v{create_mirror_view_and_copy(DMS(), dpde_hv)};
  dev_v cv_v{create_mirror_view_and_copy(DMS(), cv_hv)};
  dev_frac_v frac_mass_v{create_mirror_view_and_copy(DMS(), frac_mass_hv)};
  dev_frac_v frac_vol_v{create_mirror_view_and_copy(DMS(), frac_vol_hv)};
  dev_frac_v frac_ie_v{create_mirror_view_and_copy(DMS(), frac_ie_hv)};
  dev_frac_v frac_bmod_v, frac_dpde_v, frac_cv_v;
  if (do_frac_bmod) frac_bmod_v = create_mirror_view_and_copy(DMS(), frac_bmod_hv);
  if (do_frac_dpde) frac_dpde_v = create_mirror_view_and_copy(DMS(), frac_dpde_hv);
  if (do_frac_cv) frac_cv_v = create_mirror_view_and_copy(DMS(), frac_cv_hv);
  // array of eos's
  const auto eos_nmat{*std::max_element(eos_offsets, eos_offsets + nmat)};
  Kokkos::View<EOS *, Llft, HS, Unmgd> eos_hv(eos, eos_nmat);
  Kokkos::View<EOS *, Llft> eos_v{create_mirror_view_and_copy(DMS(), eos_hv)};
  // TODO: make this a scatter view
  constexpr auto at_int_full{0 | Kokkos::Atomic};
#ifdef KOKKOS_ENABLE_SERIAL
  constexpr auto at_int{std::is_same<DES, Kokkos::Serial>::value ? 0 : at_int_full};
#else  // KOKKOS_ENABLE_SERIAL
  // assume atomics are required for correctness if serial not even enabled
  constexpr auto at_int{at_int_full};
#endif // KOKKOS_ENABLE_SERIAL
  // set up scratch arrays
  constexpr auto KGlobal = Kokkos::Experimental::UniqueTokenScope::Global;
  Kokkos::Experimental::UniqueToken<DES, KGlobal> tokens{};
  using VAWI = Kokkos::ViewAllocateWithoutInitializing;

  const bool small_loop{tokens.size() > ncell};
  const decltype(tokens)::size_type scratch_size{std::min(tokens.size(), ncell)};
  ScratchV<int> pte_mats(VAWI("PTE::scratch mats"), scratch_size, nmat);
  ScratchV<int> pte_idxs(VAWI("PTE::scratch idxs"), scratch_size, nmat);
  ScratchV<double> mass_pte(VAWI("PTE::scratch mass"), scratch_size, nmat);
  ScratchV<double> sie_pte(VAWI("PTE::scratch sie"), scratch_size, nmat);
  ScratchV<double> vfrac_pte(VAWI("PTE::scratch vfrac"), scratch_size, nmat);
  ScratchV<double> temp_pte(VAWI("PTE::scratch temp"), scratch_size, nmat);
  ScratchV<double> press_pte(VAWI("PTE::scratch press"), scratch_size, nmat);
  ScratchV<double> rho_pte(VAWI("PTE::scratch rho"), scratch_size, nmat);
  // declare solver scratch, allocate in the case statement
  int pte_solver_scratch_size{};
  ScratchV<double> solver_scratch;
  // declare init and final functors
  auto input_int_enum = static_cast<input_condition>(input_int);
  init_functor i_func;
  final_functor f_func(spvol_v, temp_v, press_v, sie_v, bmod_v, cv_v, dpde_v, pte_mats,
                       press_pte, vfrac_pte, temp_pte, sie_pte, frac_mass_v, frac_ie_v,
                       frac_vol_v, vol_v, eos_v, pte_idxs, rho_pte, frac_bmod_v,
                       frac_cv_v, frac_dpde_v, nmat, do_frac_bmod, do_frac_cv,
                       do_frac_dpde);
  // only initialize init functor when needed
  if (input_int_enum != input_condition::P_T_INPUT) {
    i_func = init_functor(frac_mass_v, pte_idxs, eos_offsets_v, frac_vol_v, frac_ie_v,
                          pte_mats, vfrac_pte, sie_pte, temp_pte, press_pte, rho_pte,
                          spvol_v, temp_v, press_v, sie_v, nmat, mass_frac_cutoff);
  }

  // create helper lambdas to reduce code duplication
  Kokkos::View<int, MemoryTraits<at_int>> res("PTE::num fails");
  Kokkos::View<int, MemoryTraits<at_int>> n_solves("PTE::num solves");
  const std::string perf_nums =
      "[" + std::to_string(nmat) + "," + std::to_string(ncell) + "]";
  switch (input_int_enum) {
  case input_condition::RHO_T_INPUT: {
    // T-rho input
    // set frac_vol = 1/nmat
    // set rho_i = nmat / spvol * frac_mass_i
    // iterate PTE solver to obtain internal energies
    // that results in the input T
    pte_solver_scratch_size = PTESolverFixedTRequiredScratch(nmat);
    solver_scratch = ScratchV<double>(VAWI("PTE::scratch solver"), scratch_size,
                                      pte_solver_scratch_size);
    const std::string rt_name = "PTE::solve (rho,T) input" + perf_nums;
    singularity::get_sg_eos_rho_t(rt_name.c_str(), ncell, offsets_v, eos_v, press_v,
                                  pmax_v, sie_v, frac_mass_v, pte_idxs, pte_mats,
                                  press_pte, vfrac_pte, rho_pte, sie_pte, temp_pte,
                                  solver_scratch, tokens, small_loop, i_func, f_func);
    break;
  }
  case input_condition::RHO_P_INPUT: {
    // rho-P input
    // set frac_vol = 1/nmat
    // set rho_i = nmat / spvol * frac_mass_i
    // iterate PTE solver to obtain internal energies
    // that results in the input P
    pte_solver_scratch_size = PTESolverRhoTRequiredScratch(nmat);
    solver_scratch = ScratchV<double>(VAWI("PTE::scratch solver"), scratch_size,
                                      pte_solver_scratch_size);
    const std::string rp_name = "PTE::solve (rho,P) input" + perf_nums;
    singularity::get_sg_eos_rho_p(rp_name.c_str(), ncell, offsets_v, eos_v, press_v,
                                  pmax_v, sie_v, frac_mass_v, pte_idxs, pte_mats,
                                  press_pte, vfrac_pte, rho_pte, sie_pte, temp_pte,
                                  solver_scratch, tokens, small_loop, i_func, f_func);
    break;
  }
  case input_condition::P_T_INPUT: {
    // P-T input
    const int pte_solver_scratch_size = nmat * MAX_NUM_LAMBDAS;
    solver_scratch = ScratchV<double>(VAWI("PTE::scratch solver"), scratch_size,
                                      pte_solver_scratch_size);
    const std::string pt_name = "PTE::solve (P,T) input" + perf_nums;
    singularity::get_sg_eos_p_t(pt_name.c_str(), ncell, nmat, offsets_v, eos_offsets_v,
                                eos_v, press_v, pmax_v, vol_v, spvol_v, sie_v, temp_v,
                                frac_mass_v, pte_idxs, pte_mats, press_pte, vfrac_pte,
                                rho_pte, sie_pte, temp_pte, solver_scratch, tokens,
                                small_loop, f_func);
    break;
  }
  case input_condition::NORM_RHO_E_INPUT:
    // rho-sie input
    // no break so fallthrough to case 1
  case input_condition::RHO_E_INPUT: {
    // rho-sie input
    pte_solver_scratch_size = PTESolverRhoTRequiredScratch(nmat);
    solver_scratch = ScratchV<double>(VAWI("PTE::scratch solver"), scratch_size,
                                      pte_solver_scratch_size);
    const std::string re_name = "PTE::solve (rho,e) input" + perf_nums;
    singularity::get_sg_eos_rho_e(re_name.c_str(), ncell, offsets_v, eos_v, press_v,
                                  pmax_v, sie_v, pte_idxs, press_pte, vfrac_pte, rho_pte,
                                  sie_pte, temp_pte, solver_scratch, tokens, small_loop,
                                  i_func, f_func);
    break;
  }
  }

  Kokkos::fence();
  // copy results back into local values
  // there is lots of room for performance optimization
  // in terms of when to copy and when not necessary
  // this is the return value (number of solve failures)
  // deep_copy(ret, res);
  // copy pressure, this is not needed in all cases
  if (!p_is_inp) {
    deep_copy(DES(), press_hv, press_v);
  }
  // return max pressure, this may be needed
  deep_copy(DES(), pmax_hv, pmax_v);
  // I don't think the volume is necessary
  deep_copy(DES(), vol_hv, vol_v);
  // specific volume, copy-back not needed in all cases
  if (!r_is_inp) {
    deep_copy(DES(), spvol_hv, spvol_v);
  }
  // internal energy, copy-back not needed in all cases
  if (!s_is_inp) {
    deep_copy(DES(), sie_hv, sie_v);
  }
  // temperature, copy-back not needed in all cases
  if (!t_is_inp) {
    deep_copy(DES(), temp_hv, temp_v);
  }
  // bulk modulus, alwasy copy-back
  deep_copy(DES(), bmod_hv, bmod_v);
  // dpde, always copy-back
  deep_copy(DES(), dpde_hv, dpde_v);
  // specific heat, always copy-back
  deep_copy(DES(), cv_hv, cv_v);
  // volume fractions, always copy-back (maybe not for pure cells)
  deep_copy(DES(), frac_vol_hv, frac_vol_v);
  // component internal energies, always copy-back (maybe not for pure cells)
  deep_copy(DES(), frac_ie_hv, frac_ie_v);
  // optionally copy-back the component bmod, dpde, and cv
  if (do_frac_bmod) {
    deep_copy(DES(), frac_bmod_hv, frac_bmod_v);
  }
  if (do_frac_dpde) {
    deep_copy(DES(), frac_dpde_hv, frac_dpde_v);
  }
  if (do_frac_cv) {
    deep_copy(DES(), frac_cv_hv, frac_cv_v);
  }
#endif // PORTABILITY_STRATEGY_KOKKOS
  return ret;
}
