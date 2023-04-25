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

#include <map>
#include <ports-of-call/portability.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/get_sg_eos.hpp>
#include <singularity-eos/eos/get_sg_eos_lambdas.hpp>
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
    double *frac_bmod, double *frac_dpde, double *frac_cv) {
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
  const auto output = thermalqs::all_values - input;
  const bool p_is_inp{static_cast<bool>(input & thermalqs::pressure)};
  const bool r_is_inp{static_cast<bool>(input & thermalqs::density)};
  const bool t_is_inp{static_cast<bool>(input & thermalqs::temperature)};
  const bool s_is_inp{static_cast<bool>(input & thermalqs::specific_internal_energy)};
  constexpr const Real min_guess_density{1.e-20};
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

  // create helper lambdas to reduce code duplication
  const auto init_lambda = SG_GET_SG_EOS_INIT_LAMBDA_DECL;
  const auto final_lambda = SG_GET_SG_EOS_FINAL_LAMBDA_DECL;
  Kokkos::View<int, MemoryTraits<at_int>> res("PTE::num fails");
  Kokkos::View<int, MemoryTraits<at_int>> n_solves("PTE::num solves");
  const std::string perf_nums =
      "[" + std::to_string(nmat) + "," + std::to_string(ncell) + "]";
  auto input_int_enum = static_cast<input_condition>(input_int);
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
    singularity::get_sg_eos_rho_t(
        rt_name.c_str(), ncell, nmat, offsets_v, eos_offsets_v, eos_v, press_v, pmax_v,
        vol_v, spvol_v, sie_v, temp_v, bmod_v, dpde_v, cv_v, frac_mass_v, frac_vol_v,
        frac_ie_v, frac_bmod_v, frac_dpde_v, frac_cv_v, pte_idxs, pte_mats, press_pte,
        vfrac_pte, rho_pte, sie_pte, temp_pte, solver_scratch, tokens, small_loop,
        do_frac_bmod, do_frac_dpde, do_frac_cv);
    //    portableFor(
    //        rt_name.c_str(), 0, ncell, PORTABLE_LAMBDA(const int &iloop) {
    //          // cell offset
    //          const int i{offsets_v(iloop) - 1};
    //          // get "thread-id" like thing with optimization
    //          // for small loops
    //          const int32_t token{tokens.acquire()};
    //          const int32_t tid{small_loop ? iloop : token};
    //          double mass_sum{0.0};
    //          int npte{0};
    //          // initialize values for solver / lookup
    //          init_lambda(i, tid, mass_sum, npte, 1.0, 0.0, 0.0);
    //          // calculate pte condition (lookup for 1 mat cell)
    //          Real sie_tot_true{0.0};
    //          const int neq = npte;
    //          singularity::mix_impl::CacheAccessor cache(&solver_scratch(tid, 0) +
    //                                                     neq * (neq + 4) + 2 * npte);
    //          if (npte > 1) {
    //            // create solver lambda
    //            // eos accessor
    //            singularity::EOSAccessor_ eos_inx(eos_v, &pte_idxs(tid, 0));
    //            PTESolverFixedT<singularity::EOSAccessor_, Real *, Real **> method(
    //                npte, eos_inx, 1.0, temp_pte(tid, 0), &rho_pte(tid, 0),
    //                &vfrac_pte(tid, 0), &sie_pte(tid, 0), &temp_pte(tid, 0),
    //                &press_pte(tid, 0), cache, &solver_scratch(tid, 0));
    //            const bool res_{PTESolver(method)};
    //            // calculate total internal energy
    //            for (int mp = 0; mp < npte; ++mp) {
    //              const int m = pte_mats(tid, mp);
    //              sie_tot_true += sie_pte(tid, mp) * frac_mass_v(i, m);
    //            }
    //          } else {
    //            // pure cell (nmat = 1)
    //            // calculate sie from single eos
    //            sie_pte(tid, 0) = eos_v(pte_idxs(tid, 0))
    //                                  .InternalEnergyFromDensityTemperature(
    //                                      rho_pte(tid, 0), temp_pte(tid, 0), cache[0]);
    //            // set total sie to material 0 value
    //            sie_tot_true = sie_pte(tid, 0);
    //            // set pressure
    //            press_pte(tid, 0) = eos_v(pte_idxs(tid, 0))
    //                                    .PressureFromDensityTemperature(
    //                                        rho_pte(tid, 0), temp_pte(tid, 0),
    //                                        cache[0]);
    //          }
    //          // total sie is known
    //          sie_v(i) = sie_tot_true;
    //          // assign quantities
    //          final_lambda(i, tid, npte, mass_sum, 0.0, 0.0, 1.0, cache);
    //          // assign max pressure
    //          pmax_v(i) = press_v(i) > pmax_v(i) ? press_v(i) : pmax_v(i);
    //          // release the token used for scratch arrays
    //          tokens.release(token);
    //        });
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
    singularity::get_sg_eos_rho_p(
        rp_name.c_str(), ncell, nmat, offsets_v, eos_offsets_v, eos_v, press_v, pmax_v,
        vol_v, spvol_v, sie_v, temp_v, bmod_v, dpde_v, cv_v, frac_mass_v, frac_vol_v,
        frac_ie_v, frac_bmod_v, frac_dpde_v, frac_cv_v, pte_idxs, pte_mats, press_pte,
        vfrac_pte, rho_pte, sie_pte, temp_pte, solver_scratch, tokens, small_loop,
        do_frac_bmod, do_frac_dpde, do_frac_cv);
    //    portableFor(
    //        rp_name.c_str(), 0, ncell, PORTABLE_LAMBDA(const int &iloop) {
    //          // cell offset
    //          const int i{offsets_v(iloop) - 1};
    //          // get "thread-id" like thing with optimization
    //          // for small loops
    //          const int32_t token{tokens.acquire()};
    //          const int32_t tid{small_loop ? iloop : token};
    //          double mass_sum{0.0};
    //          int npte{0};
    //          // initialize values for solver / lookup
    //          init_lambda(i, tid, mass_sum, npte, 0.0, 0.0, 1.0);
    //          Real sie_tot_true{0.0};
    //          const int neq = npte + 1;
    //          singularity::mix_impl::CacheAccessor cache(&solver_scratch(tid, 0) +
    //                                                     neq * (neq + 4) + 2 * npte);
    //          if (npte > 1) {
    //            // create solver lambda
    //            // eos accessor
    //            singularity::EOSAccessor_ eos_inx(eos_v, &pte_idxs(tid, 0));
    //            PTESolverFixedP<singularity::EOSAccessor_, Real *, Real **> method(
    //                npte, eos_inx, 1.0, press_pte(tid, 0), &rho_pte(tid, 0),
    //                &vfrac_pte(tid, 0), &sie_pte(tid, 0), &temp_pte(tid, 0),
    //                &press_pte(tid, 0), cache[0], &solver_scratch(tid, 0));
    //            const bool res_{PTESolver(method)};
    //            // calculate total sie
    //            for (int mp = 0; mp < npte; ++mp) {
    //              const int m = pte_mats(tid, mp);
    //              sie_tot_true += sie_pte(tid, mp) * frac_mass_v(i, m);
    //            }
    //          } else {
    //            // pure cell (nmat = 1)
    //            // calculate sie from single eos
    //            auto p_from_t = [&](const Real &t_i) {
    //              return eos_v(pte_idxs(tid, 0))
    //                  .PressureFromDensityTemperature(rho_pte(tid, 0), t_i, cache[0]);
    //            };
    //            // calculate sie root bounds
    //            Real r_rho{}, r_temp{}, r_sie{}, r_press{}, r_cv{}, r_bmod{};
    //            Real r_dpde{}, r_dvdt{};
    //            eos_v(pte_idxs(tid, 0))
    //                .ValuesAtReferenceState(r_rho, r_temp, r_sie, r_press, r_cv, r_bmod,
    //                                        r_dpde, r_dvdt);
    //            Real temp_i;
    //            RootFinding1D::findRoot(p_from_t, press_v(i), r_temp, 1.e-10 * r_temp,
    //                                    1.e10 * r_temp, 1.e-12, 1.e-12, temp_i);
    //            sie_pte(tid, 0) = eos_v(pte_idxs(tid, 0))
    //                                  .InternalEnergyFromDensityTemperature(rho_pte(tid,
    //                                  0),
    //                                                                        temp_i,
    //                                                                        cache[0]);
    //            sie_tot_true = sie_pte(tid, 0);
    //            // set temperature
    //            temp_pte(tid, 0) = temp_i;
    //          }
    //          // sie calculate explicitly already
    //          sie_v(i) = sie_tot_true;
    //          // assign remaining outputs
    //          final_lambda(i, tid, npte, mass_sum, 1.0, 0.0, 0.0, cache);
    //          // assign max pressure
    //          pmax_v(i) = press_v(i) > pmax_v(i) ? press_v(i) : pmax_v(i);
    //          // release the token used for scratch arrays
    //          tokens.release(token);
    //        });
    break;
  }
  case input_condition::P_T_INPUT: {
    // P-T input
    const int pte_solver_scratch_size = nmat * MAX_NUM_LAMBDAS;
    solver_scratch = ScratchV<double>(VAWI("PTE::scratch solver"), scratch_size,
                                      pte_solver_scratch_size);
    const std::string pt_name = "PTE::solve (P,T) input" + perf_nums;
    singularity::get_sg_eos_p_t(
        pt_name.c_str(), ncell, nmat, offsets_v, eos_offsets_v, eos_v, press_v, pmax_v,
        vol_v, spvol_v, sie_v, temp_v, bmod_v, dpde_v, cv_v, frac_mass_v, frac_vol_v,
        frac_ie_v, frac_bmod_v, frac_dpde_v, frac_cv_v, pte_idxs, pte_mats, press_pte,
        vfrac_pte, rho_pte, sie_pte, temp_pte, solver_scratch, tokens, small_loop,
        do_frac_bmod, do_frac_dpde, do_frac_cv);
    //    portableFor(
    //        pt_name.c_str(), 0, ncell, PORTABLE_LAMBDA(const int &iloop) {
    //          // cell offset
    //          const int i{offsets_v(iloop) - 1};
    //          // get "thread-id" like thing with optimization
    //          // for small loops
    //          const int32_t token{tokens.acquire()};
    //          const int32_t tid{small_loop ? iloop : token};
    //          // caching mechanism
    //          singularity::mix_impl::CacheAccessor cache(&solver_scratch(tid, 0));
    //          double mass_sum{0.0};
    //          // normalize mass fractions
    //          // first find the mass sum
    //          // also set idxs as the decrement of the eos offsets
    //          // to take into account 1 based indexing in fortran
    //          for (int m = 0; m < nmat; ++m) {
    //            mass_sum += frac_mass_v(i, m);
    //            pte_idxs(tid, m) = eos_offsets_v(m) - 1;
    //            pte_mats(tid, m) = m;
    //            temp_pte(tid, m) = temp_v(i) * ev2k;
    //            press_pte(tid, m) = press_v(i);
    //          }
    //          for (int m = 0; m < nmat; ++m) {
    //            frac_mass_v(i, m) /= mass_sum;
    //          }
    //          // do r-e of pt for each mat
    //          singularity::EOSAccessor_ eos_inx(eos_v, &pte_idxs(tid, 0));
    //          Real vfrac_tot{0.0};
    //          Real sie_tot{0.0};
    //          for (int m = 0; m < nmat; ++m) {
    //            // obtain rho and sie from P-T
    //            eos_inx[m].DensityEnergyFromPressureTemperature(
    //                press_pte(tid, m), temp_pte(tid, m), cache[m], rho_pte(tid, m),
    //                sie_pte(tid, m));
    //            // assign volume fractions
    //            // this is a physical volume
    //            vfrac_pte(tid, m) = frac_mass_v(i, m) / rho_pte(tid, m) * mass_sum;
    //            vfrac_tot += vfrac_pte(tid, m);
    //            // add internal energy component
    //            sie_tot += sie_pte(tid, m) * frac_mass_v(i, m);
    //          }
    //          // assign volume, etc.
    //          // total sie is known
    //          sie_v(i) = sie_tot;
    //          vol_v(i) = vfrac_tot;
    //          spvol_v(i) = vol_v(i) / mass_sum;
    //          for (int m = 0; m < nmat; ++m) {
    //            vfrac_pte(tid, m) /= vfrac_tot;
    //          }
    //          // assign remaining outputs
    //          final_lambda(i, tid, nmat, mass_sum, 0.0, 0.0, 0.0, cache);
    //          // assign max pressure
    //          pmax_v(i) = press_v(i) > pmax_v(i) ? press_v(i) : pmax_v(i);
    //          // release the token used for scratch arrays
    //          tokens.release(token);
    //        });
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
    singularity::get_sg_eos_rho_p(
        re_name.c_str(), ncell, nmat, offsets_v, eos_offsets_v, eos_v, press_v, pmax_v,
        vol_v, spvol_v, sie_v, temp_v, bmod_v, dpde_v, cv_v, frac_mass_v, frac_vol_v,
        frac_ie_v, frac_bmod_v, frac_dpde_v, frac_cv_v, pte_idxs, pte_mats, press_pte,
        vfrac_pte, rho_pte, sie_pte, temp_pte, solver_scratch, tokens, small_loop,
        do_frac_bmod, do_frac_dpde, do_frac_cv);
    //    portableFor(
    //        re_name.c_str(), 0, ncell, PORTABLE_LAMBDA(const int &iloop) {
    //          // cell offset
    //          const int i{offsets_v(iloop) - 1};
    //          // get "thread-id" like thing with optimization
    //          // for small loops
    //          const int32_t token{tokens.acquire()};
    //          const int32_t tid{small_loop ? iloop : token};
    //          double mass_sum{0.0};
    //          int npte{0};
    //          // initialize values for solver / lookup
    //          init_lambda(i, tid, mass_sum, npte, 0.0, 1.0, 0.0);
    //          // get cache from offsets into scratch
    //          const int neq = npte + 1;
    //          singularity::mix_impl::CacheAccessor cache(&solver_scratch(tid, 0) +
    //                                                     neq * (neq + 4) + 2 * npte);
    //          if (npte > 1) {
    //            // create solver lambda
    //            // eos accessor
    //            singularity::EOSAccessor_ eos_inx(eos_v, &pte_idxs(tid, 0));
    //            // reset inputs
    //            PTESolverRhoT<singularity::EOSAccessor_, Real *, Real **> method(
    //                npte, eos_inx, 1.0, sie_v(i), &rho_pte(tid, 0), &vfrac_pte(tid, 0),
    //                &sie_pte(tid, 0), &temp_pte(tid, 0), &press_pte(tid, 0), cache,
    //                &solver_scratch(tid, 0));
    //            const bool res_{PTESolver(method)};
    //          } else {
    //            // pure cell (nmat = 1)
    //            // calculate sie from single eos
    //            Real dpdr_m, dtdr_m, dpde_m, dtde_m;
    //            eos_v(pte_idxs(tid, 0))
    //                .PTofRE(rho_pte(tid, 0), sie_pte(tid, 0), cache[0], press_pte(tid,
    //                0),
    //                        temp_pte(tid, 0), dpdr_m, dpde_m, dtdr_m, dtde_m);
    //          }
    //          // assign outputs
    //          final_lambda(i, tid, npte, mass_sum, 1.0, 0.0, 1.0, cache);
    //          // assign max pressure
    //          pmax_v(i) = press_v(i) > pmax_v(i) ? press_v(i) : pmax_v(i);
    //          // release the token used for scratch arrays
    //          tokens.release(token);
    //        });
    break;
  }
  }

  Kokkos::fence();
  // copy results back into local values
  // there is lots of room for performance optimization
  // in terms of when to copy and when not necessary
  // this is the return value (number of solve failures)
  deep_copy(ret, res);
  // copy pressure, this is not needed in all cases
  deep_copy(press_hv, press_v);
  // return max pressure, this may be needed
  deep_copy(pmax_hv, pmax_v);
  // I don't think the volume is necessary
  deep_copy(vol_hv, vol_v);
  // specific volume, copy-back not needed in all cases
  deep_copy(spvol_hv, spvol_v);
  // internal energy, copy-back not needed in all cases
  deep_copy(sie_hv, sie_v);
  // temperature, copy-back not needed in all cases
  deep_copy(temp_hv, temp_v);
  // bulk modulus, alwasy copy-back
  deep_copy(bmod_hv, bmod_v);
  // dpde, always copy-back
  deep_copy(dpde_hv, dpde_v);
  // specific heat, always copy-back
  deep_copy(cv_hv, cv_v);
  // volume fractions, always copy-back (maybe not for pure cells)
  deep_copy(frac_vol_hv, frac_vol_v);
  // component internal energies, always copy-back (maybe not for pure cells)
  deep_copy(frac_ie_hv, frac_ie_v);
  // optionally copy-back the component bmod, dpde, and cv
  if (do_frac_bmod) {
    deep_copy(frac_bmod_hv, frac_bmod_v);
  }
  if (do_frac_dpde) {
    deep_copy(frac_dpde_hv, frac_dpde_v);
  }
  if (do_frac_cv) {
    deep_copy(frac_cv_hv, frac_cv_v);
  }
#endif // PORTABILITY_STRATEGY_KOKKOS
  return ret;
}
