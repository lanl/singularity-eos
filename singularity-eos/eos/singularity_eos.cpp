//------------------------------------------------------------------------------
// © 2021-2022. Triad National Security, LLC. All rights reserved.  This
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

#include <algorithm>
#include <cassert>
#include <map>
#include <singularity-eos/eos/eos_builder.hpp>
#include <singularity-eos/eos/singularity_eos.hpp>

using namespace singularity;

template <typename T,
          typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
constexpr bool isfinite(const T &a) {
  return (a == a) && ((a == 0) || (a != 2 * a));
}

int init_sg_eos(const int nmat, EOS *&eos) {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  if (!Kokkos::is_initialized()) Kokkos::initialize();
#endif // PORTABILITY_STRATEGY_KOKKOS
  EOS *eos_p = new EOS[nmat];
  eos = eos_p;
  return 0;
}

#define SGAPPLYMOD(A)                                                                    \
  EOSBuilder::applyShiftAndScale(A, enabled[0] == 1, enabled[1] == 1, vals[0], vals[1])

constexpr const int def_en[2] = {0, 0};
constexpr const double def_v[2] = {0.0, 0.0};

int init_sg_IdealGas(const int matindex, EOS *eos, const double gm1, const double Cv,
                     int const *const enabled, double const *const vals) {
  assert(matindex >= 0);
  EOS eos_ = SGAPPLYMOD(IdealGas(gm1, Cv));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}

int init_sg_IdealGas(const int matindex, EOS *eos, const double gm1, const double Cv) {
  return init_sg_IdealGas(matindex, eos, gm1, Cv, def_en, def_v);
}

int init_sg_Gruneisen(const int matindex, EOS *eos, const double C0, const double s1,
                      const double s2, const double s3, const double G0, const double b,
                      const double rho0, const double T0, const double P0,
                      const double Cv, int const *const enabled,
                      double const *const vals) {
  assert(matindex >= 0);
  EOS eos_ = SGAPPLYMOD(Gruneisen(C0, s1, s2, s3, G0, b, rho0, T0, P0, Cv));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}

int init_sg_Gruneisen(const int matindex, EOS *eos, const double C0, const double s1,
                      const double s2, const double s3, const double G0, const double b,
                      const double rho0, const double T0, const double P0,
                      const double Cv) {
  return init_sg_Gruneisen(matindex, eos, C0, s1, s2, s3, G0, b, rho0, T0, P0, Cv, def_en,
                           def_v);
}

int init_sg_JWL(const int matindex, EOS *eos, const double A, const double B,
                const double R1, const double R2, const double w, const double rho0,
                const double Cv, int const *const enabled, double const *const vals) {
  assert(matindex >= 0);
  EOS eos_ = SGAPPLYMOD(JWL(A, B, R1, R2, w, rho0, Cv));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}

int init_sg_JWL(const int matindex, EOS *eos, const double A, const double B,
                const double R1, const double R2, const double w, const double rho0,
                const double Cv) {
  return init_sg_JWL(matindex, eos, A, B, R1, R2, w, rho0, Cv, def_en, def_v);
}

int init_sg_DavisProducts(const int matindex, EOS *eos, const double a, const double b,
                          const double k, const double n, const double vc,
                          const double pc, const double Cv, const double E0,
                          int const *const enabled, double const *const vals) {
  assert(matindex >= 0);
  EOS eos_ = SGAPPLYMOD(DavisProducts(a, b, k, n, vc, pc, Cv, E0));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}

int init_sg_DavisProducts(const int matindex, EOS *eos, const double a, const double b,
                          const double k, const double n, const double vc,
                          const double pc, const double Cv, const double E0) {
  return init_sg_DavisProducts(matindex, eos, a, b, k, n, vc, pc, Cv, E0, def_en, def_v);
}

int init_sg_DavisReactants(const int matindex, EOS *eos, const double rho0,
                           const double e0, const double P0, const double T0,
                           const double A, const double B, const double C,
                           const double G0, const double Z, const double alpha,
                           const double Cv0, int const *const enabled,
                           double const *const vals) {
  assert(matindex >= 0);
  EOS eos_ = SGAPPLYMOD(DavisReactants(rho0, e0, P0, T0, A, B, C, G0, Z, alpha, Cv0));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}

int init_sg_DavisReactants(const int matindex, EOS *eos, const double rho0,
                           const double e0, const double P0, const double T0,
                           const double A, const double B, const double C,
                           const double G0, const double Z, const double alpha,
                           const double Cv0) {
  return init_sg_DavisReactants(matindex, eos, rho0, e0, P0, T0, A, B, C, G0, Z, alpha,
                                Cv0, def_en, def_v);
}

#ifdef SPINER_USE_HDF
int init_sg_SpinerDependsRhoT(const int matindex, EOS *eos, const char *filename,
                              const int matid, int const *const enabled,
                              double const *const vals) {
  assert(matindex >= 0);
  EOS eos_ = SGAPPLYMOD(SpinerEOSDependsRhoT(std::string(filename), matid));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}

int init_sg_SpinerDependsRhoT(const int matindex, EOS *eos, const char *filename,
                              const int matid) {
  return init_sg_SpinerDependsRhoT(matindex, eos, filename, matid, def_en, def_v);
}

int init_sg_SpinerDependsRhoSie(const int matindex, EOS *eos, const char *filename,
                                const int matid, int const *const enabled,
                                double const *const vals) {
  assert(matindex >= 0);
  EOS eos_ = SGAPPLYMOD(SpinerEOSDependsRhoSie(std::string(filename), matid));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}
int init_sg_SpinerDependsRhoSie(const int matindex, EOS *eos, const char *filename,
                                const int matid) {
  return init_sg_SpinerDependsRhoSie(matindex, eos, filename, matid, def_en, def_v);
}
#endif

#ifdef SINGULARITY_USE_EOSPAC
int init_sg_eospac(const int matindex, EOS *eos, const int id, int const *const enabled,
                   double const *const vals) {
  assert(matindex >= 0);
  EOS eos_ = SGAPPLYMOD(EOSPAC(id));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}
int init_sg_eospac(const int matindex, EOS *eos, const int id) {
  return init_sg_eospac(matindex, eos, id, def_en, def_v);
}
#endif // SINGULARITY_USE_EOSPAC

#undef SGAPPLYMOD

#ifdef PORTABILITY_STRATEGY_KOKKOS
using Lrgt = Kokkos::LayoutRight;
template <typename T>
using ScratchV = Kokkos::View<T **, Lrgt>;
#endif // PORTABILITY_STRATEGY_KOKKOS

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
    double *frac_mass, double *frac_vol, double *frac_sie,
    // optional per material quantities
    double *frac_bmod, double *frac_dpde, double *frac_cv) {
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
  const auto input{EAPInputToBD.at(input_int)};
  const auto output = thermalqs::all_values - input;
  const bool p_is_inp{static_cast<bool>(input & thermalqs::pressure)};
  const bool r_is_inp{static_cast<bool>(input & thermalqs::density)};
  const bool t_is_inp{static_cast<bool>(input & thermalqs::temperature)};
  const bool s_is_inp{static_cast<bool>(input & thermalqs::specific_internal_energy)};
  constexpr const Real min_guess_density{1.e-20};
#ifdef PORTABILITY_STRATEGY_KOKKOS
  // convenience aliases
  using Llft = Kokkos::LayoutLeft;
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
  constexpr const unsigned int ra_u{Kokkos::Unmanaged | Kokkos::RandomAccess};
  using VAWI = Kokkos::ViewAllocateWithoutInitializing;
  using Kokkos::deep_copy;
  static constexpr const double ev2k = 1.160451930280894026e4;
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
  host_frac_v frac_sie_hv(frac_sie, cell_dim, nmat);
  host_frac_v frac_bmod_hv, frac_dpde_hv, frac_cv_hv;
  if (do_frac_bmod) frac_bmod_hv = host_frac_v(frac_bmod, cell_dim, nmat);
  if (do_frac_dpde) frac_dpde_hv = host_frac_v(frac_dpde, cell_dim, nmat);
  if (do_frac_cv) frac_cv_hv = host_frac_v(frac_cv, cell_dim, nmat);

  // get device views if necessary
  Kokkos::View<const int *, Llft, MemoryTraits<ra>> offsets_v{
      create_mirror_view_and_copy(DMS(), offsets_hv)};
  Kokkos::View<const int *, Llft, MemoryTraits<ra>> eos_offsets_v{
      create_mirror_view_and_copy(DMS(), eos_offsets_hv)};
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
  dev_frac_v frac_sie_v{create_mirror_view_and_copy(DMS(), frac_sie_hv)};
  dev_frac_v frac_bmod_v, frac_dpde_v, frac_cv_v;
  if (do_frac_bmod) frac_bmod_v = create_mirror_view_and_copy(DMS(), frac_bmod_hv);
  if (do_frac_dpde) frac_dpde_v = create_mirror_view_and_copy(DMS(), frac_dpde_hv);
  if (do_frac_cv) frac_cv_v = create_mirror_view_and_copy(DMS(), frac_cv_hv);
  // array of eos's
  const auto eos_nmat{*std::max_element(eos_offsets, eos_offsets + nmat)};
  Kokkos::View<EOS *, Llft, HS, Unmgd> eos_hv(eos, eos_nmat);
  Kokkos::View<EOS *, Llft> eos_v{create_mirror_view_and_copy(DMS(), eos_hv)};
  double **lambda_map = (double **)PORTABLE_MALLOC(sizeof(double *) * nmat);
  portableFor(
      "filling lambda ptrs", 0, nmat, PORTABLE_LAMBDA(const int &i) {
        // const int lambda_size = eos_v(i).nlambda();
        lambda_map[i] = nullptr; //(double*)malloc(sizeof(double)*lambda_size);
      });
  // TODO: make this a scatter view
  constexpr auto at_int_full{0 | Kokkos::Atomic};
#ifdef KOKKOS_ENABLE_SERIAL
  constexpr auto at_int{std::is_same<DES, Kokkos::Serial>::value ? 0 : at_int_full};
#else  // KOKKOS_ENABLE_SERIAL
  // assume atomics are required for correctness if serial not even enabled
  constexpr auto at_int{at_int_full};
#endif // KOKKOS_ENABLE_SERIAL
  Kokkos::View<int *> nmat_cell_init(VAWI("nmat_cell"), ncell);
  // minimum mass fraction necessary to participate
  constexpr const Real min_frac{1.e-4};
  int nmat_local;
  Kokkos::parallel_reduce(
      "prePTE", ncell,
      PORTABLE_LAMBDA(const int &iloop, int &maxmat) {
        const int i{offsets_v(iloop) - 1};
        nmat_cell_init(iloop) = 0;
        Real mass_sum = 0.0;
        for (int m = 0; m < nmat; ++m) {
          mass_sum += frac_mass_v(i, m);
        }
        for (int m = 0; m < nmat; ++m) {
          if (frac_mass_v(i, m) / mass_sum > min_frac) nmat_cell_init(iloop) += 1;
        }
        const int nmat_cell_loc = nmat_cell_init(iloop);
        maxmat = nmat_cell_loc > maxmat ? nmat_cell_loc : maxmat;
      },
      Kokkos::Max<int>(nmat_local));
  Kokkos::View<const int *, MemoryTraits<ra>> nmat_cell(nmat_cell_init);
  // set up scratch arrays
  constexpr auto KGlobal = Kokkos::Experimental::UniqueTokenScope::Global;
  Kokkos::Experimental::UniqueToken<DES, KGlobal> tokens{};

  const bool small_loop{tokens.size() > ncell};
  const decltype(tokens)::size_type scratch_size{std::min(tokens.size(), ncell)};
  ScratchV<int> pte_mats("PTE::scratch mats", scratch_size, nmat_local);
  ScratchV<int> pte_idxs("PTE::scratch idxs", scratch_size, nmat_local);
  ScratchV<double> mass_pte("PTE::scratch mass", scratch_size, nmat_local);
  ScratchV<double> sie_pte("PTE::scratch sie", scratch_size, nmat_local);
  ScratchV<double> vfrac_pte("PTE::scratch vfrac", scratch_size, nmat_local);
  ScratchV<double> temp_pte("PTE::scratch temp", scratch_size, nmat_local);
  ScratchV<double> press_pte("PTE::scratch press", scratch_size, nmat_local);
  ScratchV<double> rho_pte("PTE::scratch rho", scratch_size, nmat_local);
  const int pte_solver_scratch_size = PTESolverRhoTRequiredScratch(nmat_local);
  ScratchV<double> solver_scratch("PTE::scratch solver", scratch_size, pte_solver_scratch_size);
  // ScratchV<EOS> eos_pte("PTE::scratch eos", scratch_size, nmat);

  Kokkos::View<int, MemoryTraits<at_int>> res("PTE::num fails");
  Kokkos::View<int, MemoryTraits<at_int>> n_solves("PTE::num solves");
  // loop over cells
  portableFor(
      "PTE", 0, ncell, PORTABLE_LAMBDA(const int &iloop) {
        // cell offset
        const int i{offsets_v(iloop) - 1};
        // get "thread-id" like thing with optimization
        // for small loops
        const int32_t token{tokens.acquire()};
        const int32_t tid{small_loop ? iloop : token};
        // caching mechanism
        Real cache[MAX_NUM_LAMBDAS];
        // calculate total mass in the cell
        Real mass_sum{0.0};
        for (int m = 0; m < nmat; ++m)
          mass_sum += frac_mass_v(i, m);
        // validate specific volume
        spvol_v(i) = r_is_inp ? spvol_v(i) : vol_v(i) / mass_sum;
        // calculate specific internal energy if not already known
        if (!s_is_inp) {
          Real eng_sum{0.0};
          Real vol_sum{0.0};
          for (int m = 0; m < nmat; ++m) {
            Real rho = 1.0 / spvol_v(i);
            Real T = temp_v(i) * ev2k;
            Real press_mat = press_v(i);
            Real spie = sie_v(i);
            Real cv_mat, bmod_mat;
            int mid = eos_offsets_v(m) - 1;
            eos_v(mid).FillEos(rho, T, spie, press_mat, cv_mat, bmod_mat, output, cache);
            eng_sum += spie * frac_mass_v(i, m);
            vol_sum += frac_mass_v(i, m) / rho;
            if (vol_sum < 0.0)
              printf("Neg. Volume:\ninput: %i mat: %i r: %e t: %e s: %e p: %e\n",
                     input_int, mid + 1, rho, T, spie, press_mat);
            // if (!r_is_inp) printf("mat: %i rho lookup: %e\n", eos_offsets_v(m), rho);
          }
          spvol_v(i) = r_is_inp ? spvol_v(i) : vol_sum / mass_sum;
          sie_v(i) = eng_sum / mass_sum;
        }
        // calculate component volumes and energies
        // figure out participating materials, get everything
        // set up to feed in to solver
        int npte = 0;
        double vsum_nopte = 0.0;
        double esum_nopte = 0.0;
        double rhoavg_pte = 0.0;
        for (int m = 0; m < nmat; ++m) {
          const bool something = frac_mass_v(i, m) / mass_sum > min_frac;
          frac_sie_v(i, m) = sie_v(i) * frac_mass_v(i, m);
          frac_vol_v(i, m) = spvol_v(i) * frac_mass_v(i, m);
          if (something) {
            // vfrac true
            int midx = eos_offsets_v(m) - 1;
            pte_mats(tid, npte) = midx;
            pte_idxs(tid, npte) = m;
            mass_pte(tid, npte) = frac_mass_v(i, m);
            sie_pte(tid, npte) = frac_sie_v(i, m) / frac_mass_v(i, m);
            vfrac_pte(tid, npte) = frac_vol_v(i, m);
            npte += 1;
          } else {
            vsum_nopte += frac_vol_v(i, m);
            esum_nopte += frac_sie_v(i, m);
          }
        }
        // call solver to get correct volume fractions
        // component masses, total energy, and total volume are conserved
        const double vfrac_tot = r_is_inp ? spvol_v(i) * mass_sum : vol_v(i);
        const double sie_tot = sie_v(i); // - esum_nopte / rhoavg_pte;
        // calculate initial guess for mixed cells
        if (npte > 1) {
          n_solves() += 1;
          // find the most common material
          int common_mat_id = 0;
          Real max_mass = 0.0;
          for (int m = 0; m < nmat; ++m) {
            if (frac_mass_v(i, m) > max_mass) {
              max_mass = frac_mass_v(i, m);
              common_mat_id = m;
            }
          }
          int common_mid = eos_offsets_v(common_mat_id) - 1;
          // eos lookup of most common material
          Real common_rho{1.0 / spvol_v(i)};
          Real common_temp{0.};
          Real common_sie{frac_sie_v(i, common_mat_id) / frac_mass_v(i, common_mat_id)};
          Real common_press{0.};
          Real common_cv{0.0};
          Real common_bmod{0.0};
          Real ref_p, ref_dpde, ref_dvdt, ref_sie;
          eos_v(common_mid)
              .ValuesAtReferenceState(common_rho, common_temp, ref_sie, ref_p, common_cv,
                                      common_bmod, ref_dpde, ref_dvdt);
          const auto pt_output =
              thermalqs::all_values -
              (thermalqs::density | thermalqs::specific_internal_energy);
          eos_v(common_mid)
              .FillEos(common_rho, common_temp, common_sie, common_press, common_cv,
                       common_bmod, pt_output, cache);
          common_press = common_press < 0.0 ? press_v(i) : common_press;
          // PT lookup of all other materials at PT of common material
          const auto guess_output =
              thermalqs::all_values - (thermalqs::pressure | thermalqs::temperature);
          for (int m = 0; m < npte; ++m) {
            int mid = pte_mats(tid, m);
            int midx = pte_idxs(tid, m);
            Real matrho{0.0}, mat_sie{0.0};
            press_pte(tid, m) = common_press;
            rho_pte(tid, m) = common_rho;
            temp_pte(tid, m) = common_temp;
            sie_pte(tid, m) = common_sie;
            vfrac_pte(tid, m) = frac_mass_v(i, midx) / common_rho;
            if (mid != common_mid) {
              // reuse common cv and bmod b/c the values aren't needed
              eos_v(mid).FillEos(matrho, common_temp, mat_sie, common_press, common_cv,
                                 common_bmod, guess_output, cache);
              // get rid of 0 density issues
              Real guess_density = matrho;
              if (matrho < min_guess_density) {
                Real ref_T, ref_dpde, ref_dvdt;
                eos_v(mid).ValuesAtReferenceState(guess_density, ref_T, ref_sie, ref_p,
                                                  common_cv, common_bmod, ref_dpde,
                                                  ref_dvdt);
                press_pte(tid, m) = ref_p;
                temp_pte(tid, m) = ref_T;
                // sie_pte(tid, m) = ref_sie;
              }
              vfrac_pte(tid, m) = frac_mass_v(i, midx) / guess_density;
              rho_pte(tid, m) = guess_density;
              sie_pte(tid, m) = mat_sie;
            }
          }
          // loop over materials to calc temperature, per material quantities
          // Real sum_mu_cvt {0.0}, sum_mu_e0 {0.0}, sum_mu_cv {0.0};
          // for (int m = 0; m < npte; ++m) {
          //  const int mid = pte_mats(tid, m);
          //  const int midx = pte_idxs(tid, m);
          //  const Real mu = frac_mass_v(i, midx)/mass_sum;
          //  Real ref_rho {0.0}, ref_T {0.0}, ref_sie {0.0}, ref_p{0.0};
          //  Real ref_cv {0.0}, ref_bmod {0.0}, ref_dpde {0.0}, ref_dvdt {0.0};
          //  eos_v(mid).ValuesAtReferenceState(ref_rho, ref_T, ref_sie, ref_p,
          //                                    ref_cv, ref_bmod, ref_dpde, ref_dvdt);
          //  sum_mu_cvt += mu*ref_cv*ref_T;
          //  sum_mu_e0 += mu*ref_sie;
          //  sum_mu_cv += mu*ref_cv;
          //}
          // const Real temp_est {(sie_v(i) + sum_mu_cvt - sum_mu_e0) / sum_mu_cv};
          // const bool temp_est_pos {temp_est > 0.0};
          // Real vol_sum{0.0}, eng_sum{0.0};
          // for (int m = 0; m < npte; ++m) {
          //  const int mid = pte_mats(tid, m);
          //  const int midx = pte_idxs(tid, m);
          //  const Real mu = frac_mass_v(i, midx)/mass_sum;
          //  Real ref_rho {0.0}, ref_T {0.0}, ref_sie {0.0}, ref_p{0.0};
          //  Real ref_cv {0.0}, ref_bmod {0.0}, ref_dpde {0.0}, ref_dvdt {0.0};
          //  eos_v(mid).ValuesAtReferenceState(ref_rho, ref_T, ref_sie, ref_p,
          //                                    ref_cv, ref_bmod, ref_dpde, ref_dvdt);
          //  //const Real delta_T {temp_est_pos ? temp_est - ref_T : 0.0};
          //  //sie_pte(tid, m) = ref_sie;//temp_est_pos ? ref_cv*delta_T + ref_sie :
          //  mu*sie_v(i);
          //  //const Real guess_rho = (1.0 - ref_rho*ref_dvdt*delta_T) <= 0.0
          //  //                     ? ref_rho
          //  //                     : ref_rho*(1.0 - ref_rho*ref_dvdt*delta_T);
          //  vfrac_pte(tid, m) = frac_mass_v(i, midx) / ref_rho;
          //  rho_pte(tid, m) = ref_rho;
          //  temp_pte(tid, m) = ref_T;//temp_v(i);
          //  press_pte(tid, m) = ref_p;
          //  //sie_pte(tid, m) = ref_sie;
          //  //printf("ref vals:\nv r t p s\n%e %e %e %e %e\n",
          //  //       vfrac_pte(tid, m), ref_rho, ref_T, ref_p, ref_sie);
          //  //mass_pte(tid, m) = ref_rho;
          //  //rho_pte(tid, m) = 1.0/spvol_v(i);
          //  //const auto guess_output = thermalqs::all_values
          //  //                        - (thermalqs::pressure | thermalqs::temperature);
          //  //Real mat_sie {frac_sie_v(i, midx) / frac_mass_v(i, midx)};
          //  //Real mat_press {press_v(i)};
          //  //Real mat_cv, mat_bmod, mat_dpde;
          //  //temp_pte(tid, m) = temp_v(i);
          //  //eos_v(mid).FillEos(rho_pte(tid, m), temp_pte(tid, m), mat_sie,
          //  //                   mat_press, mat_cv, mat_bmod, mat_dpde);
          //  //vfrac_pte(tid, m) = frac_mass_v(i, midx) / rho_pte(tid, m);
          //  //press_pte(tid, m) = press_v(i);
          //  //sie_pte(tid, m) = mat_sie;
          //  //vfrac_pte(tid, m) = rho_pte(tid, m) / frac_mass_v(i, midx);
          //  //vol_sum += vfrac_pte(tid, m);
          //  //eng_sum += mu*sie_pte(tid, m);
          //}
        }
        int niter;
        // TODO: this struct declaration should probably be moved elsewhere
        struct EOSAccessor_ {
          PORTABLE_INLINE_FUNCTION
          EOSAccessor_(const Kokkos::View<EOS *, Llft> &eos_v, int *mats)
            : eos_v_(eos_v), mats_(mats) {}
          PORTABLE_INLINE_FUNCTION
          EOS &operator[] (const int m) const {
            return eos_v_(mats_[m]);
          }
          Kokkos::View<EOS *, Llft> eos_v_;
          int *mats_;
        };
        EOSAccessor_ eos_inx(eos_v, &pte_mats(tid, 0));
        PTESolverRhoT<EOSAccessor_, Real *, Real **>
          method(npte, eos_inx, vfrac_tot, sie_tot, &rho_pte(tid, 0), &vfrac_pte(tid, 0),
                 &sie_pte(tid, 0), &temp_pte(tid, 0), &press_pte(tid, 0), lambda_map,
                 &solver_scratch(tid, 0));
        const bool res_{PTESolver(method)};
        //const bool res_{pte_closure_josh_offset(
        //    npte, eos_v.data(), vfrac_tot, sie_tot, &pte_mats(tid, 0), &rho_pte(tid, 0),
        //    &vfrac_pte(tid, 0), &sie_pte(tid, 0), &temp_pte(tid, 0), &press_pte(tid, 0),
        //    lambda_map)};
        // assign local values to global
        press_v(i) = p_is_inp ? press_v(i) : 0.0;
        sie_v(i) = s_is_inp ? sie_v(i) : 0.0;
        temp_v(i) = t_is_inp ? temp_v(i) : 0.0;
        bmod_v(i) = 0.0;
        cv_v(i) = 0.0;
        dpde_v(i) = 0.0;
        Real vol_sum = 0.0;
        for (int m = 0; m < npte; ++m) {
          int mid = pte_mats(tid, m);
          int midx = pte_idxs(tid, m);
          Real rho = mass_pte(tid, m) / vfrac_pte(tid, m);
          Real spie = sie_pte(tid, m);
          Real press_mat = press_v(i);
          Real T = temp_v(i) * ev2k;
          Real cv_mat = 0.0;
          Real bmod_mat = 0.0;
          eos_v(mid).FillEos(rho, T, spie, press_mat, cv_mat, bmod_mat, output, cache);
          frac_sie_v(i, midx) = spie * frac_mass_v(i, midx);
          frac_vol_v(i, midx) = frac_mass_v(i, midx) / rho;
          vol_sum += frac_vol_v(i, midx);
          temp_v(i) += t_is_inp ? 0.0 : frac_vol_v(i, midx) / ev2k * T;
          press_v(i) += p_is_inp ? 0.0 : frac_vol_v(i, midx) * press_mat;
          // optional cv and bmod
          const Real dpde_mat =
              rho * eos_v(mid).GruneisenParamFromDensityInternalEnergy(rho, spie, cache);
          sie_v(i) += s_is_inp ? 0.0 : frac_sie_v(i, midx);
          bmod_v(i) += frac_vol_v(i, midx) * bmod_mat;
          cv_v(i) += frac_mass_v(i, midx) * cv_mat;
          dpde_v(i) += frac_vol_v(i, midx) * dpde_mat;
          if (do_frac_bmod) {
            frac_bmod_v(i, midx) = bmod_mat;
          }
          if (do_frac_cv) {
            frac_cv_v(i, midx) = cv_mat;
          }
          if (do_frac_dpde) {
            frac_dpde_v(i, midx) = dpde_mat;
          }
        }
        press_v(i) /= p_is_inp ? 1.0 : vol_sum;
        sie_v(i) /= s_is_inp ? 1.0 : mass_sum;
        temp_v(i) /= t_is_inp ? 1.0 : vol_sum;
        bmod_v(i) /= vol_sum;
        cv_v(i) *= ev2k / mass_sum;
        dpde_v(i) /= vol_sum;
        spvol_v(i) = r_is_inp ? spvol_v(i) : vol_sum / mass_sum;
        pmax_v(i) = press_v(i) > pmax_v(i) ? press_v(i) : pmax_v(i);
        // if failed, increment return value by 1
        // use if rather than ternary to reduce contention
        if (!res_) {
          res() += 1;
        }
        // release token
        tokens.release(token);
      });
  Kokkos::fence();
  // copy results back into local values
  deep_copy(ret, res);
  deep_copy(press_hv, press_v);
  deep_copy(pmax_hv, pmax_v);
  deep_copy(vol_hv, vol_v);
  deep_copy(spvol_hv, spvol_v);
  deep_copy(sie_hv, sie_v);
  deep_copy(temp_hv, temp_v);
  deep_copy(bmod_hv, bmod_v);
  deep_copy(dpde_hv, dpde_v);
  deep_copy(cv_hv, cv_v);
  deep_copy(frac_mass_hv, frac_mass_v);
  deep_copy(frac_vol_hv, frac_vol_v);
  deep_copy(frac_sie_hv, frac_sie_v);
  if (do_frac_bmod) {
    deep_copy(frac_bmod_hv, frac_bmod_v);
  }
  if (do_frac_dpde) {
    deep_copy(frac_dpde_hv, frac_dpde_v);
  }
  if (do_frac_cv) {
    deep_copy(frac_cv_hv, frac_cv_v);
  }
  portableFor(
      "freeing lambda ptrs", 0, nmat,
      PORTABLE_LAMBDA(const int &i) { free(lambda_map[i]); });
  PORTABLE_FREE(lambda_map);
  int tot_solves{0};
  deep_copy(tot_solves, n_solves);
  if (ret > 0) printf("%i/%i solves failed\n", ret, tot_solves);
#endif // PORTABILITY_STRATEGY_KOKKOS
  return ret;
}

int finalize_sg_eos(const int nmat, EOS *&eos, const int own_kokkos) {
  {
    for (int i = 0; i < nmat; ++i)
      eos[i].Finalize();
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif // PORTABILITY_STRATEGY_KOKKOS
  }
  delete[] eos;
  eos = nullptr;
#ifdef PORTABILITY_STRATEGY_KOKKOS
  if (own_kokkos != 0) Kokkos::finalize();
#endif
  return 0;
}
