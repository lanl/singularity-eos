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

// apply everything but ramp in order to possibly calculate the
// SAP ramp parameters from p-alhpa ramp parameters
#define SGAPPLYMODSIMPLE(A)                                                              \
  EOSBuilder::applyShiftAndScale(A, enabled[0] == 1, enabled[1] == 1, vals[0], vals[1])

#define SGAPPLYMOD(A)                                                                    \
  EOSBuilder::applyShiftAndScaleAndBilinearRamp(                                         \
      A, enabled[0] == 1, enabled[1] == 1, enabled[2] == 1 || enabled[3] == 1, vals[0],  \
      vals[1], vals[2], vals[3], vals[4], vals[5])

int def_en[4] = {0, 0, 0, 0};
double def_v[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

int init_sg_IdealGas(const int matindex, EOS *eos, const double gm1, const double Cv,
                     int const *const enabled, double *const vals) {
  assert(matindex >= 0);
  EOS eosi = SGAPPLYMODSIMPLE(IdealGas(gm1, Cv));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
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
                      const double Cv, int const *const enabled, double *const vals) {
  assert(matindex >= 0);
  EOS eosi = SGAPPLYMODSIMPLE(Gruneisen(C0, s1, s2, s3, G0, b, rho0, T0, P0, Cv));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
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
                const double Cv, int const *const enabled, double *const vals) {
  assert(matindex >= 0);
  EOS eosi = SGAPPLYMODSIMPLE(JWL(A, B, R1, R2, w, rho0, Cv));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
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
                          int const *const enabled, double *const vals) {
  assert(matindex >= 0);
  EOS eosi = SGAPPLYMODSIMPLE(DavisProducts(a, b, k, n, vc, pc, Cv, E0));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
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
                           double *const vals) {
  assert(matindex >= 0);
  EOS eosi =
      SGAPPLYMODSIMPLE(DavisReactants(rho0, e0, P0, T0, A, B, C, G0, Z, alpha, Cv0));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
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
                              double *const vals) {
  assert(matindex >= 0);
  EOS eosi = SGAPPLYMODSIMPLE(SpinerEOSDependsRhoT(std::string(filename), matid));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
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
                                double *const vals) {
  assert(matindex >= 0);
  EOS eosi = SGAPPLYMODSIMPLE(SpinerEOSDependsRhoSie(std::string(filename), matid));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
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
                   double *const vals) {
  assert(matindex >= 0);
  EOS eosi = SGAPPLYMODSIMPLE(EOSPAC(id));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
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
using Llft = Kokkos::LayoutLeft;
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
// EAP centric arguments and function signature

#include <signal.h>
#include <execinfo.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <algorithm>

static inline void printBacktrace() {
  static std::vector<std::string> bt_vec;
    void *stackBuffer[64];
    int numAddresses = backtrace((void**) &stackBuffer, 64); 
    char **addresses = backtrace_symbols(stackBuffer, numAddresses);
    std::stringstream s;
    
    for( int i = 0 ; i < numAddresses ; ++i ) {
      s << "[" << i << "]: " << addresses[i] << "\n";
      //printf("[%2d]: %s\n", i, addresses[i]); 
    }
    const bool stack_exists = std::find(bt_vec.begin(), bt_vec.end(), s.str())
      != bt_vec.end();
    if(!stack_exists) {
      bt_vec.push_back(s.str());
      std::cout << s.str();
    }
    free(addresses); 
} 

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
  //printBacktrace();
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

  // add characterization of solver convergence
  Kokkos::View<int> min_its("min_its"), max_its("max_its"), sum_its("sum_its");
  bool do_its_analysis = input_int != 0 || (input_int == 0 && ncell > 1);
  do_its_analysis = do_its_analysis && input_int != -1;
  do_its_analysis = do_its_analysis && nmat > 1;
  if (do_its_analysis) {
    deep_copy(min_its, int(1e6));
  }
  Kokkos::View<int, MemoryTraits<at_int>> res("PTE::num fails");
  Kokkos::View<int, MemoryTraits<at_int>> n_solves("PTE::num solves");
  const std::string perf_nums =
      "[" + std::to_string(nmat) + "," + std::to_string(ncell) + "]";
  switch (input_int) {
  case -3: {
    // T-rho input
    // set frac_vol = 1/nmat
    // set rho_i = nmat / spvol * frac_mass_i
    // iterate PTE solver to obtain internal energies
    // that results in the input T
    pte_solver_scratch_size = PTESolverFixedTRequiredScratch(nmat);
    solver_scratch = ScratchV<double>(VAWI("PTE::scratch solver"), scratch_size,
                                      pte_solver_scratch_size);
    const std::string rt_name = "PTE::solve (rho,T) input" + perf_nums;
    portableFor(
        rt_name.c_str(), 0, ncell, PORTABLE_LAMBDA(const int &iloop) {
          // cell offset
          const int i{offsets_v(iloop) - 1};
          // get "thread-id" like thing with optimization
          // for small loops
          const int32_t token{tokens.acquire()};
          const int32_t tid{small_loop ? iloop : token};
          double mass_sum{0.0};
          // normalize mass fractions
          // first find the mass sum
          // also set idxs as the decrement of the eos offsets
          // to take into account 1 based indexing in fortran
          for (int m = 0; m < nmat; ++m) {
            mass_sum += frac_mass_v(i, m);
            pte_idxs(tid, m) = eos_offsets_v(m) - 1;
            frac_vol_v(i, m) = 0.0;
          }
          for (int m = 0; m < nmat; ++m) {
            frac_mass_v(i, m) /= mass_sum;
          }
          // set inputs
          int npte = 0;
          for (int m = 0; m < nmat; ++m) {
            if (frac_mass_v(i, m) > 1.e-12) {
              pte_idxs(tid, npte) = eos_offsets_v(m) - 1;
              pte_mats(tid, npte) = m;
              npte += 1;
            }
            vfrac_pte(tid, m) = 0.0;
            sie_pte(tid, m) = 0.0;
            temp_pte(tid, m) = 0.0;
            press_pte(tid, m) = 0.0;
          }

          for (int mp = 0; mp < npte; ++mp) {
            const int m = pte_mats(tid, mp);
            rho_pte(tid, mp) = npte / spvol_v(i) * frac_mass_v(i, m);
            vfrac_pte(tid, mp) = 1.0 / npte;
            temp_pte(tid, mp) = temp_v(i) * ev2k;
          }
          Real sie_tot_true{0.0};
          const int neq = npte;
          singularity::mix_impl::CacheAccessor cache(&solver_scratch(tid, 0) +
                                                     neq * (neq + 4) + 2 * npte);
          if (npte > 1) {
            // create solver lambda
            // eos accessor
            singularity::EOSAccessor_ eos_inx(eos_v, &pte_idxs(tid, 0));
            PTESolverFixedT<singularity::EOSAccessor_, Real *, Real **> method(
                npte, eos_inx, 1.0, temp_pte(tid, 0), &rho_pte(tid, 0),
                &vfrac_pte(tid, 0), &sie_pte(tid, 0), &temp_pte(tid, 0),
                &press_pte(tid, 0), cache, &solver_scratch(tid, 0));
            const bool res_{PTESolver(method)};
            // calculate total internal energy
            for (int mp = 0; mp < npte; ++mp) {
              const int m = pte_mats(tid, mp);
              sie_tot_true += sie_pte(tid, mp) * frac_mass_v(i, m);
            }
	    min_its() = std::min(min_its(), method.Niter());
	    max_its() = std::max(max_its(), method.Niter());
	    sum_its() += method.Niter();
          } else {
            // pure cell (nmat = 1)
            // calculate sie from single eos
            sie_pte(tid, 0) = eos_v(pte_idxs(tid, 0))
                                  .InternalEnergyFromDensityTemperature(
                                      rho_pte(tid, 0), temp_pte(tid, 0), cache[0]);
            // set total sie to material 0 value
            sie_tot_true = sie_pte(tid, 0);
            // set pressure
            press_pte(tid, 0) = eos_v(pte_idxs(tid, 0))
                                    .PressureFromDensityTemperature(
                                        rho_pte(tid, 0), temp_pte(tid, 0), cache[0]);
          }
          // since rho_ptes are correct, now do rho-T lookups for remaining quantities
          // first zero out quantities to be averaged
          press_v(i) = 0.0;
          bmod_v(i) = 0.0;
          cv_v(i) = 0.0;
          dpde_v(i) = 0.0;
          // total sie is known
          sie_v(i) = sie_tot_true;
          for (int mp = 0; mp < npte; ++mp) {
            const int m = pte_mats(tid, mp);
            // pressure contribution from material m
            press_v(i) += press_pte(tid, mp) * vfrac_pte(tid, mp);
	    const Real press_eos = eos_v(pte_idxs(tid, mp))
	      .PressureFromDensityInternalEnergy(rho_pte(tid, mp), sie_pte(tid, mp), cache[mp]);
	    Real rel_press = std::abs(press_pte(tid, mp) - press_eos);
	    rel_press /= std::abs(press_pte(tid, mp) + press_eos);
	    //if (rel_press > 1.e-3) {
	      printf("pressures %i: %e %e, rho: %e sie:%e\n", pte_idxs(tid, mp),
		     press_pte(tid, mp), press_eos, rho_pte(tid, mp), sie_pte(tid, mp));
	      //}
            // assign per material specific internal energy
            frac_sie_v(i, m) = sie_pte(tid, mp);
            // assign volume fraction based on pte calculation
            frac_vol_v(i, m) = vfrac_pte(tid, mp) * vol_v(i);
            // calculate bulk modulus for material m
            const Real bmod_m = eos_v(pte_idxs(tid, mp))
                                    .BulkModulusFromDensityTemperature(
                                        rho_pte(tid, mp), temp_pte(tid, mp), cache[mp]);
            // add bmod contribution from material m
            bmod_v(i) += bmod_m * vfrac_pte(tid, mp);
            // calculate specific heat for material m
            const Real cv_m =
                ev2k * eos_v(pte_idxs(tid, mp))
                           .SpecificHeatFromDensityTemperature(
                               rho_pte(tid, mp), temp_pte(tid, mp), cache[mp]);
            // add mass weighted contribution specific heat for material m
            cv_v(i) += cv_m * frac_mass_v(i, m);
            // calculate gruneisen parameter for material m
            const Real dpde_m = eos_v(pte_idxs(tid, mp))
                                    .GruneisenParamFromDensityTemperature(
                                        rho_pte(tid, mp), temp_pte(tid, mp), cache[mp]);
            // add gruneisen param contribution from material m
            dpde_v(i) += dpde_m * vfrac_pte(tid, mp);
            // optionally assign per material quantities to per material arrays
            if (do_frac_bmod) {
              frac_bmod_v(i, m) = bmod_m;
            }
            if (do_frac_cv) {
              frac_cv_v(i, m) = cv_m;
            }
            if (do_frac_dpde) {
              frac_dpde_v(i, m) = dpde_m;
            }
          }
          for (int m = 0; m < nmat; ++m) {
            frac_mass_v(i, m) *= mass_sum;
          }
	  for (int mp = 0; mp < npte; ++mp) {
	    Real rel_press = std::abs(press_pte(tid, mp) - press_v(i));
	    rel_press /= std::abs(press_pte(tid, mp) + press_v(i));
	    if (rel_press > 1.e-3) {
	      printf("pressures disagree %i: %e %e\n", pte_idxs(tid, mp),
		     press_pte(tid, mp), press_v(i));
	    }
	  }
	  printf("sie total: %e\n", sie_v(i));
          // assign max pressure
          pmax_v(i) = press_v(i) > pmax_v(i) ? press_v(i) : pmax_v(i);
          // release the token used for scratch arrays
          tokens.release(token);
        });
    break;
  }
  case -2: {
    // rho-P input
    // set frac_vol = 1/nmat
    // set rho_i = nmat / spvol * frac_mass_i
    // iterate PTE solver to obtain internal energies
    // that results in the input P
    pte_solver_scratch_size = PTESolverRhoTRequiredScratch(nmat);
    solver_scratch = ScratchV<double>(VAWI("PTE::scratch solver"), scratch_size,
                                      pte_solver_scratch_size);
    const std::string rp_name = "PTE::solve (rho,P) input" + perf_nums;
    portableFor(
        rp_name.c_str(), 0, ncell, PORTABLE_LAMBDA(const int &iloop) {
          // cell offset
          const int i{offsets_v(iloop) - 1};
          // get "thread-id" like thing with optimization
          // for small loops
          const int32_t token{tokens.acquire()};
          const int32_t tid{small_loop ? iloop : token};
          double mass_sum{0.0};
          // normalize mass fractions
          // first find the mass sum
          // also set idxs as the decrement of the eos offsets
          // to take into account 1 based indexing in fortran
          for (int m = 0; m < nmat; ++m) {
            mass_sum += frac_mass_v(i, m);
            pte_idxs(tid, m) = eos_offsets_v(m) - 1;
            frac_vol_v(i, m) = 0.0;
          }
          for (int m = 0; m < nmat; ++m) {
            frac_mass_v(i, m) /= mass_sum;
          }
          // set inputs
          int npte = 0;
          for (int m = 0; m < nmat; ++m) {
            if (frac_mass_v(i, m) > 1.e-12) {
              pte_idxs(tid, npte) = eos_offsets_v(m) - 1;
              pte_mats(tid, npte) = m;
              npte += 1;
            }
            vfrac_pte(tid, m) = 0.0;
            sie_pte(tid, m) = 0.0;
            press_pte(tid, m) = press_v(i);
            temp_pte(tid, m) = 0.0;
          }

          for (int mp = 0; mp < npte; ++mp) {
            const int m = pte_mats(tid, mp);
            rho_pte(tid, mp) = npte / spvol_v(i) * frac_mass_v(i, m);
            vfrac_pte(tid, mp) = 1.0 / npte;
          }
          Real sie_tot_true{0.0};
          const int neq = npte + 1;
          singularity::mix_impl::CacheAccessor cache(&solver_scratch(tid, 0) +
                                                     neq * (neq + 4) + 2 * npte);
          if (npte > 1) {
            // create solver lambda
            // eos accessor
            singularity::EOSAccessor_ eos_inx(eos_v, &pte_idxs(tid, 0));
            PTESolverFixedP<singularity::EOSAccessor_, Real *, Real **> method(
                npte, eos_inx, 1.0, press_pte(tid, 0), &rho_pte(tid, 0),
                &vfrac_pte(tid, 0), &sie_pte(tid, 0), &temp_pte(tid, 0),
                &press_pte(tid, 0), cache[0], &solver_scratch(tid, 0));
            const bool res_{PTESolver(method)};
            // calculate total sie
            for (int mp = 0; mp < npte; ++mp) {
              const int m = pte_mats(tid, mp);
              sie_tot_true += sie_pte(tid, mp) * frac_mass_v(i, m);
            }
	    min_its() = std::min(min_its(), method.Niter());
	    max_its() = std::max(max_its(), method.Niter());
	    sum_its() += method.Niter();
          } else {
            // pure cell (nmat = 1)
            // calculate sie from single eos
            auto p_from_t = [&](const Real &t_i) {
              return eos_v(pte_idxs(tid, 0))
                  .PressureFromDensityTemperature(rho_pte(tid, 0), t_i, cache[0]);
            };
            // calculate sie root bounds
            Real r_rho{}, r_temp{}, r_sie{}, r_press{}, r_cv{}, r_bmod{};
            Real r_dpde{}, r_dvdt{};
            eos_v(pte_idxs(tid, 0))
                .ValuesAtReferenceState(r_rho, r_temp, r_sie, r_press, r_cv, r_bmod,
                                        r_dpde, r_dvdt);
            RootFinding1D::RootCounts co{};
            Real temp_i;
            RootFinding1D::findRoot(p_from_t, press_v(i), r_temp, 1.e-10 * r_temp,
                                    1.e10 * r_temp, 1.e-12, 1.e-12, temp_i, co);
            sie_pte(tid, 0) = eos_v(pte_idxs(tid, 0))
                                  .InternalEnergyFromDensityTemperature(rho_pte(tid, 0),
                                                                        temp_i, cache[0]);
            sie_tot_true = sie_pte(tid, 0);
            // set temperature
            temp_pte(tid, 0) = temp_i;
          }
          // since rho_ptes are correct, now do rho-T lookups for remaining quantities
          // first zero out quantities to be averaged
          temp_v(i) = 0.0;
          bmod_v(i) = 0.0;
          cv_v(i) = 0.0;
          dpde_v(i) = 0.0;
          // total sie is known
          sie_v(i) = sie_tot_true;
          for (int mp = 0; mp < npte; ++mp) {
            const int m = pte_mats(tid, mp);
            temp_v(i) += temp_pte(tid, mp) * vfrac_pte(tid, mp);
            // assign per material specific internal energy
            frac_sie_v(i, m) = sie_pte(tid, mp);
            // assign volume fraction based on pte calculation
            frac_vol_v(i, m) = vfrac_pte(tid, mp) * vol_v(i);
            // calculate bulk modulus for material m
            const Real bmod_m = eos_v(pte_idxs(tid, mp))
                                    .BulkModulusFromDensityTemperature(
                                        rho_pte(tid, mp), temp_pte(tid, mp), cache[mp]);
            // add bmod contribution from material m
            bmod_v(i) += bmod_m * vfrac_pte(tid, mp);
            // calculate specific heat for material m
            const Real cv_m =
                ev2k * eos_v(pte_idxs(tid, mp))
                           .SpecificHeatFromDensityTemperature(
                               rho_pte(tid, mp), temp_pte(tid, mp), cache[mp]);
            // add mass weighted contribution specific heat for material m
            cv_v(i) += cv_m * frac_mass_v(i, m);
            // calculate gruneisen parameter for material m
            const Real dpde_m = eos_v(pte_idxs(tid, mp))
                                    .GruneisenParamFromDensityTemperature(
                                        rho_pte(tid, mp), temp_pte(tid, mp), cache[mp]);
            // add gruneisen param contribution from material m
            dpde_v(i) += dpde_m * vfrac_pte(tid, mp);
            // optionally assign per material quantities to per material arrays
            if (do_frac_bmod) {
              frac_bmod_v(i, m) = bmod_m;
            }
            if (do_frac_cv) {
              frac_cv_v(i, m) = cv_m;
            }
            if (do_frac_dpde) {
              frac_dpde_v(i, m) = dpde_m;
            }
          }
          for (int m = 0; m < nmat; ++m) {
            frac_mass_v(i, m) *= mass_sum;
          }
          // convert to ev
          temp_v(i) /= ev2k;
          // assign max pressure
          pmax_v(i) = press_v(i) > pmax_v(i) ? press_v(i) : pmax_v(i);
          // release the token used for scratch arrays
          tokens.release(token);
        });
    break;
  }
  case -1: {
    // P-T input
    const std::string pt_name = "PTE::solve (P,T) input" + perf_nums;
    portableFor(
        pt_name.c_str(), 0, ncell, PORTABLE_LAMBDA(const int &iloop) {
          // cell offset
          const int i{offsets_v(iloop) - 1};
          // get "thread-id" like thing with optimization
          // for small loops
          const int32_t token{tokens.acquire()};
          const int32_t tid{small_loop ? iloop : token};
          // caching mechanism
          Real cache[MAX_NUM_LAMBDAS];
          double mass_sum{0.0};
          // normalize mass fractions
          // first find the mass sum
          // also set idxs as the decrement of the eos offsets
          // to take into account 1 based indexing in fortran
          for (int m = 0; m < nmat; ++m) {
            mass_sum += frac_mass_v(i, m);
            pte_idxs(tid, m) = eos_offsets_v(m) - 1;
            temp_pte(tid, m) = temp_v(i) * ev2k;
            press_pte(tid, m) = press_v(i);
          }
          for (int m = 0; m < nmat; ++m) {
            frac_mass_v(i, m) /= mass_sum;
          }

          // do r-e of pt for each mat
          singularity::EOSAccessor_ eos_inx(eos_v, &pte_idxs(tid, 0));
          Real vfrac_tot{0.0};
          Real sie_tot{0.0};
          for (int m = 0; m < nmat; ++m) {
            // obtain rho and sie from P-T
            eos_inx[m].DensityEnergyFromPressureTemperature(
                press_pte(tid, m), temp_pte(tid, m), cache, rho_pte(tid, m),
                sie_pte(tid, m));
            // assign volume fractions
            // this is a physical volume
            vfrac_pte(tid, m) = frac_mass_v(i, m) / rho_pte(tid, m) * mass_sum;
            vfrac_tot += vfrac_pte(tid, m);
            // add internal energy component
            sie_tot += sie_pte(tid, m) * frac_mass_v(i, m);
          }
          // since rho_ptes are correct, now do rho-T lookups for remaining quantities
          // first zero out quantities to be averaged
          bmod_v(i) = 0.0;
          cv_v(i) = 0.0;
          dpde_v(i) = 0.0;
          // total sie is known
          sie_v(i) = sie_tot;
          vol_v(i) = vfrac_tot;
          spvol_v(i) = vol_v(i) / mass_sum;
          for (int m = 0; m < nmat; ++m) {
            // assign per material specific internal energy
            frac_sie_v(i, m) = sie_pte(tid, m);
            // apply normalization to to vfrac
            vfrac_pte(tid, m) /= vfrac_tot;
            // assign volume fraction based on pte calculation
            frac_vol_v(i, m) = vfrac_pte(tid, m) * vol_v(i);
            // calculate bulk modulus for material m
            const Real bmod_m = eos_v(pte_idxs(tid, m))
                                    .BulkModulusFromDensityTemperature(
                                        rho_pte(tid, m), temp_pte(tid, m), cache);
            // add bmod contribution from material m
            bmod_v(i) += bmod_m * vfrac_pte(tid, m);
            // calculate specific heat for material m
            const Real cv_m = ev2k * eos_v(pte_idxs(tid, m))
                                         .SpecificHeatFromDensityTemperature(
                                             rho_pte(tid, m), temp_pte(tid, m), cache);
            // add mass weighted contribution specific heat for material m
            cv_v(i) += cv_m * frac_mass_v(i, m);
            frac_mass_v(i, m) *= mass_sum;
            // calculate gruneisen parameter for material m
            const Real dpde_m = eos_v(pte_idxs(tid, m))
                                    .GruneisenParamFromDensityTemperature(
                                        rho_pte(tid, m), temp_pte(tid, m), cache);
            // add gruneisen param contribution from material m
            dpde_v(i) += dpde_m * vfrac_pte(tid, m);
            // optionally assign per material quantities to per material arrays
            if (do_frac_bmod) {
              frac_bmod_v(i, m) = bmod_m;
            }
            if (do_frac_cv) {
              frac_cv_v(i, m) = cv_m;
            }
            if (do_frac_dpde) {
              frac_dpde_v(i, m) = dpde_m;
            }
          }
          // assign max pressure
          pmax_v(i) = press_v(i) > pmax_v(i) ? press_v(i) : pmax_v(i);
          // release the token used for scratch arrays
          tokens.release(token);
        });
    break;
  }
  case 0:
    // rho-sie input
    // no break so fallthrough to case 1
  case 1: {
    // rho-sie input
    pte_solver_scratch_size = PTESolverRhoTRequiredScratch(nmat);
    solver_scratch = ScratchV<double>(VAWI("PTE::scratch solver"), scratch_size,
                                      pte_solver_scratch_size);
    const std::string re_name = "PTE::solve (rho,e) input" + perf_nums;
    portableFor(
        re_name.c_str(), 0, ncell, PORTABLE_LAMBDA(const int &iloop) {
          // cell offset
          const int i{offsets_v(iloop) - 1};
          // get "thread-id" like thing with optimization
          // for small loops
          const int32_t token{tokens.acquire()};
          const int32_t tid{small_loop ? iloop : token};
          double mass_sum{0.0};
          // normalize mass fractions
          // first find the mass sum
          // also set idxs as the decrement of the eos offsets
          // to take into account 1 based indexing in fortran
          for (int m = 0; m < nmat; ++m) {
            mass_sum += frac_mass_v(i, m);
            pte_idxs(tid, m) = eos_offsets_v(m) - 1;
          }
          for (int m = 0; m < nmat; ++m) {
            frac_mass_v(i, m) /= mass_sum;
          }
          // set inputs
          int npte = 0;
          for (int m = 0; m < nmat; ++m) {
            if (frac_mass_v(i, m) > 1.e-12) {
              pte_idxs(tid, npte) = eos_offsets_v(m) - 1;
              pte_mats(tid, npte) = m;
              npte += 1;
            }
            vfrac_pte(tid, m) = 0.0;
            sie_pte(tid, m) = 0.0;
            temp_pte(tid, m) = 0.0;
            press_pte(tid, m) = 0.0;
          }
          for (int mp = 0; mp < npte; ++mp) {
            const int m = pte_mats(tid, mp);
            rho_pte(tid, mp) = npte / spvol_v(i) * frac_mass_v(i, m);
            vfrac_pte(tid, mp) = 1.0 / npte;
            sie_pte(tid, mp) = sie_v(i) * frac_mass_v(i, m);
            temp_pte(tid, mp) = 0.0;
            press_pte(tid, mp) = 0.0;
          }
          const int neq = npte + 1;
          singularity::mix_impl::CacheAccessor cache(&solver_scratch(tid, 0) +
                                                     neq * (neq + 4) + 2 * npte);
          if (npte > 1) {
            // create solver lambda
            // eos accessor
            singularity::EOSAccessor_ eos_inx(eos_v, &pte_idxs(tid, 0));
            // reset inputs
            PTESolverRhoT<singularity::EOSAccessor_, Real *, Real **> method(
                npte, eos_inx, 1.0, sie_v(i), &rho_pte(tid, 0), &vfrac_pte(tid, 0),
                &sie_pte(tid, 0), &temp_pte(tid, 0), &press_pte(tid, 0), cache,
                &solver_scratch(tid, 0));
            const bool res_{PTESolver(method)};
	    min_its() = std::min(min_its(), method.Niter());
	    max_its() = std::max(max_its(), method.Niter());
	    sum_its() += method.Niter();
          } else {
            // pure cell (nmat = 1)
            // calculate sie from single eos
            Real dpdr_m, dtdr_m, dpde_m, dtde_m;
            eos_v(pte_idxs(tid, 0))
                .PTofRE(rho_pte(tid, 0), sie_pte(tid, 0), cache[0], press_pte(tid, 0),
                        temp_pte(tid, 0), dpdr_m, dpde_m, dtdr_m, dtde_m);
          }
          // since rho_ptes are correct, now do rho-T lookups for remaining quantities
          // first zero out quantities to be averaged
          press_v(i) = 0.0;
          temp_v(i) = 0.0;
          bmod_v(i) = 0.0;
          cv_v(i) = 0.0;
          dpde_v(i) = 0.0;
          for (int mp = 0; mp < npte; ++mp) {
            const int m = pte_mats(tid, mp);
            // pressure contribution from material m
            press_v(i) += press_pte(tid, mp) * vfrac_pte(tid, mp);
	    const Real press_eos = eos_v(pte_idxs(tid, mp))
	      .PressureFromDensityInternalEnergy(rho_pte(tid, mp), sie_pte(tid, mp), cache[mp]);
	    Real rel_press = std::abs(press_pte(tid, mp) - press_eos);
	    rel_press /= std::abs(press_pte(tid, mp) + press_eos);
	    if (rel_press > 1.e-3) {
	      printf("p disagree %i: %e %e, rho: %e sie:%e\n", pte_idxs(tid, mp),
		     press_pte(tid, mp), press_eos, rho_pte(tid, mp), sie_pte(tid, mp));
	    }
            // temperature averaging
            temp_v(i) += temp_pte(tid, mp) * vfrac_pte(tid, mp);
            // assign per material specific internal energy
            frac_sie_v(i, m) = sie_pte(tid, mp);
            // assign volume fraction based on pte calculation
            frac_vol_v(i, m) = vfrac_pte(tid, mp) * vol_v(i);
            // calculate bulk modulus for material m
            const Real bmod_m = eos_v(pte_idxs(tid, mp))
                                    .BulkModulusFromDensityTemperature(
                                        rho_pte(tid, mp), temp_pte(tid, mp), cache[mp]);
            // add bmod contribution from material m
            bmod_v(i) += bmod_m * vfrac_pte(tid, mp);
            // calculate specific heat for material m
            const Real cv_m =
                ev2k * eos_v(pte_idxs(tid, mp))
                           .SpecificHeatFromDensityTemperature(
                               rho_pte(tid, mp), temp_pte(tid, mp), cache[mp]);
            // add mass weighted contribution specific heat for material m
            cv_v(i) += cv_m * frac_mass_v(i, m);
            // calculate gruneisen parameter for material m
            const Real dpde_m = eos_v(pte_idxs(tid, mp))
                                    .GruneisenParamFromDensityTemperature(
                                        rho_pte(tid, mp), temp_pte(tid, mp), cache[mp]);
            // add gruneisen param contribution from material m
            dpde_v(i) += dpde_m * vfrac_pte(tid, mp);
            // optionally assign per material quantities to per material arrays
            if (do_frac_bmod) {
              frac_bmod_v(i, m) = bmod_m;
            }
            if (do_frac_cv) {
              frac_cv_v(i, m) = cv_m;
            }
            if (do_frac_dpde) {
              frac_dpde_v(i, m) = dpde_m;
            }
          }
          for (int m = 0; m < nmat; ++m) {
            frac_mass_v(i, m) *= mass_sum;
          }
	  for (int mp = 0; mp < npte; ++mp) {
	    Real rel_press = std::abs(press_pte(tid, mp) - press_v(i));
	    rel_press /= std::abs(press_pte(tid, mp) + press_v(i));
	    if (rel_press > 1.e-3 || press_v(i) < 0.0) {
	      printf("--p disagree %i: %e %e\n", pte_idxs(tid, mp),
		     press_pte(tid, mp), press_v(i));
	    }
	  }
          // convert temperature to ev from kelvin
          temp_v(i) /= ev2k;
          // assign max pressure
          pmax_v(i) = press_v(i) > pmax_v(i) ? press_v(i) : pmax_v(i);
          // release the token used for scratch arrays
          tokens.release(token);
        });
    break;
  }
  }

  Kokkos::fence();
  if(do_its_analysis && false) {
    printf("\n---\nmin_its = %i, max_its = %i, avg_its = %e\n---\n",
	   min_its(), max_its(), (double)(sum_its()) / ncell);
  }
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
  deep_copy(frac_sie_hv, frac_sie_v);
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
//  if(input_int != 0) {
//    printf("inp: %i p: %.3e t: %.3e r: %.3e s: %.3e\n", input_int, press_hv(0),
//    temp_hv(0),
//	   1.0/spvol_hv(0), sie_hv(0));
//    for (int m = 0; m < nmat; ++m) {
//      printf("mf: %.3e vf: %.3e\n", frac_mass_hv(0,m), frac_vol_hv(0,m));
//    }
//  }
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
