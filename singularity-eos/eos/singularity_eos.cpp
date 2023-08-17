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

#include <cassert>
#include <ports-of-call/portability.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_builder.hpp>
#include <singularity-eos/eos/singularity_eos.hpp>

using namespace singularity;

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


int init_sg_SAP_Polynomial(const int matindex, EOS *eos, const double rho0,
                           const double a0, const double a1, const double a2c,
                           const double a2e, const double a3, const double b0,
                           const double b1, const double b2c, const double b2e,
                           const double b3, int const *const enabled,
                           double *const vals) {
  assert(matindex >= 0);
  EOS eosi =
      SGAPPLYMODSIMPLE(SAP_Polynomial(rho0, a0, a1, a2c, a2e, a3, b0, b1, b2c, b2e, b3));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
  EOS eos_ = SGAPPLYMOD(SAP_Polynomial(rho0, a0, a1, a2c, a2e, a3, b0, b1, b2c, b2e, b3));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}

int init_sg_StiffGas(const int matindex, EOS *eos, const double gm1, const double Cv,
                     const double Pinf, const double qq, int const *const enabled,
                     double *const vals) {
  assert(matindex >= 0);
  EOS eosi = SGAPPLYMODSIMPLE(StiffGas(gm1, Cv, Pinf, qq));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
  EOS eos_ = SGAPPLYMOD(StiffGas(gm1, Cv, Pinf, qq));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}

int init_sg_StiffGas(const int matindex, EOS *eos, const double gm1, const double Cv,
                     const double Pinf, const double qq) {
  return init_sg_StiffGas(matindex, eos, gm1, Cv, Pinf, qq, def_en, def_v);
}

int init_sg_NobleAbel(const int matindex, EOS *eos, const double gm1, const double Cv,
                      const double bb, const double qq, int const *const enabled,
                      double *const vals) {
  assert(matindex >= 0);
  EOS eosi = SGAPPLYMODSIMPLE(NobleAbel(gm1, Cv, bb, qq));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
  EOS eos_ = SGAPPLYMOD(NobleAbel(gm1, Cv, bb, qq));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}

int init_sg_NobleAbel(const int matindex, EOS *eos, const double gm1, const double Cv,
                      const double bb, const double qq) {
  return init_sg_NobleAbel(matindex, eos, gm1, Cv, bb, qq, def_en, def_v);
}

#ifdef SINGULARITY_USE_SPINER_WITH_HDF5

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
int init_sg_eospac(const int matindex, EOS *eos, const int id, double *const eospac_vals,
                   int const *const enabled, double *const vals) {

  using namespace EospacWrapper;
  assert(matindex >= 0);
  bool invert_at_setup = eospac_vals[0];
  double insert_data = eospac_vals[1];
  eospacMonotonicity monotonicity = static_cast<eospacMonotonicity>(eospac_vals[2]);
  bool apply_smoothing = eospac_vals[3];
  eospacSplit apply_splitting = static_cast<eospacSplit>(eospac_vals[4]);
  bool linear_interp = eospac_vals[5];

  EOS eosi = SGAPPLYMODSIMPLE(
      EOSPAC(id, invert_at_setup = invert_at_setup, insert_data = insert_data,
             monotonicity = monotonicity, apply_smoothing = apply_smoothing,
             apply_splitting = apply_splitting, linear_interp = linear_interp));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
  EOS eos_ = SGAPPLYMOD(
      EOSPAC(id, invert_at_setup = invert_at_setup, insert_data = insert_data,
             monotonicity = monotonicity, apply_smoothing = apply_smoothing,
             apply_splitting = apply_splitting, linear_interp = linear_interp));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}
int init_sg_eospac(const int matindex, EOS *eos, const int id,
                   double *const eospac_vals) {
  return init_sg_eospac(matindex, eos, id, eospac_vals, def_en, def_v);
}
#endif // SINGULARITY_USE_EOSPAC

int init_sg_SAP_Polynomial(const int matindex, EOS *eos, const double rho0,
                           const double a0, const double a1, const double a2c,
                           const double a2e, const double a3, const double b0,
                           const double b1, const double b2c, const double b2e,
                           const double b3, int const *const enabled,
                           double *const vals) {
  assert(matindex >= 0);
  EOS eosi =
      SGAPPLYMODSIMPLE(SAP_Polynomial(rho0, a0, a1, a2c, a2e, a3, b0, b1, b2c, b2e, b3));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
  EOS eos_ = SGAPPLYMOD(SAP_Polynomial(rho0, a0, a1, a2c, a2e, a3, b0, b1, b2c, b2e, b3));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}

int init_sg_StiffGas(const int matindex, EOS *eos, const double gm1, const double Cv,
                     const double Pinf, const double qq, int const *const enabled,
                     double *const vals) {
  assert(matindex >= 0);
  EOS eosi = SGAPPLYMODSIMPLE(StiffGas(gm1, Cv, Pinf, qq));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
  EOS eos_ = SGAPPLYMOD(StiffGas(gm1, Cv, Pinf, qq));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}

int init_sg_StiffGas(const int matindex, EOS *eos, const double gm1, const double Cv,
                     const double Pinf, const double qq) {
  return init_sg_StiffGas(matindex, eos, gm1, Cv, Pinf, qq, def_en, def_v);
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
