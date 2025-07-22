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

#include <cassert>
#include <ports-of-call/portability.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_builder.hpp>
#include <singularity-eos/eos/singularity_eos.hpp>
#include <singularity-eos/eos/singularity_eos_init_utils.hpp>

/*
===============================================
2D Lambda indexer class. Usage is as follows:

Assuming lambda is a std::array<Real, nCell>, instantiate this class
as

   idx = lambdaIndexer2D(lambda.data(), n)

We can now use the [] operator as follows:

   idx[i]

which will return the memory address of the element n*i of the array lambda
===============================================
*/

class lambdaIndexer2D {
 public:
  lambdaIndexer2D(int n, double *data) : n_(n), data_(data) {}

  PORTABLE_FORCEINLINE_FUNCTION
  double *operator[](int i) const { return &(data_[n_ * i]); }

 private:
  int n_;
  double *data_;
};

namespace singularity {
int def_en[4] = {0, 0, 0, 0};
double def_v[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
} // namespace singularity

using namespace singularity;

int init_sg_eos(const int nmat, EOS *&eos) {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  if (!Kokkos::is_initialized()) Kokkos::initialize();
#endif // PORTABILITY_STRATEGY_KOKKOS
  EOS *eos_p = new EOS[nmat];
  eos = eos_p;
  return 0;
}

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
                          const double pc, const double Cv, int const *const enabled,
                          double *const vals) {
  assert(matindex >= 0);
  EOS eosi = SGAPPLYMODSIMPLE(DavisProducts(a, b, k, n, vc, pc, Cv));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
  EOS eos_ = SGAPPLYMOD(DavisProducts(a, b, k, n, vc, pc, Cv));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}

int init_sg_DavisProducts(const int matindex, EOS *eos, const double a, const double b,
                          const double k, const double n, const double vc,
                          const double pc, const double Cv) {
  return init_sg_DavisProducts(matindex, eos, a, b, k, n, vc, pc, Cv, def_en, def_v);
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
#if SINGULARITY_USE_V_AND_V_EOS
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
#else
  PORTABLE_THROW_OR_ABORT("SAP Polynomial not currently supported. Please build with "
                          "-DSINGULARITY_USE_V_AND_V_EOS");
  return 1;
#endif // SINGULARITY_USE_V_AND_V_EOS
}

int init_sg_StiffGas(const int matindex, EOS *eos, const double gm1, const double Cv,
                     const double Pinf, const double qq, int const *const enabled,
                     double *const vals) {
#if SINGULARITY_USE_V_AND_V_EOS
  assert(matindex >= 0);
  EOS eosi = SGAPPLYMODSIMPLE(StiffGas(gm1, Cv, Pinf, qq));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
  EOS eos_ = SGAPPLYMOD(StiffGas(gm1, Cv, Pinf, qq));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
#else
  PORTABLE_THROW_OR_ABORT("Stiff Gas not currently supported. Please build with "
                          "-DSINGULARITY_USE_V_AND_V_EOS");
  return 1;
#endif // SINGULARITY_USE_V_AND_V_EOS
}

int init_sg_StiffGas(const int matindex, EOS *eos, const double gm1, const double Cv,
                     const double Pinf, const double qq) {
  return init_sg_StiffGas(matindex, eos, gm1, Cv, Pinf, qq, def_en, def_v);
}

int init_sg_NobleAbel(const int matindex, EOS *eos, const double gm1, const double Cv,
                      const double bb, const double qq, int const *const enabled,
                      double *const vals) {
#if SINGULARITY_USE_V_AND_V_EOS
  assert(matindex >= 0);
  EOS eosi = SGAPPLYMODSIMPLE(NobleAbel(gm1, Cv, bb, qq));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
  EOS eos_ = SGAPPLYMOD(NobleAbel(gm1, Cv, bb, qq));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
#else
  PORTABLE_THROW_OR_ABORT("NobleAbel not currently supported. Please build with "
                          "-DSINGULARITY_USE_V_AND_V_EOS");
  return 1;
#endif // SINGULARITY_USE_V_AND_V_EOS
}

int init_sg_NobleAbel(const int matindex, EOS *eos, const double gm1, const double Cv,
                      const double bb, const double qq) {
  return init_sg_NobleAbel(matindex, eos, gm1, Cv, bb, qq, def_en, def_v);
}

int init_sg_CarnahanStarling(const int matindex, EOS *eos, const double gm1,
                             const double Cv, const double bb, const double qq,
                             int const *const enabled, double *const vals) {
#if SINGULARITY_USE_V_AND_V_EOS
  assert(matindex >= 0);
  EOS eosi = SGAPPLYMODSIMPLE(CarnahanStarling(gm1, Cv, bb, qq));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
  EOS eos_ = SGAPPLYMOD(CarnahanStarling(gm1, Cv, bb, qq));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
#else
  PORTABLE_THROW_OR_ABORT("Carnahan Starling not currently supported. Please build with "
                          "-DSINGULARITY_USE_V_AND_V_EOS");
  return 1;
#endif // SINGULARITY_USE_V_AND_V_EOS
}

int init_sg_CarnahanStarling(const int matindex, EOS *eos, const double gm1,
                             const double Cv, const double bb, const double qq) {
  return init_sg_CarnahanStarling(matindex, eos, gm1, Cv, bb, qq, def_en, def_v);
}

#ifdef SINGULARITY_USE_SPINER_WITH_HDF5

#ifdef SINGULARITY_USE_HELMHOLTZ
int init_sg_Helmholtz(const int matindex, EOS *eos, const char *filename, const bool rad,
                      const bool gas, const bool coul, const bool ion, const bool ele,
                      const bool verbose, int const *const enabled, double *const vals) {
  assert(matindex >= 0);
  EOS eosi = SGAPPLYMODSIMPLE(
      Helmholtz(std::string(filename), rad, gas, coul, ion, ele, verbose));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
  EOS eos_ =
      SGAPPLYMOD(Helmholtz(std::string(filename), rad, gas, coul, ion, ele, verbose));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}
int init_sg_Helmholtz(const int matindex, EOS *eos, const char *filename, const bool rad,
                      const bool gas, const bool coul, const bool ion, const bool ele,
                      const bool verbose) {
  return init_sg_Helmholtz(matindex, eos, filename, rad, gas, coul, ion, ele, verbose,
                           def_en, def_v);
}
#endif

int init_sg_SpinerDependsRhoT(const int matindex, EOS *eos, const char *filename,
                              const int matid, const TableSplit split,
                              int const *const enabled, double *const vals) {
  assert(matindex >= 0);
  EOS eosi = SGAPPLYMODSIMPLE(SpinerEOSDependsRhoT(std::string(filename), matid, split));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
  EOS eos_ = SGAPPLYMOD(SpinerEOSDependsRhoT(std::string(filename), matid, split));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}

int init_sg_SpinerDependsRhoT(const int matindex, EOS *eos, const char *filename,
                              const int matid, const TableSplit split) {
  return init_sg_SpinerDependsRhoT(matindex, eos, filename, matid, split, def_en, def_v);
}

int init_sg_SpinerDependsRhoT(const int matindex, EOS *eos, const char *filename,
                              const int matid) {
  return init_sg_SpinerDependsRhoT(matindex, eos, filename, matid, TableSplit::Total);
}

int init_sg_SpinerDependsRhoSie(const int matindex, EOS *eos, const char *filename,
                                const int matid, const TableSplit split,
                                int const *const enabled, double *const vals) {
  assert(matindex >= 0);
  EOS eosi =
      SGAPPLYMODSIMPLE(SpinerEOSDependsRhoSie(std::string(filename), matid, split));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
  EOS eos_ = SGAPPLYMOD(SpinerEOSDependsRhoSie(std::string(filename), matid, split));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}

int init_sg_SpinerDependsRhoSie(const int matindex, EOS *eos, const char *filename,
                                const int matid, const TableSplit split) {
  return init_sg_SpinerDependsRhoSie(matindex, eos, filename, matid, split, def_en,
                                     def_v);
}

int init_sg_SpinerDependsRhoSie(const int matindex, EOS *eos, const char *filename,
                                const int matid) {
  return init_sg_SpinerDependsRhoSie(matindex, eos, filename, matid, TableSplit::Total);
}
#endif

#ifdef SINGULARITY_USE_EOSPAC
int init_sg_eospac(const int matindex, EOS *eos, const int id, const TableSplit split,
                   double *const eospac_vals, int const *const enabled,
                   double *const vals) {

  using namespace EospacWrapper;
  assert(matindex >= 0);
  bool invert_at_setup = eospac_vals[0];
  double insert_data = eospac_vals[1];
  eospacMonotonicity monotonicity = static_cast<eospacMonotonicity>(eospac_vals[2]);
  bool apply_smoothing = eospac_vals[3];
  eospacSplit apply_splitting = static_cast<eospacSplit>(eospac_vals[4]);
  bool linear_interp = eospac_vals[5];

  EOS eosi = SGAPPLYMODSIMPLE(
      EOSPAC(id, split, invert_at_setup = invert_at_setup, insert_data = insert_data,
             monotonicity = monotonicity, apply_smoothing = apply_smoothing,
             apply_splitting = apply_splitting, linear_interp = linear_interp));
  if (enabled[3] == 1) {
    singularity::pAlpha2BilinearRampParams(eosi, vals[2], vals[3], vals[4], vals[2],
                                           vals[3], vals[4], vals[5]);
  }
  EOS eos_ = SGAPPLYMOD(
      EOSPAC(id, split, invert_at_setup = invert_at_setup, insert_data = insert_data,
             monotonicity = monotonicity, apply_smoothing = apply_smoothing,
             apply_splitting = apply_splitting, linear_interp = linear_interp));
  eos[matindex] = eos_.GetOnDevice();
  return 0;
}

int init_sg_eospac(const int matindex, EOS *eos, const int id, const TableSplit split,
                   double *const eospac_vals) {
  return init_sg_eospac(matindex, eos, id, split, eospac_vals, def_en, def_v);
}
#endif // SINGULARITY_USE_EOSPAC

int get_sg_EntropyFromDensityInternalEnergy(int matindex, EOS *eos, const double *rhos,
                                            const double *sies, double *entropies,
                                            const int len, const int stride = -1,
                                            double *lambda_data = nullptr) {
  if (stride != -1 && lambda_data != nullptr) {
    lambdaIndexer2D idx(stride, lambda_data);
    eos[matindex].EntropyFromDensityInternalEnergy(rhos, sies, entropies, len, idx);
  } else
    eos[matindex].EntropyFromDensityInternalEnergy(rhos, sies, entropies, len);

  return 0;
}
int get_sg_PressureFromDensityInternalEnergy(int matindex, EOS *eos, const double *rhos,
                                             const double *sies, double *pressures,
                                             const int len, const int stride = -1,
                                             double *lambda_data = nullptr) {
  if (stride != -1 && lambda_data != nullptr) {
    lambdaIndexer2D idx(stride, lambda_data);
    eos[matindex].PressureFromDensityInternalEnergy(rhos, sies, pressures, len, idx);
  } else
    eos[matindex].PressureFromDensityInternalEnergy(rhos, sies, pressures, len);

  return 0;
}
int get_sg_MinInternalEnergyFromDensity(int matindex, EOS *eos, const double *rhos,
                                        double *sies, const int len) {
  eos[matindex].MinInternalEnergyFromDensity(rhos, sies, len);
  return 0;
}
int get_sg_BulkModulusFromDensityInternalEnergy(int matindex, EOS *eos,
                                                const double *rhos, const double *sies,
                                                double *bmods, const int len,
                                                const int stride = -1,
                                                double *lambda_data = nullptr) {
  if (stride != -1 && lambda_data != nullptr) {
    lambdaIndexer2D idx(stride, lambda_data);
    eos[matindex].BulkModulusFromDensityInternalEnergy(rhos, sies, bmods, len, idx);
  } else
    eos[matindex].BulkModulusFromDensityInternalEnergy(rhos, sies, bmods, len);

  return 0;
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
