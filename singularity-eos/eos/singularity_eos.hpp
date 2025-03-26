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

#ifndef _SINGULARITY_EOS_EOS_SINGULARITY_EOS_HPP_
#define _SINGULARITY_EOS_EOS_SINGULARITY_EOS_HPP_

#include <singularity-eos/eos/eos.hpp>

using singularity::EOS;

#if defined(__cplusplus)
extern "C" {
#endif

int init_sg_eos(const int nmat, EOS *&eos);

int init_sg_IdealGas(const int matindex, EOS *eos, const double gm1, const double Cv,
                     int const *const enabled, double *const vals);

int init_sg_JWL(const int matindex, EOS *eos, const double A, const double B,
                const double R1, const double R2, const double w, const double rho0,
                const double Cv, int const *const enabled, double *const vals);

int init_sg_Gruneisen(const int matindex, EOS *eos, const double C0, const double s1,
                      const double s2, const double s3, const double G0, const double b,
                      const double rho0, const double T0, const double P0,
                      const double Cv, int const *const enabled, double *const vals);

int init_sg_DavisProducts(const int matindex, EOS *eos, const double a, const double b,
                          const double k, const double n, const double vc,
                          const double pc, const double Cv, int const *const enabled,
                          double *const vals);

int init_sg_DavisReactants(const int matindex, EOS *eos, const double rho0,
                           const double e0, const double P0, const double T0,
                           const double A, const double B, const double C,
                           const double G0, const double Z, const double alpha,
                           const double Cv0, int const *const enabled,
                           double *const vals);

int init_sg_SAP_Polynomial(const int matindex, EOS *eos, const double rho0,
                           const double a0, const double a1, const double a2c,
                           const double a2e, const double a3, const double b0,
                           const double b1, const double b2c, const double b2e,
                           const double b3, int const *const enabled, double *const vals);

int init_sg_StiffGas(const int matindex, EOS *eos, const double gm1, const double Cv,
                     const double Pinf, const double qq, int const *const enabled,
                     double *const vals);

int init_sg_NobleAbel(const int matindex, EOS *eos, const double gm1, const double Cv,
                      const double bb, const double qq, int const *const enabled,
                      double *const vals);

int init_sg_CarnahanStarling(const int matindex, EOS *eos, const double gm1,
                             const double Cv, const double bb, const double qq,
                             int const *const enabled, double *const vals);

#ifdef SINGULARITY_USE_SPINER_WITH_HDF5

#ifdef SINGULARITY_USE_HELMHOLTZ
int init_sg_Helmholtz(const int matindex, EOS *eos, const char *filename, const bool rad,
                      const bool gas, const bool coul, const bool ion, const bool ele,
                      const bool verbose, int const *const enabled, double *const vals);
#endif

int init_sg_SpinerDependsRhoT(const int matindex, EOS *eos, const char *filename,
                              const int id, int const *const enabled, double *const vals);

int init_sg_SpinerDependsRhoSie(const int matindex, EOS *eos, const char *filename,
                                const int id, int const *const enabled,
                                double *const vals);
#endif

#ifdef SINGULARITY_USE_EOSPAC
// capitalize? eospaceos Eospac Eospaceos EOSPAC EOSPACeos?
int init_sg_eospac(const int matindex, EOS *eos, const int id, double *const eospac_vals,
                   int const *const enabled, double *const vals);
#endif // SINGULARITY_USE_EOSPAC

int get_sg_EntropyFromDensityInternalEnergy(int matindex, EOS *eos, const double *rhos,
                                             const double *sies, double *entropies,
                                             const int len, const int stride,
                                             double *lambda_data);

int get_sg_PressureFromDensityInternalEnergy(int matindex, EOS *eos, const double *rhos,
                                             const double *sies, double *pressures,
                                             const int len, const int stride,
                                             double *lambda_data);

int get_sg_MinInternalEnergyFromDensity(int matindex, EOS *eos, const double *rhos,
                                        double *sies, const int len);

int get_sg_BulkModulusFromDensityInternalEnergy(int matindex, EOS *eos,
                                                const double *rhos, const double *sies,
                                                double *bmods, const int len,
                                                const int stride, double *lambda_data);

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
    double mass_frac_cutoff);

int finalize_sg_eos(const int nmat, EOS *&eos, const int own_kokkos = 0);

#if defined(__cplusplus)
}
#endif

// outside of C scope, provide overloads
int init_sg_IdealGas(const int matindex, EOS *eos, const double gm1, const double Cv);

int init_sg_JWL(const int matindex, EOS *eos, const double A, const double B,
                const double R1, const double R2, const double w, const double rho0,
                const double Cv);

int init_sg_Gruneisen(const int matindex, EOS *eos, const double C0, const double s1,
                      const double s2, const double s3, const double G0, const double b,
                      const double rho0, const double T0, const double P0,
                      const double Cv);

int init_sg_DavisProducts(const int matindex, EOS *eos, const double a, const double b,
                          const double k, const double n, const double vc,
                          const double pc, const double Cv);

int init_sg_DavisReactants(const int matindex, EOS *eos, const double rho0,
                           const double e0, const double P0, const double T0,
                           const double A, const double B, const double C,
                           const double G0, const double Z, const double alpha,
                           const double Cv0);

int init_sg_StiffGas(const int matindex, EOS *eos, const double gm1, const double Cv,
                     const double Pinf, const double qq);

int init_sg_NobleAbel(const int matindex, EOS *eos, const double gm1, const double Cv,
                      const double bb, const double qq);

int init_sg_CarnahanStarling(const int matindex, EOS *eos, const double gm1,
                             const double Cv, const double bb, const double qq);

#ifdef SINGULARITY_USE_SPINER_WITH_HDF5

#ifdef SINGULARITY_USE_HELMHOLTZ
int init_sg_Helmholtz(const int matindex, EOS *eos, const char *filename, const bool rad,
                      const bool gas, const bool coul, const bool ion, const bool ele,
                      const bool verbose);
#endif

int init_sg_SpinerDependsRhoT(const int matindex, EOS *eos, const char *filename,
                              const int id);

int init_sg_SpinerDependsRhoSie(const int matindex, EOS *eos, const char *filename,
                                const int id);
#endif

#ifdef SINGULARITY_USE_EOSPAC
// capitalize? eospaceos Eospac Eospaceos EOSPAC EOSPACeos?
int init_sg_eospac(const int matindex, EOS *eos, const int id, const double *eospac_vals);
#endif // SINGULARITY_USE_EOSPAC

#endif // EOS_SINGULARITY_EOS_HPP_
