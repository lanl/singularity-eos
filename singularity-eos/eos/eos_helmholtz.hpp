//------------------------------------------------------------------------------
// The Helmholtz EOS was first presented in: Timmes and Swesty, ApJS
// 126:501-516 (2000).  Reference implementation in fortan is here:
// https://cococubed.com/code_pages/eos.shtml
// Reference implementation is CC Frank TImmes and Douglas Swesty, 1999.
// Original work is open-sourced under the CC-By license
// https://creativecommons.org/licenses/by/4.0/
//------------------------------------------------------------------------------
// © 2023-2025. Triad National Security, LLC. All rights reserved.  This
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

#ifndef _SINGULARITY_EOS_EOS_HELMHOLTZ_HPP_
#define _SINGULARITY_EOS_EOS_HELMHOLTZ_HPP_

// TODO(JMM): Currently spiner and HDf5 are entangled. But really this
// model should depend only on spiner and NOT HDF5.
#ifdef SINGULARITY_USE_SPINER_WITH_HDF5
#ifdef SINGULARITY_USE_HELMHOLTZ

#include <cstdio>

#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// ports of call
#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>

// singularity-eos
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/hermite.hpp>
#include <singularity-eos/base/indexable_types.hpp>
#include <singularity-eos/base/math_utils.hpp>
#include <singularity-eos/base/robust_utils.hpp>
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>
#include <singularity-eos/base/spiner_table_utils.hpp>
#include <singularity-eos/eos/eos_base.hpp>

// spiner
#include <spiner/databox.hpp>

/*
 * The singularity-eos implementation of the Helmhotlz equation of state
 * provided by Timmes, and Swesty
 * Astrophysical Journal Supplements, Volume 126, Issue 2, Page 501 (2000).
 *
 * Port primarily by Jonah Miller with assistance from Philipp
 * Edelmann, Sam Jones, and the 2023 Co-Design Summer School.
 *
 * The Helmholtz EOS is a three part thermodynamically consistent EOS
 * for a hot, ionized gas. It consists of a thermal radiation term:
 *
 * P = sigma T^4
 *
 * an ions term, treated as an ideal gas:
 *
 * P = (gamma - 1) rho e
 *
 * and a degenerate electron term. Additionally, coulomb force
 * corrections can be applied on top of the full model.
 *
 * This multi-component model depends on the relative abundances of
 * electrons and ions, as well as the atomic mass and charge of the
 * ions. As such, the Helmholtz EOS requires two additional indepenent
 * variables, the average atomic mass, Abar, and the average atomic
 * number, Zbar. These are passed in through the lambda pointer. As
 * with the other tabulated EOS's, the log of the temperature is also
 * stored in the lambda pointer as a cache for root finding.
 *
 * From the density, Abar, and Zbar, several convenient derived
 * quantities can be comnputed, which are used internally. These are:
 *
 * - the ratio of electron number to Baryon number, called electron
 *   fraction, or Ye,
 * - The ion number fraction Ytot
 * - The electron number density De
 * - The log10 of De, lDe
 * - Additionally a convenience ration, ywot, is used to set atomic
 *   mass for the EOS if the gas is not ionized.
 *
 * The degenerate electron term is computed via thermodynamic
 * derivatives of the Helmholtz free energy (hence the name Helmholtz
 * EOS). The free energy is pre-computed via integrals over the Fermi
 * sphere and tabulated in a file, helmholtz/helm_table.dat, provided
 * from
 *
 * https://cococubed.com/code_pages/eos.shtml
 *
 * The table is a simple small ascii file. To ensure thermodyanic
 * consistency, the table is interpolated using either biquintic or
 * bicubic Hermite polynomials, which are sufficiently high order that
 * their high-order derivatives match the underlying data.
 *
 * The implication of interpolating from the free energy is that each
 * EOS evaluation provides ALL relevant EOS data and thermodynamic
 * derivatives. Thus the per-quantity EOS calls are relatively
 * inefficient, and it is instead better to use the FillEos call to
 * get the entire model at once.
 *
 * For this reason, each internal EOS function fills 5 arrays, each
 * containing relevant thermodynamic quantities and their derivatives
 * with respect to the independent variables, density, temperature,
 * abar, and zbar (in that order):
 * - p = pressure
 * - e = specific internal energy
 * - s = entropy
 * - etaele = checmical potential of electrons
 * - xne = number density of electron + positron pairs
 *
 * From a design standpoint, I have split each piece of the EOS off
 * and only combined them in the final "Helmholtz" class. The
 * component classes are HelmRad, HelmIons, HelmCoul, and
 * HelmElectrons. Right now, they are hardcoded, but I chose this
 * modular design with the intent, that these classes might be
 * hot-swapped with more sophisticated treatments down the line.
 *
 * Note that while singularity-eos has an ideal gas EOS, it's not
 * written in terms of fundamental quantities like atomic mass, which
 * made the translation cumbersom, so I just ported the reference
 * implementation.
 *
 * Note that e.g. the Gruneisenparameter is defined differently
 * compared to other EOSs. Here the Gruneisenparameter is the
 * Gamma3 of Cox & Giuli 1968 (Princiiples of Stellar Structure),
 * c&g in the following. I.e.
 * Gamma3 - 1 = (d ln T / d ln rho)|ad
 *
 * Some important formulas to be used when using this EOS:
 * - the temperature and density exponents (c&g 9.81 9.82)
 * - the specific heat at constant volume (c&g 9.92)
 * - the third adiabatic exponent (c&g 9.93)
 * - the first adiabatic exponent (c&g 9.97)
 * - the second adiabatic exponent (c&g 9.105)
 * - the specific heat at constant pressure (c&g 9.98)
 * - and relativistic formula for the sound speed (c&g 14.29)
 *
 */

namespace singularity {
using namespace eos_base;

// TODO(JMM): Maybe want to move these utility functions into something like an
// ASCII-utils file. Worth considering at some later date.
namespace HelmUtils {
using DataBox = Spiner::DataBox<>;

// Components of the arrays returned by internal routines:
// Variable, derivs w.r.t density, temperature, abar, zbar
constexpr std::size_t NDERIV = 5;
enum DERIV { VAL = 0, DDR = 1, DDT = 2, DDA = 3, DDZ = 4 };

// Tail-recursive resize tables
inline void ResizeTables(int n1, Real r1min, Real r1max, int n0, Real r0min, Real r0max,
                         DataBox &db) {
  db.resize(n1, n0); // set shape and log bounds
  db.setRange(1, r1min, r1max, n1);
  db.setRange(0, r0min, r0max, n0);
}
template <typename... Args>
inline void ResizeTables(int n1, Real r1min, Real r1max, int n0, Real r0min, Real r0max,
                         DataBox &head, Args &&...tail) {
  ResizeTables(n1, r1min, r1max, n0, r0min, r0max, head);
  ResizeTables(n1, r1min, r1max, n0, r0min, r0max, std::forward<Args>(tail)...);
}
// Tail-recursive read one i,j from ASCII text file
template <typename Msg_t, typename T>
inline void Read(std::ifstream &file, Msg_t &error_msg, T &var) {
  if (!(file >> var)) {
    file.close(); // is this needed?
    PORTABLE_ALWAYS_THROW_OR_ABORT(error_msg);
  }
}
template <typename Msg_t, typename T, typename... Args>
inline void Read(std::ifstream &file, Msg_t &error_msg, T &head, Args &&...tail) {
  Read(file, error_msg, head);
  Read(file, error_msg, std::forward<Args>(tail)...);
}
// Read all i,j from text file
template <typename... Args>
inline void SetTablesFromFile(std::ifstream &file, int n1, Real r1min, Real r1max, int n0,
                              Real r0min, Real r0max, Args &&...tables) {
  ResizeTables(n1, r1min, r1max, n0, r0min, r0max, std::forward<Args>(tables)...);
  for (int j = 0; j < n1; ++j) {
    for (int i = 0; i < n0; ++i) {
      std::stringstream error_msg;
      error_msg << "Error reading the Helmholtz free energy table at j = " << j
                << ", i = " << i << std::endl;
      Read(file, error_msg, tables(j, i)...);
    }
  }
}
} // namespace HelmUtils

class HelmElectrons {
  friend class table_utils::SpinerTricks<HelmElectrons>;
  using SpinerTricks = table_utils::SpinerTricks<HelmElectrons>;

 public:
  // may change with time
  using DataBox = HelmUtils::DataBox;
  static constexpr std::size_t NDERIV = HelmUtils::NDERIV;

  HelmElectrons() = default;
  inline HelmElectrons(const std::string &filename) {
    InitDataFile_(filename);
    CheckParams();
  }

  inline HelmElectrons GetOnDevice();
  inline void Finalize();
  inline void CheckParams() const {
    // better than nothing...
    PORTABLE_ALWAYS_REQUIRE(rho_.size() == NRHO, "Density grid correct");
    PORTABLE_ALWAYS_REQUIRE(T_.size() == NTEMP, "Temperature grid correct");
  }

  std::size_t DynamicMemorySizeInBytes() const {
    return SpinerTricks::DynamicMemorySizeInBytes(this);
  }
  std::size_t DumpDynamicMemory(char *dst) {
    return SpinerTricks::DumpDynamicMemory(dst, this);
  }
  std::size_t SetDynamicMemory(char *src,
                               const SharedMemSettings &stngs = DEFAULT_SHMEM_STNGS) {
    return SpinerTricks::SetDynamicMemory((stngs.data == nullptr) ? src : stngs.data,
                                          this);
  }

  PORTABLE_INLINE_FUNCTION
  void GetFromDensityTemperature(Real rho, Real lT, Real Ye, Real Ytot, Real De, Real lDe,
                                 Real pele[NDERIV], Real eele[NDERIV], Real sele[NDERIV],
                                 Real etaele[NDERIV], Real xne[NDERIV],
                                 bool only_e = false) const;

  // We COULD just expose the const vars under the hood, but I think
  // this is safer if, for example, the table size ever changes under
  // the hood.
  PORTABLE_FORCEINLINE_FUNCTION
  Real lTMin() const { return lTMin_; }
  PORTABLE_FORCEINLINE_FUNCTION
  Real lTMax() const { return lTMax_; }
  PORTABLE_FORCEINLINE_FUNCTION
  Real lRhoMin() const { return lRhoMin_; }
  PORTABLE_FORCEINLINE_FUNCTION
  Real lRhoMax() const { return lRhoMax_; }
  PORTABLE_FORCEINLINE_FUNCTION
  Real rhoMin() const { return rhoMin_; }
  PORTABLE_FORCEINLINE_FUNCTION
  Real rhoMax() const { return rhoMax_; }
  PORTABLE_FORCEINLINE_FUNCTION
  std::size_t numTemp() const { return NTEMP; }
  PORTABLE_FORCEINLINE_FUNCTION
  std::size_t numRho() const { return NRHO; }
  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumDensity() const { return rhoMin(); }
  PORTABLE_FORCEINLINE_FUNCTION
  Real MinimumTemperature() const { return TMin_; }
  PORTABLE_FORCEINLINE_FUNCTION
  Real MaximumDensity() const { return rhoMax(); }

 private:
  inline void InitDataFile_(const std::string &filename);

  // rho and T caches (to go between log/linear scale)
  DataBox rho_, T_;
  // Free energy and derivatives
  // Convention:
  // fd_ = df/drho
  // fdd_ = d^2 f/d rho^2
  // fdt_ = d^2 f/d rho dT
  // etc
  DataBox f_, fd_, ft_, fdd_, ftt_, fdt_, fddt_, fdtt_, fddtt_;
  // derivatives
  DataBox dpdf_, dpdfd_, dpdft_, dpdfdt_;
  // chemical potential
  DataBox ef_, efd_, eft_, efdt_;
  // number density
  DataBox xf_, xfd_, xft_, xfdt_;

#define DBLIST                                                                           \
  &rho_, &T_, &f_, &fd_, &ft_, &fdd_, &ftt_, &fdt_, &fddt_, &fdtt_, &fddtt_, &dpdf_,     \
      &dpdfd_, &dpdft_, &dpdfdt_, &ef_, &efd_, &eft_, &efdt_, &xf_, &xfd_, &xft_, &xfdt_
  auto GetDataBoxPointers_() const { return std::vector<const DataBox *>{DBLIST}; }
  auto GetDataBoxPointers_() { return std::vector<DataBox *>{DBLIST}; }
#undef DBLIST

  DataStatus memoryStatus_ = DataStatus::Deallocated;

  static constexpr std::size_t NTEMP = 101;
  static constexpr std::size_t NRHO = 271;

  // JMM: Following trick only works because the powers are even powers of 10.
  // This is to get compile-time values for log10/powers of 10
  // integer powers
  static constexpr int ilTMin_ = 3;
  static constexpr int ilTMax_ = 13;
  static constexpr int ilRhoMin_ = -12;
  static constexpr int ilRhoMax_ = 15;
  // real values to auto-cast correctly
  static constexpr Real lTMin_ = ilTMin_;
  static constexpr Real lTMax_ = ilTMax_;
  static constexpr Real lRhoMin_ = ilRhoMin_;
  static constexpr Real lRhoMax_ = ilRhoMax_;

  static constexpr Real dlT_ = (lTMax_ - lTMin_) / (static_cast<Real>(NTEMP) - 1.0);
  static constexpr Real dlRho_ = (lRhoMax_ - lRhoMin_) / (static_cast<Real>(NRHO) - 1.0);

  static constexpr Real TMin_ = math_utils::ipow10<ilTMin_>();
  static constexpr Real TMax_ = math_utils::ipow10<ilTMax_>();
  static constexpr Real rhoMin_ = math_utils::ipow10<ilRhoMin_>();
  static constexpr Real rhoMax_ = math_utils::ipow10<ilRhoMax_>();
};

// Currently trivial. But we may want to swap this out at some point,
// which is why I made it a class. We may also want to write a full
// radiation EOS at some point.
class HelmRad {
 public:
  static constexpr Real STEFAN_BOLTZMANN_CONSTANT = 5.67040047374e-5; /* g / K^4 s^3 */
  static constexpr Real CL = 2.99792458e10; // speed of light. cm/s
  static constexpr Real SIG_O_C = STEFAN_BOLTZMANN_CONSTANT / CL;
  static constexpr Real PREFACTOR = (4. / 3.) * SIG_O_C;
  static constexpr std::size_t NDERIV = HelmUtils::NDERIV;

  HelmRad() = default;
  HelmRad GetOnDevice() { return *this; }
  void Finalize() {}

  PORTABLE_INLINE_FUNCTION
  void GetFromDensityTemperature(const Real rho, const Real temp, Real prad[NDERIV],
                                 Real erad[NDERIV], Real srad[NDERIV]) const;
};

// TODO(JMM): Use singularity's built in ideal gas EOS
// rather than this hardcoded one
class HelmIon {
 public:
  static constexpr std::size_t NDERIV = HelmUtils::NDERIV;
  static constexpr Real KB = 1.3806504e-16; // Boltzmann constant in cgs
  static constexpr Real NA = 6.02214199e23; // Avogadro's number. 1/mol
  static constexpr Real KBNA = KB * NA;
  static constexpr Real KBi = 1.0 / KB;
  static constexpr Real UNIFIED_ATOMIC_MASS = 1.660538782e-24; /* g */
  static constexpr Real PLANCK_H = 6.62606896e-27;             /* g cm^2 / s */
  // 1.5 * ln((2.0 * M_PI * UNIFIED_ATOMIC_MASS * KB) / (PLANCK_H * PLANCK_H));
  static constexpr Real LSWOT15 = 46.682612633059801;

  HelmIon() = default;
  HelmIon GetOnDevice() { return *this; }
  void Finalize() {}

  PORTABLE_INLINE_FUNCTION
  void GetFromDensityTemperature(const Real rho, const Real temp, const Real y,
                                 const Real ytot, const Real xni, const Real dxnidd,
                                 const Real dxnida, Real pion[NDERIV], Real eion[NDERIV],
                                 Real sion[NDERIV]) const;
};

// Coulomb corrections. Again extraneous to make it a class. But
// perhaps that will make it easier to swap out down the line.
class HelmCoulomb {
 public:
  static constexpr Real LN10 = 2.30258509299405e+00; // ln(10)
  static constexpr Real ELECTRON_CHARGE_ESU = 4.80320427e-10;
  static constexpr Real KB = 1.3806504e-16; // Boltzmann constant in cgs
  static constexpr Real NA = 6.02214199e23; // Avogadro's number. 1/mol
  static constexpr Real KBNA = KB * NA;
  static constexpr std::size_t NDERIV = HelmUtils::NDERIV;

  HelmCoulomb() = default;
  HelmCoulomb GetOnDevice() { return *this; }
  void Finalize() {}

  PORTABLE_INLINE_FUNCTION
  void GetFromDensityTemperature(const Real rho, const Real temp, const Real ytot,
                                 const Real abar, const Real zbar, const Real xni,
                                 const Real dxnidd, const Real dxnida, Real pcoul[NDERIV],
                                 Real ecoul[NDERIV], Real scoul[NDERIV]) const;

 private:
  template <int n>
  PORTABLE_INLINE_FUNCTION void butterworth(Real freq, double cfreq, Real &g,
                                            Real &dgdf) const {
    g = robust::ratio(1.0, 1.0 + math_utils::pow<2 * n>(robust::ratio(freq, cfreq)));
    dgdf = robust::ratio(-math_utils::pow<2>(g) * 2 * n, math_utils::pow<2 * n>(cfreq)) *
           math_utils::pow<2 * n - 1>(freq);
  }
};

// The reason to separate out the full Helmholtz EOS and the
// electrons bit is that this should make things easier to
// mix/match/compose in the future. For example if more advanced ion
// corrections are desired.
class Helmholtz : public EosBase<Helmholtz> {
 public:
  static constexpr std::size_t NDERIV = HelmUtils::NDERIV;
  // These are the indexes in the lambda array
  struct Lambda {
    enum Index {
      Abar = 0, // Average atomic mass
      Zbar = 1, // Average atomic number
      lT = 2    // log10 temperature. used for root finding.
    };
  };
  // Options struct. You can create one of these and modify it to set
  // options at initialization
  struct Options {
    Options() = default;
    Options(const bool rad, const bool gas, const bool coul, const bool ion,
            const bool electron)
        : ENABLE_RAD(rad), ENABLE_GAS(gas), ENABLE_COULOMB_CORRECTIONS(coul),
          GAS_IONIZED(ion), GAS_DEGENERATE(electron) {}
    Options(const bool rad, const bool gas, const bool coul, const bool ion,
            const bool electron, const bool verbose)
        : ENABLE_RAD(rad), ENABLE_GAS(gas), ENABLE_COULOMB_CORRECTIONS(coul),
          GAS_IONIZED(ion), GAS_DEGENERATE(electron), VERBOSE(verbose) {}
    Options(const bool rad, const bool gas, const bool coul, const bool ion,
            const bool electron, const bool verbose, const bool newton_raphson)
        : ENABLE_RAD(rad), ENABLE_GAS(gas), ENABLE_COULOMB_CORRECTIONS(coul),
          GAS_IONIZED(ion), GAS_DEGENERATE(electron), VERBOSE(verbose),
          USE_NEWTON_RAPHSON(newton_raphson) {}
    bool ENABLE_RAD = true;
    bool ENABLE_GAS = true;
    bool ENABLE_COULOMB_CORRECTIONS = true;
    bool GAS_IONIZED = true;
    bool GAS_DEGENERATE = true;
    bool VERBOSE = false;
    bool USE_NEWTON_RAPHSON = true;
  };

  Helmholtz() = default;
  Helmholtz(const std::string &filename) : electrons_(filename) {}
  Helmholtz(const std::string &filename, Options options)
      : electrons_(filename), options_(options) {}
  Helmholtz(const std::string &filename, const bool rad, const bool gas, const bool coul,
            const bool ion, const bool ele)
      : electrons_(filename), options_(rad, gas, coul, ion, ele) {}
  Helmholtz(const std::string &filename, const bool rad, const bool gas, const bool coul,
            const bool ion, const bool ele, const bool verbose)
      : electrons_(filename), options_(rad, gas, coul, ion, ele, verbose) {}
  Helmholtz(const std::string &filename, const bool rad, const bool gas, const bool coul,
            const bool ion, const bool ele, const bool verbose, const bool newton_raphson)
      : electrons_(filename),
        options_(rad, gas, coul, ion, ele, verbose, newton_raphson) {}

  PORTABLE_INLINE_FUNCTION void CheckParams() const { electrons_.CheckParams(); }

  PORTABLE_INLINE_FUNCTION int nlambda() const noexcept { return 3; }
  template <typename T>
  static inline constexpr bool NeedsLambda() {
    return std::is_same<T, IndexableTypes::MeanAtomicMass>::value ||
           std::is_same<T, IndexableTypes::MeanAtomicNumber>::value ||
           std::is_same<T, IndexableTypes::LogTemperature>::value;
  }
  static constexpr unsigned long PreferredInput() {
    return thermalqs::density | thermalqs::temperature;
  }

  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    printf("Helmholtz Parameters:\n"
           "ENABLE_RAD = %d\n"
           "ENABLE_GAS = %d\n"
           "ENABLE_COULOMB_CORRECTIONS = %d\n"
           "GAS_IONIZED = %d\n"
           "GAS_DEGENERATE = %d\n",
           options_.ENABLE_RAD, options_.ENABLE_GAS, options_.ENABLE_COULOMB_CORRECTIONS,
           options_.GAS_IONIZED, options_.GAS_DEGENERATE);
  }

  inline Helmholtz GetOnDevice() {
    Helmholtz other;
    other.rad_ = rad_.GetOnDevice();
    other.ions_ = ions_.GetOnDevice();
    other.coul_ = coul_.GetOnDevice();
    other.electrons_ = electrons_.GetOnDevice();
    other.options_ = options_;
    return other;
  }
  inline void Finalize() {
    rad_.Finalize();
    ions_.Finalize();
    coul_.Finalize();
    electrons_.Finalize();
  }
  std::size_t DynamicMemorySizeInBytes() const {
    return electrons_.DynamicMemorySizeInBytes();
  }
  std::size_t DumpDynamicMemory(char *dst) { return electrons_.DumpDynamicMemory(dst); }
  std::size_t SetDynamicMemory(char *src,
                               const SharedMemSettings &stngs = DEFAULT_SHMEM_STNGS) {
    return electrons_.SetDynamicMemory(src);
  }

  PORTABLE_INLINE_FUNCTION
  void GetMassFractions(const Real rho, const Real temp, const Real ytot, Real &xni,
                        Real &dxnidd, Real &dxnida) const {
    xni = ions_.NA * ytot * rho;
    dxnidd = ions_.NA * ytot;
    dxnida = -xni * ytot;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    Real rl = rho;
    Real el = sie;
    Real temperature, p, cv, bmod;
    FillEos(rl, temperature, el, p, cv, bmod, thermalqs::temperature, lambda);
    return temperature;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    Real rl = rho;
    Real tl = temperature;
    Real sie, p, cv, bmod;
    FillEos(rl, tl, sie, p, cv, bmod, thermalqs::specific_internal_energy, lambda);
    return sie;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    Real rl = rho;
    Real tl = temperature;
    Real sie, p, cv, bmod;
    FillEos(rl, tl, sie, p, cv, bmod, thermalqs::pressure, lambda);
    return p;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    Real rl = rho;
    Real el = sie;
    Real temperature, p, cv, bmod;
    FillEos(rl, temperature, el, p, cv, bmod,
            thermalqs::pressure | thermalqs::temperature, lambda);
    return p;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real
  EntropyFromDensityTemperature(const Real rho, const Real temperature,
                                Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    Real p[NDERIV], e[NDERIV], s[NDERIV], etaele[NDERIV], nep[NDERIV];
    GetFromDensityTemperature_(rho, temperature, lambda, p, e, s, etaele, nep);
    return s[HelmUtils::VAL];
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    Real p[NDERIV], e[NDERIV], s[NDERIV], etaele[NDERIV], nep[NDERIV];
    GetFromDensityInternalEnergy_(rho, sie, lambda, p, e, s, etaele, nep);
    return s[HelmUtils::VAL];
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    Real rl = rho;
    Real tl = temperature;
    Real sie, p, cv, bmod;
    FillEos(rl, tl, sie, p, cv, bmod, thermalqs::specific_heat, lambda);
    return cv;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    Real rl = rho;
    Real el = sie;
    Real temperature, p, cv, bmod;
    FillEos(rl, temperature, el, p, cv, bmod,
            thermalqs::specific_heat | thermalqs::temperature, lambda);
    return cv;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    Real rl = rho;
    Real tl = temperature;
    Real sie, p, cv, bmod;
    FillEos(rl, tl, sie, p, cv, bmod, thermalqs::bulk_modulus, lambda);
    return bmod;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    Real rl = rho;
    Real el = sie;
    Real temperature, p, cv, bmod;
    FillEos(rl, temperature, el, p, cv, bmod,
            thermalqs::bulk_modulus | thermalqs::temperature, lambda);
    return bmod;
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      const Real rho, const Real temperature,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    using namespace HelmUtils;
    Real p[NDERIV], e[NDERIV], s[NDERIV], etaele[NDERIV], nep[NDERIV];
    GetFromDensityTemperature_(rho, temperature, lambda, p, e, s, etaele, nep);
    Real gamma3 = ComputeGamma3_(rho, temperature, p, e);
    return gamma3 - 1.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      const Real rho, const Real sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    using namespace HelmUtils;
    Real p[NDERIV], e[NDERIV], s[NDERIV], etaele[NDERIV], nep[NDERIV];
    Real abar = IndexerUtils::Get<IndexableTypes::MeanAtomicMass>(lambda, Lambda::Abar);
    Real zbar = IndexerUtils::Get<IndexableTypes::MeanAtomicNumber>(lambda, Lambda::Zbar);
    Real ytot, ye, ywot, De, lDe;
    GetElectronDensities_(rho, abar, zbar, ytot, ye, ywot, De, lDe);
    Real lT = lTFromRhoSie_(rho, sie, abar, zbar, ye, ytot, ywot, De, lDe, lambda);
    Real T = math_utils::pow10(lT);
    GetFromDensityLogTemperature_(rho, T, abar, zbar, ye, ytot, ywot, De, lDe, p, e, s,
                                  etaele, nep);
    Real gamma3 = ComputeGamma3_(rho, T, p, e);
    return gamma3 - 1.0;
  }

  PORTABLE_INLINE_FUNCTION
  Real MeanAtomicMass() const {
    PORTABLE_THROW_OR_ABORT("For Helmholtz EOS, mean atomic mass is an input!\n");
    return 1.0;
  }
  PORTABLE_INLINE_FUNCTION
  Real MeanAtomicNumber() const {
    PORTABLE_THROW_OR_ABORT("For Helmholtz EOS, mean atomic number is an input!\n");
    return 1.0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicMassFromDensityTemperature(
      const Real rho, const Real T,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    using namespace HelmUtils;
    return IndexerUtils::Get<IndexableTypes::MeanAtomicMass>(lambda, Lambda::Abar);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MeanAtomicNumberFromDensityTemperature(
      const Real rho, const Real T,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    using namespace HelmUtils;
    return IndexerUtils::Get<IndexableTypes::MeanAtomicNumber>(lambda, Lambda::Zbar);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
          const unsigned long output,
          Indexer_t &&lambda = static_cast<Real *>(nullptr)) const;

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                         Real &bmod, Real &dpde, Real &dvdt,
                         Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    // JMM: Conditions for an oxygen burning shell in a stellar
    // core. Not sure if that's the best choice.
    rho = 1e7;
    temp = 1.5e9;
    FillEos(rho, temp, sie, press, cv, bmod,
            thermalqs::specific_internal_energy | thermalqs::pressure |
                thermalqs::specific_heat | thermalqs::bulk_modulus,
            lambda);
    dpde = GruneisenParamFromDensityTemperature(rho, temp, lambda) * rho;
    dvdt = 0;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(const Real press, const Real temp,
                                       Indexer_t &&lambda, Real &rho, Real &sie) const {
    using RootFinding1D::regula_falsi;
    using RootFinding1D::Status;
    PORTABLE_REQUIRE(temp > = 0, "Non-negative temperature required");
    auto PofRT = [&](const Real r) {
      return PressureFromDensityTemperature(r, temp, lambda);
    };
    auto status = regula_falsi(PofRT, press, 1e7, electrons_.rhoMin(),
                               electrons_.rhoMax(), 1.0e-8, 1.0e-8, rho);
    if (status != Status::SUCCESS) {
      PORTABLE_THROW_OR_ABORT("Helmholtz::DensityEnergyFromPressureTemperature: "
                              "Root find failed to find a solution given P, T\n");
    }
    if (rho < 0.) {
      PORTABLE_THROW_OR_ABORT("Helmholtz::DensityEnergyFromPressureTemperature: "
                              "Root find produced negative energy solution\n");
    }
    sie = InternalEnergyFromDensityTemperature(rho, temp, lambda);
  }
  static std::string EosType() { return std::string("Helmholtz"); }
  static std::string EosPyType() { return EosType(); }

  SG_ADD_BASE_CLASS_USINGS(Helmholtz)
 private:
  PORTABLE_INLINE_FUNCTION
  Real ComputeGamma1_(const Real rho, const Real T, const Real p[NDERIV],
                      const Real e[NDERIV]) const {
    using namespace HelmUtils;
    // Gamma1 (c&g 9.97)
    const Real chit = robust::ratio(T, p[0]) * p[DDT];
    const Real chid = robust::ratio(p[DDR] * rho, p[0]);
    const Real x = robust::ratio(p[0] * chit, rho * T * e[DDT]);
    return chit * x + chid;
  }

  PORTABLE_INLINE_FUNCTION
  Real ComputeGamma3_(const Real rho, const Real T, const Real p[NDERIV],
                      const Real e[NDERIV]) const {
    using namespace HelmUtils;
    // Gamma3 (c&g 9.93)
    const Real chit = robust::ratio(T, p[0]) * p[DDT];
    const Real x = robust::ratio(p[0] * chit, rho * T * e[DDT]);
    return x + 1.0;
  }

  PORTABLE_INLINE_FUNCTION
  void GetElectronDensities_(const Real rho, const Real abar, const Real zbar, Real &ytot,
                             Real &ye, Real &ywot, Real &De, Real &lDe) const {
    ytot = robust::ratio(1.0, abar);
    ye = zbar * ytot;
    // TODO(JMM): should we be passing around logrho maybe?
    ywot = std::log(robust::ratio(abar * abar * std::sqrt(abar), rho * ions_.NA));
    De = rho * ye;
    lDe = std::log10(De);
  }

  PORTABLE_INLINE_FUNCTION
  Real lTAnalytic_(const Real rho, const Real e, const Real ni, const Real ne) const {
    return std::log10((2.0 / 3.0) * robust::ratio(e * rho, ni + ne) * ions_.KBi);
  }

  template <typename Indexer_t>
  PORTABLE_INLINE_FUNCTION void
  GetFromDensityTemperature_(const Real rho, const Real temperature, Indexer_t &&lambda,
                             Real p[NDERIV], Real e[NDERIV], Real s[NDERIV],
                             Real etaele[NDERIV], Real nep[NDERIV]) const {
    Real abar = IndexerUtils::Get<IndexableTypes::MeanAtomicMass>(lambda, Lambda::Abar);
    Real zbar = IndexerUtils::Get<IndexableTypes::MeanAtomicNumber>(lambda, Lambda::Zbar);
    Real lT = std::log10(temperature);
    IndexerUtils::Get<IndexableTypes::LogTemperature>(lambda, Lambda::lT) = lT;
    Real ytot, ye, ywot, De, lDe;
    GetElectronDensities_(rho, abar, zbar, ytot, ye, ywot, De, lDe);
    GetFromDensityLogTemperature_(rho, temperature, abar, zbar, ye, ytot, ywot, De, lDe,
                                  p, e, s, etaele, nep);
  }

  template <typename Indexer_t>
  PORTABLE_INLINE_FUNCTION void
  GetFromDensityInternalEnergy_(const Real rho, const Real sie, Indexer_t &&lambda,
                                Real p[NDERIV], Real e[NDERIV], Real s[NDERIV],
                                Real etaele[NDERIV], Real nep[NDERIV]) const {
    Real abar = IndexerUtils::Get<IndexableTypes::MeanAtomicMass>(lambda, Lambda::Abar);
    Real zbar = IndexerUtils::Get<IndexableTypes::MeanAtomicNumber>(lambda, Lambda::Zbar);
    Real ytot, ye, ywot, De, lDe;
    GetElectronDensities_(rho, abar, zbar, ytot, ye, ywot, De, lDe);
    Real lT = lTFromRhoSie_(rho, sie, abar, zbar, ye, ytot, ywot, De, lDe, lambda);
    Real T = math_utils::pow10(lT);
    GetFromDensityLogTemperature_(rho, T, abar, zbar, ye, ytot, ywot, De, lDe, p, e, s,
                                  etaele, nep);
  }

  PORTABLE_INLINE_FUNCTION
  void GetFromDensityLogTemperature_(
      const Real rho, const Real T, const Real abar, const Real zbar, const Real ye,
      const Real ytot, const Real ywot, const Real De, const Real lDe, Real p[NDERIV],
      Real e[NDERIV], Real s[NDERIV], Real etaele[NDERIV], Real nep[NDERIV],
      // TODO(JMM): Decide which of the quantities below to keep
      const bool only_e = false) const;

  template <typename Indexer_t>
  PORTABLE_INLINE_FUNCTION Real lTFromRhoSie_(const Real rho, const Real e,
                                              const Real abar, const Real zbar,
                                              const Real ye, const Real ytot,
                                              const Real ywot, const Real De,
                                              const Real lDe, Indexer_t &&lambda) const;

  static constexpr Real ROOT_THRESH = 1e-14;
  static constexpr Real HELM_EOS_EPS = 1e-10;
  Options options_;
  HelmRad rad_;
  HelmIon ions_;
  HelmCoulomb coul_;
  HelmElectrons electrons_;
};

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION void
Helmholtz::FillEos(Real &rho, Real &temp, Real &energy, Real &press, Real &cv, Real &bmod,
                   const unsigned long output, Indexer_t &&lambda) const {
  using namespace HelmUtils;
  bool need_temp = (output & thermalqs::temperature);
  bool need_sie = (output & thermalqs::specific_internal_energy);
  bool need_rho = (output & thermalqs::density);
  PORTABLE_ALWAYS_REQUIRE(!need_rho, "Density output not supported by this EOS");
  PORTABLE_ALWAYS_REQUIRE(
      !(need_temp && need_sie),
      "Either specific internal energy or temperature must be provided.");
  Real abar = IndexerUtils::Get<IndexableTypes::MeanAtomicMass>(lambda, Lambda::Abar);
  Real zbar = IndexerUtils::Get<IndexableTypes::MeanAtomicNumber>(lambda, Lambda::Zbar);
  Real ytot, ye, ywot, De, lDe, lT;
  GetElectronDensities_(rho, abar, zbar, ytot, ye, ywot, De, lDe);
  if (need_temp) {
    lT = lTFromRhoSie_(rho, energy, abar, zbar, ye, ytot, ywot, De, lDe, lambda);
    temp = math_utils::pow10(lT);
  } else {
    lT = std::log10(temp);
    IndexerUtils::Get<IndexableTypes::LogTemperature>(lambda, Lambda::lT) = lT;
  }
  Real p[NDERIV], e[NDERIV], s[NDERIV], etaele[NDERIV], nep[NDERIV];
  GetFromDensityLogTemperature_(rho, temp, abar, zbar, ye, ytot, ywot, De, lDe, p, e, s,
                                etaele, nep);
  if (output & thermalqs::specific_internal_energy) {
    energy = e[0];
  }
  if (output & thermalqs::pressure) {
    press = p[0];
  }
  if (output & thermalqs::specific_heat) {
    cv = e[DDT];
  }
  if (output & thermalqs::bulk_modulus) {
    Real gamma1 = ComputeGamma1_(rho, temp, p, e);

    // Note this is the bulk modulus, not the sound speed. To compute
    // the sound speed, you must take this value and divide by the
    // enthalpy. In particular:
    // c_s = sqrt(bmod / (w)) = sqrt(bmod / (h rho))
    // see page 108 of Rezzolla and Zanotti
    bmod = std::max(robust::EPS(), p[0] * gamma1);
  }
}

template <typename Indexer_t>
PORTABLE_INLINE_FUNCTION Real Helmholtz::lTFromRhoSie_(const Real rho, const Real e,
                                                       const Real abar, const Real zbar,
                                                       const Real ye, const Real ytot,
                                                       const Real ywot, const Real De,
                                                       const Real lDe,
                                                       Indexer_t &&lambda) const {
  using namespace HelmUtils;
  const Real abari = robust::ratio(1.0, abar);
  const Real ni = abari * rho * ions_.NA;
  const Real ne = zbar * ni;
  Real lT;
  Real T;

  if (options_.ENABLE_RAD || options_.GAS_DEGENERATE ||
      options_.ENABLE_COULOMB_CORRECTIONS) {
    Real lTguess = IndexerUtils::Get<IndexableTypes::LogTemperature>(lambda, Lambda::lT);
    if (!((electrons_.lTMin() <= lTguess) && (lTguess <= electrons_.lTMax()))) {
      lTguess = lTAnalytic_(rho, e, ni, ne);
      if (!((electrons_.lTMin() <= lTguess) && (lTguess <= electrons_.lTMax()))) {
        // Fallback initial guess
        lTguess = 8;
      }
    }
    Real Tguess = math_utils::pow10(lTguess);
    auto &copy = *this; // stupid C++17 workaround
    if (options_.USE_NEWTON_RAPHSON) {
      auto status = RootFinding1D::newton_raphson(
          [&](Real T) {
            Real p[NDERIV], e[NDERIV], s[NDERIV], etaele[NDERIV], nep[NDERIV];
            copy.GetFromDensityLogTemperature_(rho, T, abar, zbar, ye, ytot, ywot, De,
                                               lDe, p, e, s, etaele, nep, true);
            return std::make_tuple(e[VAL], e[DDT]);
          },
          e, Tguess, math_utils::pow10(electrons_.lTMin()),
          math_utils::pow10(electrons_.lTMax()), HELM_EOS_EPS, T, nullptr,
          options_.VERBOSE, false);
      if (status != RootFinding1D::Status::SUCCESS) {
        if (options_.VERBOSE) {
          printf("Newton-Raphson failed to converge, falling back to regula falsi\n");
        }
        status = RootFinding1D::regula_falsi(
            [&](Real T) {
              Real p[NDERIV], e[NDERIV], s[NDERIV], etaele[NDERIV], nep[NDERIV];
              copy.GetFromDensityLogTemperature_(rho, T, abar, zbar, ye, ytot, ywot, De,
                                                 lDe, p, e, s, etaele, nep, true);
              return e[VAL];
            },
            e, Tguess, math_utils::pow10(electrons_.lTMin()),
            math_utils::pow10(electrons_.lTMax()), ROOT_THRESH, ROOT_THRESH, T, nullptr,
            options_.VERBOSE);
        if (status != RootFinding1D::Status::SUCCESS) {
          lT = lTAnalytic_(rho, e, ni, options_.GAS_IONIZED * ne);
          T = math_utils::pow10(lT);
        }
      }
    } else {
      auto status = RootFinding1D::regula_falsi(
          [&](Real T) {
            Real p[NDERIV], e[NDERIV], s[NDERIV], etaele[NDERIV], nep[NDERIV];
            copy.GetFromDensityLogTemperature_(rho, T, abar, zbar, ye, ytot, ywot, De,
                                               lDe, p, e, s, etaele, nep, true);
            return e[VAL];
          },
          e, Tguess, math_utils::pow10(electrons_.lTMin()),
          math_utils::pow10(electrons_.lTMax()), ROOT_THRESH, ROOT_THRESH, T, nullptr,
          options_.VERBOSE);
      if (status != RootFinding1D::Status::SUCCESS) {
        lT = lTAnalytic_(rho, e, ni, options_.GAS_IONIZED * ne);
        T = math_utils::pow10(lT);
      }
    }
  } else {
    lT = lTAnalytic_(rho, e, ni, options_.GAS_IONIZED * ne);
    T = math_utils::pow10(lT);
  }
  lT = std::log10(T);
  // Make sure the result is within the table bounds
  if (lT < electrons_.lTMin()) {
    if (options_.VERBOSE) {
      printf("Temperature below table limit, setting lT = lTMin. (lT = %f)\n", lT);
    }
    lT = electrons_.lTMin();
  }
  if (lT > electrons_.lTMax()) {
    if (options_.VERBOSE) {
      printf("Temperature above table limit, setting lT = lTMax. (lT = %f)\n", lT);
    }
    lT = electrons_.lTMax();
  }
  IndexerUtils::Get<IndexableTypes::LogTemperature>(lambda, Lambda::lT) = lT;
  return lT;
}

inline void HelmElectrons::InitDataFile_(const std::string &filename) {
  using namespace HelmUtils;
  std::ifstream file(filename);
  PORTABLE_ALWAYS_REQUIRE(file.is_open(),
                          "HelmElectrons file " + filename + " not found!");
  SetTablesFromFile(file, NTEMP, lTMin_, lTMax_, NRHO, lRhoMin_, lRhoMax_, f_, fd_, ft_,
                    fdd_, ftt_, fdt_, fddt_, fdtt_, fddtt_);
  SetTablesFromFile(file, NTEMP, lTMin_, lTMax_, NRHO, lRhoMin_, lRhoMax_, dpdf_, dpdfd_,
                    dpdft_, dpdfdt_);
  SetTablesFromFile(file, NTEMP, lTMin_, lTMax_, NRHO, lRhoMin_, lRhoMax_, ef_, efd_,
                    eft_, efdt_);
  SetTablesFromFile(file, NTEMP, lTMin_, lTMax_, NRHO, lRhoMin_, lRhoMax_, xf_, xfd_,
                    xft_, xfdt_);
  file.close();

  auto &lRhoRange = f_.range(0);
  auto &lTRange = f_.range(1);
  rho_.resize(NRHO);
  for (int i = 0; i < NRHO; ++i) {
    rho_(i) = math_utils::pow10(lRhoRange.x(i));
  }
  T_.resize(NTEMP);
  for (int i = 0; i < NTEMP; ++i) {
    T_(i) = math_utils::pow10(lTRange.x(i));
  }

  memoryStatus_ = DataStatus::OnHost;
}

inline HelmElectrons HelmElectrons::GetOnDevice() {
  return SpinerTricks::GetOnDevice(this);
}

inline void HelmElectrons::Finalize() { SpinerTricks::Finalize(this); }

PORTABLE_INLINE_FUNCTION
void HelmElectrons::GetFromDensityTemperature(Real rho, Real lT, Real Ye, Real Ytot,
                                              Real De, Real lDe, Real pele[NDERIV],
                                              Real eele[NDERIV], Real sele[NDERIV],
                                              Real etaele[NDERIV], Real xne[NDERIV],
                                              bool only_e) const {
  // Bound lRho, lT
  rho = std::min(rhoMax(), std::max(rhoMin(), rho));
  De = std::min(rhoMax(), std::max(rhoMin(), De));
  lDe = std::min(lRhoMax(), std::max(lRhoMin(), lDe));
  Real T = math_utils::pow10(lT);

  // Find central indexes in table
  auto lRhoRange = f_.range(0);
  auto lTRange = f_.range(1);
  int iat = lRhoRange.index(lDe); // auto bounded
  int jat = lTRange.index(lT);

  // compute deltas
  Real dth = T_(jat + 1) - T_(jat);
  Real dt2 = dth * dth;
  Real dti = robust::ratio(1.0, dth);
  Real dt2i = dti * dti; // division is slow
  Real dd = rho_(iat + 1) - rho_(iat);
  Real dd2 = dd * dd;
  Real ddi = robust::ratio(1.0, dd);

  // contiguous cache of values for helm interp
  Real fi[36];
  fi[0] = f_(jat, iat);
  fi[1] = f_(jat, iat + 1);
  fi[2] = f_(jat + 1, iat);
  fi[3] = f_(jat + 1, iat + 1);
  fi[4] = ft_(jat, iat);
  fi[5] = ft_(jat, iat + 1);
  fi[6] = ft_(jat + 1, iat);
  fi[7] = ft_(jat + 1, iat + 1);
  fi[8] = ftt_(jat, iat);
  fi[9] = ftt_(jat, iat + 1);
  fi[10] = ftt_(jat + 1, iat);
  fi[11] = ftt_(jat + 1, iat + 1);
  fi[12] = fd_(jat, iat);
  fi[13] = fd_(jat, iat + 1);
  fi[14] = fd_(jat + 1, iat);
  fi[15] = fd_(jat + 1, iat + 1);
  fi[16] = fdd_(jat, iat);
  fi[17] = fdd_(jat, iat + 1);
  fi[18] = fdd_(jat + 1, iat);
  fi[19] = fdd_(jat + 1, iat + 1);
  fi[20] = fdt_(jat, iat);
  fi[21] = fdt_(jat, iat + 1);
  fi[22] = fdt_(jat + 1, iat);
  fi[23] = fdt_(jat + 1, iat + 1);
  fi[24] = fddt_(jat, iat);
  fi[25] = fddt_(jat, iat + 1);
  fi[26] = fddt_(jat + 1, iat);
  fi[27] = fddt_(jat + 1, iat + 1);
  fi[28] = fdtt_(jat, iat);
  fi[29] = fdtt_(jat, iat + 1);
  fi[30] = fdtt_(jat + 1, iat);
  fi[31] = fdtt_(jat + 1, iat + 1);
  fi[32] = fddtt_(jat, iat);
  fi[33] = fddtt_(jat, iat + 1);
  fi[34] = fddtt_(jat + 1, iat);
  fi[35] = fddtt_(jat + 1, iat + 1);

  // differences
  Real xt = std::max(0.0, (T - T_(jat)) * dti);
  Real xd = std::max(0.0, (De - rho_(iat)) * ddi);
  Real mxt = 1. - xt;
  Real mxd = 1. - xd;

  // Density and temperature basis functions
  Real si0t = hermite::psi0(xt);
  Real si1t = hermite::psi1(xt) * dth;
  Real si2t = hermite::psi2(xt) * dt2;

  Real si0mt = hermite::psi0(mxt);
  Real si1mt = -hermite::psi1(mxt) * dth;
  Real si2mt = hermite::psi2(mxt) * dt2;

  Real si0d = hermite::psi0(xd);
  Real si1d = hermite::psi1(xd) * dd;
  Real si2d = hermite::psi2(xd) * dd2;

  Real si0md = hermite::psi0(mxd);
  Real si1md = -hermite::psi1(mxd) * dd;
  Real si2md = hermite::psi2(mxd) * dd2;

  // derivatives of the weight functions
  Real dsi0t = hermite::dpsi0(xt) * dti;
  Real dsi1t = hermite::dpsi1(xt);
  Real dsi2t = hermite::dpsi2(xt) * dth;

  Real dsi0mt = -hermite::dpsi0(mxt) * dti;
  Real dsi1mt = hermite::dpsi1(mxt);
  Real dsi2mt = -hermite::dpsi2(mxt) * dth;

  Real dsi0d = hermite::dpsi0(xd) * ddi;
  Real dsi1d = hermite::dpsi1(xd);
  Real dsi2d = hermite::dpsi2(xd) * dd;

  Real dsi0md = -hermite::dpsi0(mxd) * ddi;
  Real dsi1md = hermite::dpsi1(mxd);
  Real dsi2md = -hermite::dpsi2(mxd) * dd;

  // second derivatives of the weight functions
  Real ddsi0t = hermite::ddpsi0(xt) * dt2i;
  Real ddsi1t = hermite::ddpsi1(xt) * dti;
  Real ddsi2t = hermite::ddpsi2(xt);

  Real ddsi0mt = hermite::ddpsi0(mxt) * dt2i;
  Real ddsi1mt = -hermite::ddpsi1(mxt) * dti;
  Real ddsi2mt = hermite::ddpsi2(mxt);

  // free energy
  Real free_en = hermite::h5(fi, si0t, si1t, si2t, si0mt, si1mt, si2mt, si0d, si1d, si2d,
                             si0md, si1md, si2md);
  // derivative with respect to temperature
  Real df_t = hermite::h5(fi, dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt, si0d, si1d,
                          si2d, si0md, si1md, si2md);
  // second derivative with respect to temperature
  Real df_tt = hermite::h5(fi, ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, si0d,
                           si1d, si2d, si0md, si1md, si2md);

  if (only_e) { // just set energy and return
    eele[0] = Ye * free_en + T * (-df_t) * Ye;
    eele[2] = T * (-df_tt) * Ye;
    return;
  }

  // Otherwise keep going
  // derivative with respect to density
  Real df_d = hermite::h5(fi, si0t, si1t, si2t, si0mt, si1mt, si2mt, dsi0d, dsi1d, dsi2d,
                          dsi0md, dsi1md, dsi2md);

  // derivative with respect to temperature and density
  Real df_dt = hermite::h5(fi, dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt, dsi0d, dsi1d,
                           dsi2d, dsi0md, dsi1md, dsi2md);

  // now get the pressure derivative with density, chemical potential, and
  // electron positron number densities
  // get the interpolation weight functions
  si0t = hermite::xpsi0(xt);
  si1t = hermite::xpsi1(xt) * dth;

  si0mt = hermite::xpsi0(mxt);
  si1mt = -hermite::xpsi1(mxt) * dth;

  si0d = hermite::xpsi0(xd);
  si1d = hermite::xpsi1(xd) * dd;

  si0md = hermite::xpsi0(mxd);
  si1md = -hermite::xpsi1(mxd) * dd;

  // derivatives of weight functions
  dsi0t = hermite::xdpsi0(xt) * dti;
  dsi1t = hermite::xdpsi1(xt);

  dsi0mt = -hermite::xdpsi0(mxt) * dti;
  dsi1mt = hermite::xdpsi1(mxt);

  dsi0d = hermite::xdpsi0(xd) * ddi;
  dsi1d = hermite::xdpsi1(xd);

  dsi0md = -hermite::xdpsi0(mxd) * ddi;
  dsi1md = hermite::xdpsi1(mxd);

  // Re-use cache
  fi[0] = dpdf_(jat, iat);
  fi[1] = dpdf_(jat, iat + 1);
  fi[2] = dpdf_(jat + 1, iat);
  fi[3] = dpdf_(jat + 1, iat + 1);
  fi[4] = dpdft_(jat, iat);
  fi[5] = dpdft_(jat, iat + 1);
  fi[6] = dpdft_(jat + 1, iat);
  fi[7] = dpdft_(jat + 1, iat + 1);
  fi[8] = dpdfd_(jat, iat);
  fi[9] = dpdfd_(jat, iat + 1);
  fi[10] = dpdfd_(jat + 1, iat);
  fi[11] = dpdfd_(jat + 1, iat + 1);
  fi[12] = dpdfdt_(jat, iat);
  fi[13] = dpdfdt_(jat, iat + 1);
  fi[14] = dpdfdt_(jat + 1, iat);
  fi[15] = dpdfdt_(jat + 1, iat + 1);

  // pressure derivative with respect to density
  pele[1] = std::max(
      0.0, Ye * hermite::h3(fi, si0t, si1t, si0mt, si1mt, si0d, si1d, si0md, si1md));

  // chemical potentials
  fi[0] = ef_(jat, iat);
  fi[1] = ef_(jat, iat + 1);
  fi[2] = ef_(jat + 1, iat);
  fi[3] = ef_(jat + 1, iat + 1);
  fi[4] = eft_(jat, iat);
  fi[5] = eft_(jat, iat + 1);
  fi[6] = eft_(jat + 1, iat);
  fi[7] = eft_(jat + 1, iat + 1);
  fi[8] = efd_(jat, iat);
  fi[9] = efd_(jat, iat + 1);
  fi[10] = efd_(jat + 1, iat);
  fi[11] = efd_(jat + 1, iat + 1);
  fi[12] = efdt_(jat, iat);
  fi[13] = efdt_(jat, iat + 1);
  fi[14] = efdt_(jat + 1, iat);
  fi[15] = efdt_(jat + 1, iat + 1);

  // electron chemical potential etaele
  etaele[0] = hermite::h3(fi, si0t, si1t, si0mt, si1mt, si0d, si1d, si0md, si1md);

  // derivative with respect to density
  Real x = hermite::h3(fi, si0t, si1t, si0mt, si1mt, dsi0d, dsi1d, dsi0md, dsi1md);
  etaele[1] = Ye * x;

  // derivative with respect to temperature
  etaele[2] = hermite::h3(fi, dsi0t, dsi1t, dsi0mt, dsi1mt, si0d, si1d, si0md, si1md);

  // derivative with respect to abar and zbar
  etaele[3] = -x * De * Ytot;
  etaele[4] = x * rho * Ytot;

  // look in the number density table only once
  fi[0] = xf_(jat, iat);
  fi[1] = xf_(jat, iat + 1);
  fi[2] = xf_(jat + 1, iat);
  fi[3] = xf_(jat + 1, iat + 1);
  fi[4] = xft_(jat, iat);
  fi[5] = xft_(jat, iat + 1);
  fi[6] = xft_(jat + 1, iat);
  fi[7] = xft_(jat + 1, iat + 1);
  fi[8] = xfd_(jat, iat);
  fi[9] = xfd_(jat, iat + 1);
  fi[10] = xfd_(jat + 1, iat);
  fi[11] = xfd_(jat + 1, iat + 1);
  fi[12] = xfdt_(jat, iat);
  fi[13] = xfdt_(jat, iat + 1);
  fi[14] = xfdt_(jat + 1, iat);
  fi[15] = xfdt_(jat + 1, iat + 1);

  // electron + positron number densities
  xne[0] = hermite::h3(fi, si0t, si1t, si0mt, si1mt, si0d, si1d, si0md, si1md);

  // derivative with respect to density
  x = std::max(0.0,
               hermite::h3(fi, si0t, si1t, si0mt, si1mt, dsi0d, dsi1d, dsi0md, dsi1md));
  xne[1] = Ye * x;

  // derivative with respect to temperature
  xne[2] = hermite::h3(fi, dsi0t, dsi1t, dsi0mt, dsi1mt, si0d, si1d, si0md, si1md);

  // derivative with respect to abar and zbar
  xne[3] = -x * De * Ytot;
  xne[4] = x * rho * Ytot;

  // the desired electron-positron thermodynamic quantities

  // dpepdd at high temperatures and low densities is below the
  // floating point limit of the subtraction of two large terms.
  // since dpresdd doesn't enter the maxwell relations at all, use the
  // bicubic interpolation done above instead of this one225
  x = De * De;
  pele[0] = x * df_d;
  pele[2] = x * df_dt;
  // pele[1]  = ye * (x * df_dd + 2.0 * din * df_d);
  Real s = pele[1] / Ye - 2.0 * De * df_d;
  pele[3] = -Ytot * (2.0 * pele[0] + s * De);
  pele[4] = rho * Ytot * (2.0 * De * df_d + s);

  x = Ye * Ye;
  sele[0] = -df_t * Ye;
  sele[2] = -df_tt * Ye;
  sele[1] = -df_dt * x;
  sele[3] = Ytot * (Ye * df_dt * De - sele[0]);
  sele[4] = -Ytot * (Ye * df_dt * rho + df_t);

  eele[0] = Ye * free_en + T * sele[0];
  eele[2] = T * sele[2];
  eele[1] = x * df_d + T * sele[1];
  eele[3] = -Ye * Ytot * (free_en + df_d * De) + T * sele[3];
  eele[4] = Ytot * (free_en + Ye * df_d * rho) + T * sele[4];
}

PORTABLE_INLINE_FUNCTION
void HelmRad::GetFromDensityTemperature(const Real rho, const Real temp,
                                        Real prad[NDERIV], Real erad[NDERIV],
                                        Real srad[NDERIV]) const {
  using namespace HelmUtils;
  const Real rhoi = robust::ratio(1.0, rho);
  const Real tempi = robust::ratio(1.0, temp);
  const Real T3 = math_utils::pow<3>(temp);
  const Real T4 = T3 * temp;
  // pressure
  Real P = PREFACTOR * T4;
  Real dPdT = 4.0 * PREFACTOR * T3;
  prad[VAL] = P;
  prad[DDR] = 0;
  prad[DDT] = dPdT;
  prad[DDA] = 0;
  prad[DDZ] = 0;

  // free energy
  const Real PoR = P * rhoi;
  Real e = 3.0 * PoR;
  Real dedR = -e * rhoi;
  Real dedT = 3.0 * dPdT * rhoi;
  erad[VAL] = e;
  erad[DDR] = dedR;
  erad[DDT] = dedT;
  erad[DDA] = 0;
  erad[DDZ] = 0;

  // entropy
  constexpr Real dPdR = 0;
  const Real dEdR = erad[DDR];
  Real s = (PoR + e) * tempi;
  Real dsdR = ((dPdR - PoR) * rhoi + dEdR) * tempi;
  Real dsdT = (dPdT * rhoi + dedT - s) * tempi;
  srad[VAL] = s;
  srad[DDR] = dsdR;
  srad[DDT] = dsdT;
  srad[DDA] = 0;
  srad[DDZ] = 0;
}

PORTABLE_INLINE_FUNCTION
void HelmIon::GetFromDensityTemperature(const Real rho, const Real temp, const Real y,
                                        const Real ytot, const Real xni,
                                        const Real dxnidd, const Real dxnida,
                                        Real pion[NDERIV], Real eion[NDERIV],
                                        Real sion[NDERIV]) const {
  using namespace HelmUtils;
  const Real kbT = KB * temp;
  // pressure
  const Real P = xni * kbT;
  const Real dPdR = dxnidd * kbT;
  const Real dPdT = xni * KB;
  const Real dPdA = dxnida * kbT;
  pion[VAL] = P;
  pion[DDR] = dPdR;
  pion[DDT] = dPdT;
  pion[DDA] = dPdA;
  pion[DDZ] = 0;

  // energy
  const Real rhoi = robust::ratio(1.0, rho);
  const Real PoR = P * rhoi;
  const Real e = 1.5 * PoR;
  const Real dedR = (1.5 * dPdR - e) * rhoi;
  const Real dedT = 1.5 * dPdT * rhoi;
  const Real dedA = 1.5 * dPdA * rhoi;
  eion[VAL] = e;
  eion[DDR] = dedR;
  eion[DDT] = dedT;
  eion[DDA] = dedA;
  eion[DDZ] = 0;

  // entropy
  const Real tempi = robust::ratio(1.0, temp);
  const Real KBNAY = KBNA * ytot;
  const Real s = (PoR + e) * tempi + KBNAY * y;
  const Real dsdR = (dPdR * rhoi - P * rhoi * rhoi + dedR) * tempi - KBNAY * rhoi;
  const Real dsdT =
      (dPdT * rhoi + dedT) * tempi - (PoR + e) * tempi * tempi + 1.5 * KBNAY * tempi;
  const Real dsdA = (dPdA * rhoi + dedA) * tempi + KBNAY * ytot * (2.5 - y);
  sion[VAL] = s;
  sion[DDR] = dsdR;
  sion[DDT] = dsdT;
  sion[DDA] = dsdA;
  sion[DDZ] = 0;
}

PORTABLE_INLINE_FUNCTION
void HelmCoulomb::GetFromDensityTemperature(const Real rho, const Real temp,
                                            const Real ytot, const Real abar,
                                            const Real zbar, const Real xni,
                                            const Real dxnidd, const Real dxnida,
                                            Real pcoul[NDERIV], Real ecoul[NDERIV],
                                            Real scoul[NDERIV]) const {
  using namespace HelmUtils;
  // fitting parameters
  constexpr Real a1 = -0.898004;
  constexpr Real b1 = 0.96786;
  constexpr Real c1 = 0.220703;
  constexpr Real d1 = -0.86097;
  constexpr Real e1 = 2.5269;
  constexpr Real a2 = 0.29561;
  constexpr Real b2 = 1.9885;
  constexpr Real c2 = 0.288675;

  // TODO(JMM): This is computed right out of the ion class. should
  // this be unified with the ion code?
  const Real kbT = KB * temp;
  const Real pion = xni * kbT;
  const Real dpiondd = dxnidd * kbT;
  const Real dpiondt = xni * KB;
  const Real dpionda = dxnida * kbT;
  constexpr Real dpiondz = 0.0;

  constexpr Real four_thirds = 4.0 / 3.0;
  constexpr Real ftpi = four_thirds * M_PI;
  Real s = ftpi * xni;
  const Real dsdd = ftpi * dxnidd;
  const Real dsda = ftpi * dxnida;

  constexpr Real one_third = 1.0 / 3.0;
  const Real si = robust::ratio(1.0, s);
  const Real lami = robust::ratio(1.0, std::cbrt(s));
  const Real lamidd = -lami * one_third * dsdd * si;
  const Real lamida = -lami * one_third * dsda * si;

  const Real plasg = math_utils::pow<2>(ELECTRON_CHARGE_ESU * zbar) / (kbT * lami);
  const Real plasg_o_lami = robust::ratio(plasg, lami);
  const Real plasgdd = -plasg_o_lami * lamidd;
  const Real plasgda = -plasg_o_lami * lamida;
  const Real plasgdt = -robust::ratio(plasg, temp);
  const Real plasgdz = 2.0 * robust::ratio(plasg, zbar);

  // TODO(JMM): Might be able to vectorize this with masking but
  // probably not worth it.
  const Real abari = robust::ratio(1, abar);
  if (plasg >= 0) {
    const Real x = std::sqrt(std::sqrt(plasg));
    Real y = KBNA * ytot;
    const Real c1_o_x = robust::ratio(c1, x);
    ecoul[VAL] = y * temp * (a1 * plasg + b1 * x + c1_o_x + d1);
    pcoul[VAL] = one_third * rho * ecoul[0];
    scoul[VAL] = -y * (3.0 * b1 * x - 5.0 * c1_o_x + d1 * (std::log(plasg) - 1.0) - e1);

    y = KBNA * temp * ytot * (a1 + robust::ratio(0.25, plasg) * (b1 * x - c1_o_x));
    ecoul[DDR] = y * plasgdd;
    ecoul[DDT] = y * plasgdt + robust::ratio(ecoul[0], temp);
    ecoul[DDA] = y * plasgda - ecoul[0] * abari;
    ecoul[DDZ] = y * plasgdz;

    y = one_third * rho;
    pcoul[DDR] = one_third * ecoul[0] + y * ecoul[1];
    pcoul[DDT] = y * ecoul[2];
    pcoul[DDA] = y * ecoul[3];
    pcoul[DDZ] = y * ecoul[4];

    y = -KBNA * plasg * abari * (0.75 * b1 * x + 1.25 * c1_o_x + d1);
    scoul[DDR] = y * plasgdd;
    scoul[DDT] = y * plasgdt;
    scoul[DDA] = y * plasgda - scoul[0] * abari;
    scoul[DDZ] = y * plasgdz;
  } else {
    const Real x = plasg * std::sqrt(std::abs(plasg));
    const Real y = std::pow(plasg, b2);
    const Real z = c2 * x - one_third * a2 * y;

    const Real plasgi = robust::ratio(1., plasg);
    const Real rhoi = robust::ratio(1.0, rho);

    pcoul[VAL] = -pion * z;
    ecoul[VAL] = 3.0 * pcoul[0] * rhoi;
    scoul[VAL] = -KBNA * abari * (c2 * x - a2 * robust::ratio(b2 - 1.0, b2) * y);

    s = 1.5 * c2 * x * plasgi - one_third * a2 * b2 * y * plasgi;
    pcoul[DDR] = -dpiondd * z - pion * s * plasgdd;
    pcoul[DDT] = -dpiondt * z - pion * s * plasgdt;
    pcoul[DDA] = -dpionda * z - pion * s * plasgda;
    pcoul[DDZ] = -dpiondz * z - pion * s * plasgdz;

    s = 3.0 * rhoi;
    ecoul[DDR] = s * pcoul[1] - ecoul[0] * rhoi;
    ecoul[DDT] = s * pcoul[2];
    ecoul[DDA] = s * pcoul[3];
    ecoul[DDZ] = s * pcoul[4];

    s = -KBNA * abari * plasgi * (1.5 * c2 * x - a2 * (b2 - 1.0) * y);
    scoul[DDR] = s * plasgdd;
    scoul[DDT] = s * plasgdt;
    scoul[DDA] = s * plasgda - scoul[0] * abari;
    scoul[DDZ] = s * plasgdz;
  }

  // butterworth bomb proofing by ^_^ : beware the butterbomb
  Real g_r, dgdf_r, g_t, dgdf_t;
  butterworth<12>(log10(temp) - 4.5, 3.0, g_r, dgdf_r);
  butterworth<12>(log10(rho) + 1.0, 6.0, g_t, dgdf_t);

  // derivatives (and conversion from logarithmic derivative)
  Real gain = (1.0 - g_t * g_r);
  Real dgaindt = -robust::ratio(g_r * dgdf_t, temp * LN10);
  Real dgaindd = -robust::ratio(g_t * dgdf_r, rho * LN10);

  // straight up gain
  pcoul[VAL] = pcoul[VAL] * gain;
  ecoul[VAL] = ecoul[VAL] * gain;
  scoul[VAL] = scoul[VAL] * gain;

  // derivatives via chain rule
  pcoul[DDR] = gain * pcoul[DDR] + robust::ratio(pcoul[VAL] * dgaindd, gain);
  pcoul[DDT] = gain * pcoul[DDT] + robust::ratio(pcoul[VAL] * dgaindt, gain);
  pcoul[DDA] = gain * pcoul[DDA];
  pcoul[DDZ] = gain * pcoul[DDZ];

  ecoul[DDR] = gain * ecoul[DDR] + robust::ratio(ecoul[VAL] * dgaindd, gain);
  ecoul[DDT] = gain * ecoul[DDT] + robust::ratio(ecoul[VAL] * dgaindt, gain);
  ecoul[DDA] = gain * ecoul[DDA];
  ecoul[DDZ] = gain * ecoul[DDZ];

  scoul[DDR] = gain * scoul[DDR] + robust::ratio(scoul[VAL] * dgaindd, gain);
  scoul[DDT] = gain * scoul[DDT] + robust::ratio(scoul[VAL] * dgaindt, gain);
  scoul[DDA] = gain * scoul[DDA];
  scoul[DDZ] = gain * scoul[DDZ];
}

PORTABLE_INLINE_FUNCTION
void Helmholtz::GetFromDensityLogTemperature_(
    const Real rho, const Real T, const Real abar, const Real zbar, const Real ye,
    const Real ytot, const Real ywot, const Real De, const Real lDe, Real p[NDERIV],
    Real e[NDERIV], Real s[NDERIV], Real etaele[NDERIV], Real nep[NDERIV],
    // TODO(JMM): Decide which of the quantities below to keep
    const bool only_e) const {
  double prad[NDERIV] = {0}, pion[NDERIV] = {0}, pele[NDERIV] = {0}, pcoul[NDERIV] = {0};
  double erad[NDERIV] = {0}, eion[NDERIV] = {0}, eele[NDERIV] = {0}, ecoul[NDERIV] = {0};
  double srad[NDERIV] = {0}, sion[NDERIV] = {0}, sele[NDERIV] = {0}, scoul[NDERIV] = {0};
  // electron chemical potential and electron + positron number density
  for (int i = 0; i < 5; ++i) {
    etaele[i] = nep[i] = 0;
  }
  Real lT = std::log10(T);

  const Real log10e = std::log10(M_E);
  const Real lnT = lT / log10e;
  if (options_.ENABLE_RAD) {
    rad_.GetFromDensityTemperature(rho, T, prad, erad, srad);
  }
  if (options_.ENABLE_GAS) {
    // If gas is not ionized, just use ideal gas for ions
    // If gas is ionized but the electrons are not degenerate, modify ideal gas
    // coefficient If gas is ionized and the electrons are degenerate, do ideal gas for
    // ions and add degerate fermi gas
    bool do_unmodified_ions =
        (options_.GAS_IONIZED && options_.GAS_DEGENERATE) || !options_.GAS_IONIZED;
    if (do_unmodified_ions) { // ionized plasma. ions and electrons treated separately
      Real xni, dxnidd, dxnida;
      GetMassFractions(rho, T, ytot, xni, dxnidd, dxnida);
      const Real y = ywot + ions_.LSWOT15 * lnT;
      ions_.GetFromDensityTemperature(rho, T, y, ytot, xni, dxnidd, dxnida, pion, eion,
                                      sion);
    } else { // modify ideal gas to include ions + electrons
      const Real abar_ion = robust::ratio(abar, zbar + 1);
      // const Real zbar_ion = zbar;
      const Real ytot_ion = robust::ratio(1.0, abar_ion);
      const Real ywot_ion = robust::ratio(
          std::log(abar_ion * abar_ion * std::sqrt(abar_ion)), rho * ions_.NA);
      const Real y_ion = ywot_ion + ions_.LSWOT15 * lnT;
      Real xni, dxnidd, dxnida;
      GetMassFractions(rho, T, ytot_ion, xni, dxnidd, dxnida);
      ions_.GetFromDensityTemperature(rho, T, y_ion, ytot_ion, xni, dxnidd, dxnida, pion,
                                      eion, sion);
    }
    if (options_.GAS_DEGENERATE) { // treat degenerate electron gas
      electrons_.GetFromDensityTemperature(rho, lT, ye, ytot, De, lDe, pele, eele, sele,
                                           etaele, nep, only_e);
      if (options_.ENABLE_COULOMB_CORRECTIONS) {
        Real xni, dxnidd, dxnida;
        GetMassFractions(rho, T, ytot, xni, dxnidd, dxnida);
        coul_.GetFromDensityTemperature(rho, T, ytot, abar, zbar, xni, dxnidd, dxnida,
                                        pcoul, ecoul, scoul);
      }
    }
  }
  for (int i = 0; i < 5; ++i) {
    e[i] = erad[i] + eion[i] + eele[i] + ecoul[i];
  }
  if (!only_e) {
    for (int i = 0; i < 5; ++i) {
      p[i] = prad[i] + pion[i] + pele[i] + pcoul[i];
      s[i] = srad[i] + sion[i] + sele[i] + scoul[i];
    }
  }
}

}; // namespace singularity

#undef ROOT_FINDER
#endif // SINGULARITY_USE_HELMHOLTZ
#endif // SINGULARITY_USE_SPINER_WITH_HDF5
#endif // _SINGULARITY_EOS_EOS_HELMHOLTZ_HPP_
