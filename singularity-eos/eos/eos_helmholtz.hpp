//------------------------------------------------------------------------------
// The Helmholtz EOS was first presented in: Timmes and Swesty, ApJS
// 126:501-516 (2000).  Reference implementation in fortan is here:
// https://cococubed.com/code_pages/eos.shtml
// Reference implementation is CC Frank TImmes and Douglas Swesty, 1999.
// Original work is open-sourced under the CC-By license
// https://creativecommons.org/licenses/by/4.0/
//------------------------------------------------------------------------------
// © 2023. Triad National Security, LLC. All rights reserved.  This
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
#ifdef SPINER_USE_HDF

#include <cstdio>

#include <fstream>
#include <sstream>
#include <string>
#include <utility>

// ports of call
#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>

// singularity-eos
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/eos_base.hpp>
#include <singularity-eos/base/hermite.hpp>
#include <singularity-eos/base/math_utils.hpp>
#include <singularity-eos/base/robust_utils.hpp>

// spiner
#include <spiner/databox.hpp>

namespace singularity {
using namespace eos_base;

class Helmholtz : public EosBase<Helmholtz> {
 public:
  // These are the indexes in the lambda array
  struct Lambda {
    enum Index {
      Ye = 0,   // electron fraction zbar/abar
      Ytot = 1, // 1.0 / abar
      De = 2,   // Electron density rho * Ye
      lDe = 3,  // log10(De)
      lT = 4 // log10 temperature. used for root finding.
    }; 
  };

  inline Helmholtz(const std::string &filename) {
    InitDataFile_(filename);
  }
  inline Helmholtz(const std::string &data_file, const std::string &species_file) {
    InitDataFile_(data_life);
    InitSpeciesFile_(species_file);
  }

  inline Helmholtz GetOnDevice();  
  inline void Finalize();

 private:
  inline void InitDataFile_(const std::string &filename);
  inline void InitSpeciesFile_(const std::string &filename) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Stub");
  }

  PORTABLE_INLINE_FUNCTION
  void GetAll_(Real rho, Real T, Real lT, Real Ye, Real Ytot, Real De, Real lDe,
               double pele[5], double eele[5], double sele[5], double etaele[5],
               double xne[5], bool only_e = false) const;

  // Tail-recursive resize tables
  inline void ResizeTables_(Spiner::DataBox &db) {
    db.resize(NTEMP, NRHO); // set shape and log bounds
    db.setRange(1, lTMin_, lTMax_, NTEMP);
    db.setRange(0, lRhoMin_, lRhoMax_, NRHO);
  }
  template <typename... Args>
  inline void ResizeTables_(Spiner::DataBox &head, Args &&...tail) {
    ResizeTables_(head);
    ResizeTables_(std::forward<Args>(tail)...);
  }
  // Tail-recursive read one i,j from text file
  template <typename Msg_t, typename T>
  inline void Read_(std::ifstream &file, Msg_t &error_msg, T &var) {
    if (!(file >> var)) {
      file.close(); // is this needed?
      PORTABLE_ALWAYS_THROW_OR_ABORT(error_msg);
    }
  }
  template <typename Msg_t, typename... Args>
  inline void Read_(std::ifstream &file, Msg_t &error_msg, T &head, Args &&...tail) {
    Read_(file, msg, head);
    Read_(file, msg, std::forward<Args>(tail)...);
  }
  // Read all i,j from text file
  template <typename... Args>
  inline void SetTablesFromFile_(std::ifstream &file, Args &&...tables) {
    ResizeTables(std::forward<Args>(tables)...);
    for (int j = 0; j < NTEMP; ++j) {
      for (int i = 0; i < NTEMP; ++i) {
        std::stringstream error_msg;
        msg << "Error reading the Helmholtz free energy table at j = " << j
            << ", i = " << i << std::endl;
        Read_(file, error_msg, tables(j, i)...);
      }
    }
  }

  // rho and T caches (to go between log/linear scale)
  Spiner::DataBox rho_, T_;
  // Free energy and derivatives
  Spiner::DataBox f_, fd_, ft_, fdd_, ftt_, fdt_, fddt_, fdtt_, fddtt_;
  // derivatives
  Spiner::DataBox dpdf_, dpdfd_, dpdft_, dpdfdt_;
  // chemical potential
  Spiner::DataBox ef_, efd_, eft_, efdt_;
  // number density
  Spiner::DataBox xf_, xfd_, xft_, xfdt_;
  // TODO(JMM): Add species file

  constexpr std::size_t NTEMP = 101;
  constexpr std::size_t NRHO = 271;

  constexpr Real lTMin_ = 3.0;
  constexpr Real lTMax_ = 13.0;
  constexpr Real lRhoMin_ = -12.0;
  constexpr Real lRhoMax_ = 15.0;

  constexpr Real dlT_ = (lTMax_ - lTMin_) / (static_cast<Real>(NTEMP) - 1.0);
  constexpr Real dlRho_ = (lRhoMax_ - lRhoMin_) / (static_cast<Real>(NRHO) - 1.0);

  const Real TMin_ = math_utils::pow10(lTMin_);
  const Real TMax_ = math_utils::pow10(lTMax_);
  const Real rhoMin_ = math_utils::pow10(lRhoMin_);
  const Real rhoMax_ = math_utils::pow10(rhoMax_);
};

inline void Helmholtz::InitDataFile_(const std::string &filename) {
  std::ifstream file(filename);
  PORTABLE_ALWAYS_REQUIRE(file.is_open(), "Helmholtz file " + filename + " not found!");
  SetTablesFromFile_(file, f_, fd_, ft_, fdd_, ftt_, fdt_, fddt_, fdtt_, fddtt_);
  SetTablesFromFile_(file, dpdf_, dpdfd_, dpdft_, dpdfdt_);
  SetTablesFromFile_(file, ef_, efd_, eft_, efdt_);
  SetTablesFromFile_(file, xf_, xfd_, xft_, xfdt_);
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
}

inline Helmholtz Helmholtz::GetOnDevice() {
  Helmholtz other;
  other.rho_ = Spiner::getOnDeviceDataBox(rho_);
  other.T_ = Spiner::getOnDeviceDataBox(T_);
  other.f_ = Spiner::getOnDeviceDataBox(f_);
  other.fd_ = Spiner::getOnDeviceDataBox(fd_);
  other.ft_ = Spiner::getOnDeviceDataBox(ft_);
  other.fdd_ = Spiner::getOnDeviceDataBox(fdd_);
  other.ftt_ = Spiner::getOnDeviceDataBox(ftt_);
  other.fdt_ = Spiner::getOnDeviceDataBox(fdt_);
  other.fdd_ = Spiner::getOnDeviceDataBox(fdd_);
  other.fdtt_ = Spiner::getOnDeviceDataBox(fdtt_);
  other.fddtt_ = Spiner::getOnDeviceDataBox(fddtt_);
  other.dpdf_ = Spiner::getOnDeviceDataBox(dpdf_);
  other.dpdfd_ = Spiner::getOnDeviceDataBox(dpdfd_);
  other.dpdft_ = Spiner::getOnDeviceDataBox(dpdft_);
  other.dpdfdt_ = Spiner::getOnDeviceDataBox(dpdfdt_);
  other.ef_ = Spiner::getOnDeviceDataBox(ef_);
  other.efd_ = Spiner::getOnDeviceDataBox(efd_);
  other.eft_ = Spiner::getOnDeviceDataBox(eftf_);
  other.efdtf_ = Spiner::getOnDeviceDataBox(efdt_);
  other.xf_ = Spiner::getOnDeviceDataBox(xf_);
  other.xfd_ = Spiner::getOnDeviceDataBox(xfd_);
  other.xft_ = Spiner::getOnDeviceDataBox(xft_);
  other.xfdt_ = Spiner::getOnDeviceDataBox(xfdt_);
}

inline void Helmholtz::Finalize() {
  rho_.finalize();
  T_.finalize();
  f_.finalize();
  fd_.finalize();
  ft_.finalize();
  fdd_.finalize();
  ftt_.finalize();
  fdt_.finalize();
  fdd_.finalize();
  fdtt_.finalize();
  fddtt_.finalize();
  dpdf_.finalize();
  dpdfd_.finalize();
  dpdft_.finalize();
  dpdfdt_.finalize();
  ef_.finalize();
  efd_.finalize();
  eft_.finalize();
  efdtf_.finalize();
  xf_.finalize();
  xfd_.finalize();
  xft_.finalize();
  xfdt_.finalize();
}

PORTABLE_INLINE_FUNCTION
void Helmholtz::GetAll_(Real rho, Real T, Real lT, Real Ye, Real Ytot, Real De, Real lDe,
                        double pele[5], double eele[5], double sele[5], double etaele[5],
                        double xne[5], bool only_e = false) const {
  // Bound lRho, lT
  rho = std::min(rhoMax_, std::max(rhoMin_, rho));
  De = std::min(rhoMax_, std::max(rhoMin_, De));
  lDe = std::min(lRhoMax_, std::max(lRhoMin_, lDe));
  T = std::min(TMax_, std::max(TMin_, T));
  lT = std::min(lTMax_, std::max(lTMin_, lT));

  // Find central indexes in table
  auto &lRhoRange = f_.range(0);
  auto &lTRange = f_.range(1);
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
  fi[0] = f_[jat][iat];
  fi[1] = f_[jat + 1][iat];
  fi[2] = f_[jat][iat + 1];
  fi[3] = f_[jat + 1][iat + 1];
  fi[4] = ft_[jat][iat];
  fi[5] = ft_[jat + 1][iat];
  fi[6] = ft_[jat][iat + 1];
  fi[7] = ft_[jat + 1][iat + 1];
  fi[8] = ftt_[jat][iat];
  fi[9] = ftt_[jat + 1][iat];
  fi[10] = ftt_[jat][iat + 1];
  fi[11] = ftt_[jat + 1][iat + 1];
  fi[12] = fd_[jat][iat];
  fi[13] = fd_[jat + 1][iat];
  fi[14] = fd_[jat][iat + 1];
  fi[15] = fd_[jat + 1][iat + 1];
  fi[16] = fdd_[jat][iat];
  fi[17] = fdd_[jat + 1][iat];
  fi[18] = fdd_[jat][iat + 1];
  fi[19] = fdd_[jat + 1][iat + 1];
  fi[20] = fdt_[jat][iat];
  fi[21] = fdt_[jat + 1][iat];
  fi[22] = fdt_[jat][iat + 1];
  fi[23] = fdt_[jat + 1][iat + 1];
  fi[24] = fddt_[jat][iat];
  fi[25] = fddt_[jat + 1][iat];
  fi[26] = fddt_[jat][iat + 1];
  fi[27] = fddt_[jat + 1][iat + 1];
  fi[28] = fdtt_[jat][iat];
  fi[29] = fdtt_[jat + 1][iat];
  fi[30] = fdtt_[jat][iat + 1];
  fi[31] = fdtt_[jat + 1][iat + 1];
  fi[32] = fddtt_[jat][iat];
  fi[33] = fddtt_[jat + 1][iat];
  fi[34] = fddtt_[jat][iat + 1];
  fi[35] = fddtt_[jat + 1][iat + 1];

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
  fi[0] = dpdf_[jat][iat];
  fi[1] = dpdf_[jat + 1][iat];
  fi[2] = dpdf_[jat][iat + 1];
  fi[3] = dpdf_[jat + 1][iat + 1];
  fi[4] = dpdft_[jat][iat];
  fi[5] = dpdft_[jat + 1][iat];
  fi[6] = dpdft_[jat][iat + 1];
  fi[7] = dpdft_[jat + 1][iat + 1];
  fi[8] = dpdfd_[jat][iat];
  fi[9] = dpdfd_[jat + 1][iat];
  fi[10] = dpdfd_[jat][iat + 1];
  fi[11] = dpdfd_[jat + 1][iat + 1];
  fi[12] = dpdfdt_[jat][iat];
  fi[13] = dpdfdt_[jat + 1][iat];
  fi[14] = dpdfdt_[jat][iat + 1];
  fi[15] = dpdfdt_[jat + 1][iat + 1];

  // pressure derivative with respect to density
  pele[1] = std::max(
      0.0, Ye * hermite::h3(fi, si0t, si1t, si0mt, si1mt, si0d, si1d, si0md, si1md));

  // chemical potentials
  fi[0] = ef_[jat][iat];
  fi[1] = ef_[jat + 1][iat];
  fi[2] = ef_[jat][iat + 1];
  fi[3] = ef_[jat + 1][iat + 1];
  fi[4] = eft_[jat][iat];
  fi[5] = eft_[jat + 1][iat];
  fi[6] = eft_[jat][iat + 1];
  fi[7] = eft_[jat + 1][iat + 1];
  fi[8] = efd_[jat][iat];
  fi[9] = efd_[jat + 1][iat];
  fi[10] = efd_[jat][iat + 1];
  fi[11] = efd_[jat + 1][iat + 1];
  fi[12] = efdt_[jat][iat];
  fi[13] = efdt_[jat + 1][iat];
  fi[14] = efdt_[jat][iat + 1];
  fi[15] = efdt_[jat + 1][iat + 1];

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
  fi[0] = xf_[jat][iat];
  fi[1] = xf_[jat + 1][iat];
  fi[2] = xf_[jat][iat + 1];
  fi[3] = xf_[jat + 1][iat + 1];
  fi[4] = xft_[jat][iat];
  fi[5] = xft_[jat + 1][iat];
  fi[6] = xft_[jat][iat + 1];
  fi[7] = xft_[jat + 1][iat + 1];
  fi[8] = xfd_[jat][iat];
  fi[9] = xfd_[jat + 1][iat];
  fi[10] = xfd_[jat][iat + 1];
  fi[11] = xfd_[jat + 1][iat + 1];
  fi[12] = xfdt_[jat][iat];
  fi[13] = xfdt_[jat + 1][iat];
  fi[14] = xfdt_[jat][iat + 1];
  fi[15] = xfdt_[jat + 1][iat + 1];

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
  // bicubic interpolation done above instead of this one
  x = De * De;
  pele[0] = x * df_d;
  pele[2] = x * df_dt;
  // pele[1]  = ye * (x * df_dd + 2.0 * din * df_d);
  s = pele[1] / Ye - 2.0 * De * df_d;
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

}; // namespace singularity

#endif // SPINER_USE_HDF
#endif // _SINGULARITY_EOS_EOS_HELMHOLTZ_HPP_
