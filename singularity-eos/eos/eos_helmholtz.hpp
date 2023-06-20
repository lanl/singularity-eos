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

// ports of call
#include <ports-of-call/portability.hpp>

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
 private:
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

}; // namespace singularity

#endif // SPINER_USE_HDF
#endif // _SINGULARITY_EOS_EOS_HELMHOLTZ_HPP_
