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
#ifndef _SINGULARITY_EOS_EOS_EOS_SPINER_SIE_TRANSFORMS_HPP_
#define _SINGULARITY_EOS_EOS_EOS_SPINER_SIE_TRANSFORMS_HPP_

#ifdef SINGULARITY_USE_SPINER_WITH_HDF5

#include <type_traits>

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
// #include <iostream> // debug
// #include <stdio.h> // debug

// ports-of-call
#include <ports-of-call/portability.hpp>

// base
#include <singularity-eos/eos/eos_spiner_common.hpp>

namespace singularity {
namespace transformations {

template <typename Data = void>
struct NullTransform {

  template <typename... Args>
  constexpr NullTransform(Args &&...) {}

  template <typename... Args>
  constexpr auto transform(const Real e, Args &&...) const {
    return e;
  }

  template <typename... Args>
  constexpr auto inverse(const Real e_transformed, Args &&...) const {
    return e_transformed;
  }
};

template <typename Data>
struct ShiftTransform {

  template <typename DataT_in>
  constexpr ShiftTransform(const DataT_in &data) : data_(data) {}

 private:
  Data data_;

 public:
  template <typename... Args>
  constexpr auto transform(const Real e, const Real rho, Args &&...) const {
    const Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);
    const Real e_cold = data_.sieCold.interpToReal(lRho);
    return e - e_cold;
  }

  template <typename... Args>
  constexpr auto inverse(const Real e_transformed, const Real rho, Args &&...) const {
    const Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);
    const Real e_cold = data_.sieCold.interpToReal(lRho);
    return e_transformed + e_cold;
  }
};

// Divide by the heat capacity
// TODO: Get working on device
template <typename Data>
struct DivideByCvTransform {

  template <typename DataT_in>
  constexpr DivideByCvTransform(const DataT_in &data) : data_(data) {}

 private:
  Data data_;
  static constexpr Real min_Cv_ = 1.0e-08;
  static constexpr Real min_iCv_ = 1.0e-08;

 public:
  template <typename... Args>
  constexpr auto transform(Real e, Real rho, Args &&...) const {
    const Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);
    const Real lE = spiner_common::to_log(e, data_.lEOffset);
    Real iCv = data_.dTdE.interpToReal(lRho, lE);
    iCv = iCv > min_iCv_ ? iCv : min_iCv_;
    return e * iCv;
  }

  template <typename... Args>
  constexpr auto inverse(Real e_transformed, Real rho, Real e_orig, Args &&...) const {
    const Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);
    const Real lE = spiner_common::to_log(e_orig, data_.lEOffset);
    Real Cv = 1. / data_.dTdE.interpToReal(lRho, lE);
    Cv = Cv > min_Cv_ ? Cv : min_Cv_;
    return e_transformed * Cv;
  }
};

// TODO: Get working on device
template <typename Data>
struct ShiftandDivideByCvTransform {

  template <typename DataT_in>
  constexpr ShiftandDivideByCvTransform(const DataT_in &data) : data_(data) {}

 private:
  Data data_;
  static constexpr Real min_Cv_ = 1.0e-08;
  static constexpr Real min_iCv_ = 1.0e-08;

 public:
  template <typename... Args>
  constexpr auto transform(Real e, Real rho, Args &&...) const {
    const Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);
    const Real e_cold = data_.sieCold.interpToReal(lRho);
    const Real lE = spiner_common::to_log(e, data_.lEOffset);
    Real e_coldshift = e - e_cold;
    Real iCv = data_.dTdE.interpToReal(lRho, lE);
    iCv = iCv > min_iCv_ ? iCv : min_iCv_;
    return e_coldshift * iCv;
  }

  template <typename... Args>
  constexpr auto inverse(Real e_transformed, Real rho, Real e_orig, Args &&...) const {
    const Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);
    const Real e_cold = data_.sieCold.interpToReal(lRho);
    const Real lE = spiner_common::to_log(e_orig, data_.lEOffset);
    Real Cv = 1. / data_.dTdE.interpToReal(lRho, lE);
    Cv = Cv > min_Cv_ ? Cv : min_Cv_;
    const Real e_inverse_cv = e_transformed * Cv;
    return e_inverse_cv + e_cold;
  }
};

// TO DO: T table would depend on rho and sie, need to change in order to work
// Divide by T^alpha
// TODO: Get working on device
template <typename Data>
struct ScaleTransform {

  template <typename DataT_in>
  constexpr ScaleTransform(const DataT_in &data) : data_(data) {}

 private:
  Data data_;
  static constexpr Real min_T_ = 1.0e-08; // Good minimum temperature?

 public:
  template <typename... Args>
  constexpr auto transform(Real e, Real rho, Args &&...) const {
    Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);
    Real T = data_.T.interpToReal(lRho); // T?
    T = T > min_T_ ? T : min_T_;
    return e / pow(T, 3);
  }

  template <typename... Args>
  constexpr auto inverse(Real e, Real rho, Args &&...) const {
    Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);
    Real T = data_.T.interpToReal(lRho);
    T = T > min_T_ ? T : min_T_;
    return e * pow(T, 3);
  }
};

// TO DO: Do to T table depending on rho and sie, this transformation does
// not work as intended
// TODO: Get working on device
template <typename Data>
struct AllTransform {

  template <typename DataT_in>
  constexpr AllTransform(const DataT_in &data) : data_(data) {}

 private:
  Data data_;
  static constexpr Real min_T_ = 1.0e-08; // Good minimum temperature?
  static constexpr Real min_Cv_ = 1.0e-08;
  static constexpr Real min_iCv_ = 1.0e-08;

 public:
  template <typename... Args>
  constexpr auto transform(Real e, Real rho, Args &&...) const {

    Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);

    Real e_cold = data_.sieCold.interpToReal(lRho);
    Real T = data_.T.interpToReal(lRho);
    T = T > min_T_ ? T : min_T_;

    Real e_transformed = e - e_cold;
    Real lE = spiner_common::to_log(e_transformed, data_.lEOffset);
    Real iCv = data_.dTdE.interpToReal(lRho, lE);
    iCv = iCv > min_iCv_ ? iCv : min_iCv_;
    e_transformed = e_transformed * iCv;

    return e_transformed / pow(T, 3);
  }

  template <typename... Args>
  constexpr auto inverse(Real e_transformed, Real rho, Real e_orig, Args &&...) const {

    Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);
    Real e_cold = data_.sieCold.interpToReal(lRho);
    Real T = data_.T.interpToReal(lRho);
    T = T > min_T_ ? T : min_T_;

    Real e_transformed_inverse = e_orig - e_cold;
    Real lE = spiner_common::to_log(e_transformed_inverse, data_.lEOffset);
    Real Cv = 1. / data_.dTdE.interpToReal(lRho, lE);

    e_transformed = e_transformed * pow(T, 3);

    Cv = Cv > min_Cv_ ? Cv : min_Cv_;
    e_transformed = e_transformed * Cv;

    return e_transformed + e_cold;
  }
};

} // namespace transformations
} // namespace singularity

#endif // SINGULARITY_USE_SPINER_WITH_HDF5
#endif // _SINGULARITY_EOS_EOS_EOS_SPINER_SIE_TRANSFORMS_HPP_
