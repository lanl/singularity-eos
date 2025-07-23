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

template <typename Data = void>
struct NullTransform {

  template <typename... Args>
  PORTABLE_INLINE_FUNCTION NullTransform(Args &&...) {}

  PORTABLE_INLINE_FUNCTION
  NullTransform() = default;

  template <typename... Args>
  PORTABLE_INLINE_FUNCTION auto transform(Real e, Args &&...) const {
    return e;
  }

  template <typename... Args>
  PORTABLE_INLINE_FUNCTION auto inverse(Real e_transformed, Args &&...) const {
    return e_transformed;
  }
};


template <typename Data>
struct ShiftTransform {
  template <typename DataT_in>
  PORTABLE_INLINE_FUNCTION ShiftTransform(const DataT_in &data)
      : data_{std::forward<DataT_in>(data)} {}

 private:
  Data data_;

 public:
  template <typename... Args>
  PORTABLE_INLINE_FUNCTION auto transform(Real e, Real rho, Args &&...) const {
    Real lRho = toLog_(rho, data_.lRhoOffset);
    Real e_cold = data_.sieCold.interpToReal(lRho);
    return e - e_cold;
  }

  template <typename... Args>
  PORTABLE_INLINE_FUNCTION auto inverse(Real e_transformed, Real lRho, Args &&...) const {
    Real e_cold = data_.sieCold.interpToReal(lRho);
    return e_transformed + e_cold;
  }
};

} //namespace singularity


#endif // SINGULARITY_USE_SPINER_WITH_HDF5
#endif // _SINGULARITY_EOS_EOS_EOS_SPINER_SIE_TRANSFORMS_HPP_
