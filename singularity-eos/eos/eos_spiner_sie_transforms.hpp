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
  PORTABLE_INLINE_FUNCTION auto transform(const Real e, Args &&...) const {
    return e;
  }

  template <typename... Args>
  PORTABLE_INLINE_FUNCTION auto inverse(const Real e_transformed, Args &&...) const {
    return e_transformed;
  }
};

template <typename Data>
struct ShiftTransform {
  PORTABLE_INLINE_FUNCTION  ShiftTransform() = default;
  
  template <typename DataT_in>
  PORTABLE_INLINE_FUNCTION ShiftTransform(const DataT_in &data) : data_(data) {}

 private:
  Data data_;

 public:
  template <typename... Args>
  PORTABLE_INLINE_FUNCTION auto transform(const Real e, const Real rho,
                                          Args &&...) const {
    const Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);
    const Real e_cold = data_.sieCold.interpToReal(lRho);
    return e - e_cold;
  }

  template <typename... Args>
  PORTABLE_INLINE_FUNCTION auto inverse(const Real e_transformed, const Real rho,
                                        Args &&...) const {
    const Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);   
    const Real e_cold = data_.sieCold.interpToReal(lRho);
    return e_transformed + e_cold;
  }
};



//Divide by the heat capacity
template<typename Data>
struct DivideByCvTransform {
  PORTABLE_INLINE_FUNCTION  DivideByCvTransform() = default;
  
  template <typename DataT_in>  
  PORTABLE_INLINE_FUNCTION DivideByCvTransform(const DataT_in &data) : data_(data) {}

  private: 
    Data data_;
    Real min_Cv_ = 1.0e-08;

  public:
    template <typename... Args>
    PORTABLE_INLINE_FUNCTION
    auto transform(Real e, Real rho, Args &&...) const {
        const Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);
	const Real lE   = spiner_common::to_log(e, data_.lEOffset);
        Real Cv = 1./data_.dTdE.interpToReal(lRho, lE);
        Cv = std::max(Cv, min_Cv_);
        return e / Cv;
    }

    template <typename... Args>
    PORTABLE_INLINE_FUNCTION
    auto inverse(Real e_transformed, Real rho, Real e_orig, Args &&...) const {
        const Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);
	const Real lE   = spiner_common::to_log(e_orig, data_.lEOffset);
	Real Cv = 1./data_.dTdE.interpToReal(lRho, lE);
        Cv = std::max(Cv, min_Cv_);
        return e_transformed * Cv;
    }
};



//Divide by T^alpha
template<typename Data>
struct ScaleTransform {
    PORTABLE_INLINE_FUNCTION  ScaleTransform() = default;
    
    template <typename DataT_in>
PORTABLE_INLINE_FUNCTION ScaleTransform(const DataT_in &data) : data_(data) {}

private:
    Data data_;
    Real min_T_ = 1.0e-08; //Good minimum temperature? 

public: 
    template <typename... Args>
    PORTABLE_INLINE_FUNCTION
    auto transform(Real e, Real rho, Args &&...) const {
        Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);
        Real T = data_.T.interpToReal(lRho); //T?
        T = std::max(T, min_T_);
        return e / pow(T, 3);
    }

    template <typename... Args>
    PORTABLE_INLINE_FUNCTION
    auto inverse(Real e, Real rho, Args &&...) const {
        Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);
        Real T = data_.T.interpToReal(lRho);
        T = std::max(T, min_T_);
        return e * pow(T, 3);
    }
};

template<typename Data>
struct AllTransform {
PORTABLE_INLINE_FUNCTION  AllTransform() = default;
   
template <typename DataT_in>
PORTABLE_INLINE_FUNCTION AllTransform(const DataT_in &data) : data_(data) {}



private:
    Data data_;
    Real min_T_ = 1.0e-08; //Good minimum temperature? 
    Real min_Cv_ = 1.0e-08;


public:
    template <typename... Args>
    PORTABLE_INLINE_FUNCTION
     auto transform(Real e, Real rho, Args &&...) const {

     Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);

     Real e_cold = data_.sieCold.interpToReal(lRho);
     Real T = data_.T.interpToReal(lRho);
     T = std::max(T, min_T_);

     Real e_transformed = e - e_cold;
     Real lE   = spiner_common::to_log(e_transformed, data_.lEOffset);
     Real Cv = 1./data_.dTdE.interpToReal(lRho, lE);
     Cv = std::max(Cv, min_Cv_);
     e_transformed = e_transformed / Cv;

     return e_transformed / pow(T,3);
     }

    template <typename... Args>
    PORTABLE_INLINE_FUNCTION
    auto inverse(Real e_transformed, Real rho, Real e_orig, Args &&...) const {
        
	Real lRho = spiner_common::to_log(rho, data_.lRhoOffset);
	Real e_cold = data_.sieCold.interpToReal(lRho);
        Real T = data_.T.interpToReal(lRho);
        T = std::max(T, min_T_);

	  
        Real e_transformed_inverse = e_orig - e_cold;
        Real lE   = spiner_common::to_log(e_transformed_inverse, data_.lEOffset);
	Real Cv = 1./data_.dTdE.interpToReal(lRho, lE);

        e_transformed = e_transformed * pow(T, 3);

        Cv = std::max(Cv, min_Cv_);
        e_transformed = e_transformed * Cv;

	return e_transformed + e_cold;
	
        
    }
};






} // namespace singularity

#endif // SINGULARITY_USE_SPINER_WITH_HDF5
#endif // _SINGULARITY_EOS_EOS_EOS_SPINER_SIE_TRANSFORMS_HPP_
