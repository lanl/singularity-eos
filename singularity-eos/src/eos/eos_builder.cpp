//------------------------------------------------------------------------------#
// © 2021. Triad National Security, LLC. All rights reserved.  This
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
//------------------------------------------------------------------------------#


#include <iostream> // debug
#include <eos/eos_builder.hpp>

#define GETBASE(T,N) T N = mpark::get<T>(base_params[#N])

namespace singularity {

EOS EOSBuilder::buildEOS(EOSBuilder::EOSType type,
					   EOSBuilder::params_t base_params,
					   EOSBuilder::modifiers_t modifiers) {
  using namespace EOSBuilder;
  using std::string;
  bool scaled  = (modifiers.count(EOSModifier::Scaled) > 0);
  bool shifted = (modifiers.count(EOSModifier::Shifted) > 0);
  bool relativistic = (modifiers.count(EOSModifier::Relativistic) > 0);
  if ((shifted || scaled || relativistic) && !(isModifiable(type))) {
    EOS_ERROR("Modifiers not supported for this EOS");
  }
  if (relativistic && (shifted || scaled)) {
    EOS_ERROR("Relativistic modifier cannot be combined with shift or scale");
  }
  Real scale = 1;
  if (scaled) {
    scale  = mpark::get<Real>(modifiers[EOSModifier::Scaled]["scale"]);
  }
  Real shift = 0;
  if (shifted) {
    shift = mpark::get<Real>(modifiers[EOSModifier::Shifted]["shift"]);
  }
  Real cl = 2.99792458e8; // speed of light in cm/s
  if (relativistic && (modifiers[EOSModifier::Relativistic].count("cl") > 0)) {
    cl = mpark::get<Real>(modifiers[EOSModifier::Relativistic]["cl"]);
  }
  if (type == EOSType::IdealGas) {
    Real gm1 = mpark::get<Real>(base_params["gm1"]);
    Real Cv = mpark::get<Real>(base_params["Cv"]);
    IdealGas g(gm1,Cv);
    if (relativistic) {
      return makeRelativistic(std::move(g),cl);
    }
    return applyScaleAndShift(std::move(g),
			      scaled,shifted,scale,shift);
  }
#ifdef SPINER_USE_HDF
  if (type == EOSType::SpinerEOSDependsRhoT
      || type == EOSType::SpinerEOSDependsRhoSie) {
    string filename = mpark::get<string>(base_params["filename"]);
    bool reproducibility_mode = mpark::get<bool>(base_params["reproducibility_mode"]);
    if (base_params.count("matid") > 0) {
      int matid = mpark::get<int>(base_params["matid"]);
      if (type == EOSType::SpinerEOSDependsRhoT) {
        SpinerEOSDependsRhoT s(filename,matid,reproducibility_mode);
        if (relativistic) {
          return makeRelativistic(std::move(s),cl);
        }
        return applyScaleAndShift(std::move(s),
                                  scaled,shifted,scale,shift);
      } 
      else {
        SpinerEOSDependsRhoSie s(filename,matid,reproducibility_mode);
        if (relativistic) {
          return makeRelativistic(std::move(s),cl);
        }
        return applyScaleAndShift(std::move(s),
                                  scaled,shifted,scale,shift);
      }
    } else {
      string materialName = mpark::get<string>(base_params["materialName"]);
      if (type == EOSType::SpinerEOSDependsRhoT) {
        SpinerEOSDependsRhoT s(filename,materialName,reproducibility_mode);
        if (relativistic) {
          return makeRelativistic(std::move(s),cl);
        }
        return applyScaleAndShift(std::move(s),
                                  scaled,shifted,scale,shift);
      }
        else {
        SpinerEOSDependsRhoSie s(filename,materialName,reproducibility_mode);
        if (relativistic) {
          return makeRelativistic(std::move(s),cl);
        }
        return applyScaleAndShift(std::move(s),
                                  scaled,shifted,scale,shift);
      }
    }
  }
  if (type == EOSType::StellarCollapse) {
    string filename = mpark::get<string>(base_params["filename"]);
    bool use_sp5 = false;
    if (base_params.count("use_sp5") > 0) {
      use_sp5 = mpark::get<bool>(base_params["use_sp5"]);
    }
    bool filter_bmod = true;
    if (base_params.count("filter_bmod") > 0) {
      filter_bmod = mpark::get<bool>(base_params["filter_bmod"]);
    }
    StellarCollapse s(filename, use_sp5, filter_bmod);
    if (relativistic) {
      return makeRelativistic(std::move(s), cl);
    }
    return applyScaleAndShift(std::move(s),
                              scaled, shifted, scale, shift);
  }
#endif
  if (type == EOSType::Gruneisen) {
    Real C0   = mpark::get<Real>(base_params["C0"]);
    Real s1   = mpark::get<Real>(base_params["s1"]);
    Real s2   = mpark::get<Real>(base_params["s2"]);
    Real s3   = mpark::get<Real>(base_params["s3"]);
    Real G0   = mpark::get<Real>(base_params["G0"]);
    Real b    = mpark::get<Real>(base_params["b"]);
    Real rho0 = mpark::get<Real>(base_params["rho0"]);
    Real T0   = mpark::get<Real>(base_params["T0"]);
    Real P0   = mpark::get<Real>(base_params["P0"]);
    Real Cv   = mpark::get<Real>(base_params["Cv"]);
    return Gruneisen(C0,s1,s2,s3,G0,b,rho0,T0,P0,Cv);
  }
  if (type == EOSType::JWL) {
    GETBASE(Real,A); // I got tired of writing this over and over
    GETBASE(Real,B);
    GETBASE(Real,R1);
    GETBASE(Real,R2);
    GETBASE(Real,w);
    GETBASE(Real,rho0);
    GETBASE(Real,Cv);
    return JWL(A,B,R1,R2,w,rho0,Cv);
  }
  if (type == EOSType::DavisReactants) {
    GETBASE(Real,rho0);
    GETBASE(Real,e0);
    GETBASE(Real,P0);
    GETBASE(Real,T0);
    GETBASE(Real,A);
    GETBASE(Real,B);
    GETBASE(Real,C);
    GETBASE(Real,G0);
    GETBASE(Real,Z);
    GETBASE(Real,alpha);
    GETBASE(Real,Cv0);
    return DavisReactants(rho0,e0,P0,T0,A,B,C,G0,Z,alpha,Cv0);
  }
  if (type == EOSType::DavisProducts) {
    GETBASE(Real,a);
    GETBASE(Real,b);
    GETBASE(Real,k);
    GETBASE(Real,n);
    GETBASE(Real,vc);
    GETBASE(Real,pc);
    GETBASE(Real,Cv);
    GETBASE(Real,E0);
    return DavisProducts(a,b,k,n,vc,pc,Cv,E0);
  }
  exit(1);
}

bool EOSBuilder::isModifiable(EOSBuilder::EOSType type) {
  return EOSBuilder::modifiable.count(type) > 0;
}

} // namespace singularity
