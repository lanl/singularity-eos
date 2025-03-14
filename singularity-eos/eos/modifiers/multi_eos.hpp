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

#ifndef _SINGULARITY_EOS_EOS_MULTI_EOS_
#define _SINGULARITY_EOS_EOS_MULTI_EOS_

#include "stdio.h"
#include <cstdlib>
#include <iostream>
#include <tuple>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/eos_error.hpp>
#include <singularity-eos/eos/eos_base.hpp>

namespace singularity {

using namespace eos_base;

template<typename ClosureMethod, typename... EOSModelsT>
class MultiEOS : public EosBase<MultiEOS<ClosureMethod, EOSModelsT...>> {
 public:
  SG_ADD_BASE_CLASS_USINGS(MultiEOS)

  static std::string EosType() {
    std::string all_models;
    // Use a fold expression to append the models and put commas between each model
    ((all_models += (all_models.empty() ? "" : ", ") + EOSModelsT::ModelType()), ...);
    return "MultiEOS<"ClosureMethod::MethodType() + ", " + all_models + ">";
  }

  static std::string EosPyType() {
    // Use a fold expression to append each model
    std::string all_models = (std::string() + ... + EOSModelsT::ModelType());
    return "MultiEOS" + ClosureMethod::MethodType() + all_models;
  }

  MultiEOS() = default;

  template<typename... EOSModelsT_>
  TwoEOS(EOSModelsT_... eos_models)
      : models_(std::forward<EOSModelsT_>(eos_models))
  {}

 private:
  std::tuple<EOSModelsT...> models_;
};

#endif // #ifndef _SINGULARITY_EOS_EOS_MULTI_EOS_
