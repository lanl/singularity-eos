//------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------

// C++ headers
#include <cmath>
#include <iostream>

// Needed to import the eos models
#include <singularity-eos/eos/eos.hpp>
// One way of initializing models with modifiers
#include <singularity-eos/eos/eos_builder.hpp>

// For simplicity. Most things live in the singularity namespace
using namespace singularity;

constexpr double SMALL = 1e-20; // to avoid dividing by zero

// Loop through a flattened list of cells and fill pressure and sound
// speed using a given EOS. We make the call two different ways, once
// with single calls, and once with the FillEos call, which should be
// more performant.
inline void PressureSoundSpeedFromDensityEnergyDensity(double *rho, // inputs
                                                       double *uu, EOS &eos,
                                                       double *P, // outputs
                                                       double *cs, const int Ncells) {
  // Allocate lambda object, which is used for caching by the tables.
  // Lambdas can do other things, such as take extra parameters as well.
  //
  // Technically this object should be allocated outside this function
  // and passed in, so that it can be persistent. But it's here for
  // example's sake.
  //
  // For parallel execution, you need one lambda vector per thread.
  // Since we're running in serial, we'll just use one.
  const int nlambda = eos.nlambda(); // get number of elements per lambda
  std::vector<double> lambda(nlambda);

  // Loop through the cells and use the two function calls
  for (int i = 0; i < Ncells; ++i) {
    double sie = uu[i] / (rho[i] + SMALL); // convert to specific internal energy
    P[i] = eos.PressureFromDensityInternalEnergy(rho[i], sie, lambda.data());
    double bmod = eos.BulkModulusFromDensityInternalEnergy(rho[i], sie, lambda.data());
    cs[i] = std::sqrt(bmod / (rho[i] + SMALL));
  }

  /*
   * request the output variables pressure and bulk modulus.
   * Use bitwise or to build the bitfield of quantities you want.
   * Available variables for input/output are:
   *
   * density, specific_internal_energy, pressure,
   * temperature, specific_heat, bulk_modulus
   *
   * Note that all quantities not requested as output are considered
   * inputs.  Moreover, most equations of state calls require density
   * and either temperature or energy as one of their inputs.
   * Therefore, to tell the machinery that we want to use energy as an
   * input, not temperature, we must also request temperature as an
   * output.
   */
  constexpr unsigned long output =
      (thermalqs::temperature | thermalqs::pressure | thermalqs::bulk_modulus);

  // Loop through cells and use the FillEos function call
  for (int i = 0; i < Ncells; ++i) {
    double eps, temp, cv;
    // FillEos is very general and is capable of modifying any of the inputs,
    // so const vars cannot be passed into it. However, it is often more performant
    // than making individual function calls.
    Real sie = uu[i] / (rho[i] + SMALL); // convert to specific internal energy
    eos.FillEos(rho[i], temp, sie, P[i], cv, cs[i], output, lambda.data());
    // convert bulk modulus to cs
    cs[i] = std::sqrt(cs[i] / (rho[i] + SMALL));
  }
}

int main() {

  // Parameters for ideal gas
  constexpr double gm1 = 0.6;
  constexpr double Cv = 2;

  // We initialize the eos two ways to show how it can be done
  // First initialize the EOS the "easy" way:
  EOS eos1 = IdealGas(gm1, Cv);

  // The EOS Builder automates building an eos with modifiers, such as
  // shift and scale by worrying about the combinatorics/switch
  // statements, etc, for you. We'll initialize a shifted/scaled EOS
  // with shift 0 and scale 1 to show you how it's done
  constexpr double shift = 0;
  constexpr double scale = 1;

  // Eos builder requires a type specification,
  // and then for each modifier, it requires a set of parameters,
  // which are implemented as a variatic dictionary.
  //
  // This is equivalent to
  // EOS eos2 = Shifted(Scaled(IdealGas(gm1, Cv), scale), shift);
  EOSBuilder::EOSType type = EOSBuilder::EOSType::IdealGas;
  EOSBuilder::modifiers_t modifiers;                               // this is a dictionary
  EOSBuilder::params_t base_params, shifted_params, scaled_params; // so are these
  base_params["Cv"].emplace<Real>(
      Cv); // base params is the parameters of the underlying eos
  base_params["gm1"].emplace<Real>(gm1);
  shifted_params["shift"].emplace<Real>(
      shift); // these are the parameters for the modifiers
  scaled_params["scale"].emplace<Real>(scale); // you use strings for the variable names
  // for each modifier put the relevant params in the modifiers object
  modifiers[EOSBuilder::EOSModifier::Shifted] = shifted_params;
  modifiers[EOSBuilder::EOSModifier::Scaled] = scaled_params;
  EOS eos2 = EOSBuilder::buildEOS(type, base_params, modifiers); // build the builder

  // If you're on device, you need to call
  // eos1.GetOnDevice();
  // to call the EOS on device
  // However, since we're doing everything on host here, we don't make this call

  // Make some arrays
  constexpr int N = 50;
  std::vector<double> rho(N);
  std::vector<double> uu(N);
  std::vector<double> P(N);
  std::vector<double> cs(N);

  // Here we fill the rho and uu arrays with something sensible.
  for (int i = 0; i < N; ++i) {
    rho[i] = 1 + 0.1 * std::sin(2 * M_PI * i / static_cast<double>(N));
    uu[i] = 1e-2 * rho[i];
  }

  // Call it!
  PressureSoundSpeedFromDensityEnergyDensity(rho.data(), uu.data(), eos1, P.data(),
                                             cs.data(), N);

  // And let's print out the final value just for fun
  std::cout << "The final values are:\n"
            << "rho = " << rho[N - 1] << "\n"
            << "uu  = " << uu[N - 1] << "\n"
            << "P   = " << P[N - 1] << "\n"
            << "cs  = " << cs[N - 1] << std::endl;

  // It's good practice to call Finalize() after you're done using an EOS.
  // This usually only does anything if you're on device, but for GPU-data it
  // will call the appropriate destructor.
  eos1.Finalize();
  eos2.Finalize();

  return 0;
}
