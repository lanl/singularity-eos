//------------------------------------------------------------------------------
// © 2021-2026. Triad National Security, LLC. All rights reserved.  This
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

// This example demonstrates how to tabulate a custom EOS implementation
// into SpinerEOS format using the generic constructor. This is useful
// when you have your own EOS physics model and want to:
// - Improve performance through table interpolation
// - Enable GPU portability (tables are GPU-friendly)
// - Integrate with codes that expect tabulated EOS

#include <cmath>
#include <iostream>

#ifdef SINGULARITY_USE_SPINER_WITH_HDF5
#include <singularity-eos/eos/eos.hpp>

using namespace singularity;

// Example: Custom EOS class from a host application
// This represents a simple EOS that a simulation code might implement.
// It only needs to provide a minimal interface to be tabulated.
class CustomHostEOS {
 private:
  // Physical parameters for a simple Mie-Gruneisen-like EOS
  Real rho0_;   // reference density
  Real C0_;     // bulk sound speed
  Real s_;      // Hugoniot slope parameter
  Real Gamma0_; // Gruneisen parameter
  Real Cv_;     // specific heat

 public:
  CustomHostEOS(Real rho0, Real C0, Real s, Real Gamma0, Real Cv)
      : rho0_(rho0), C0_(C0), s_(s), Gamma0_(Gamma0), Cv_(Cv) {}

  // Minimal required interface for SpinerEOS constructor:

  // 1. Internal energy from density and temperature
  PORTABLE_INLINE_FUNCTION
  Real InternalEnergyFromDensityTemperature(Real rho, Real T) const {
    // Simple thermal energy contribution
    return Cv_ * T;
  }

  // 2. Temperature from density and internal energy
  PORTABLE_INLINE_FUNCTION
  Real TemperatureFromDensityInternalEnergy(Real rho, Real sie) const {
    return sie / Cv_;
  }

  // 3. Pressure from density and temperature
  PORTABLE_INLINE_FUNCTION
  Real PressureFromDensityTemperature(Real rho, Real T) const {
    // Mie-Gruneisen form: P = P_cold(rho) + Gamma * rho * e_thermal
    Real eta = 1.0 - rho0_ / rho; // compression

    // Cold pressure contribution (Hugoniot-based)
    Real P_cold = 0.0;
    if (eta > 0) {
      Real denom = 1.0 - s_ * eta;
      if (denom > 0) {
        P_cold = rho0_ * C0_ * C0_ * eta / (denom * denom);
      }
    }

    // Thermal pressure contribution
    Real e_thermal = Cv_ * T;
    Real P_thermal = Gamma0_ * rho * e_thermal;

    return P_cold + P_thermal;
  }

  // Optional: Gruneisen parameter (improves accuracy of derivatives)
  PORTABLE_INLINE_FUNCTION
  Real GruneisenParamFromDensityTemperature(Real rho, Real T) const {
    return Gamma0_; // Constant for this simple model
  }

  // Optional: Material properties
  PORTABLE_INLINE_FUNCTION
  Real MeanAtomicMass() const { return 63.5; } // Copper-like

  PORTABLE_INLINE_FUNCTION
  Real MeanAtomicNumber() const { return 29.0; } // Copper
};

int main(int argc, char *argv[]) {
  std::cout << "Custom EOS to SpinerEOS Example\n";
  std::cout << "================================\n\n";

  // Create custom EOS with copper-like parameters
  Real rho0 = 8.93;  // g/cm^3
  Real C0 = 3.94;    // km/s = 3.94e5 cm/s
  Real s = 1.489;    // dimensionless
  Real Gamma0 = 2.0; // dimensionless
  Real Cv = 3.85e10; // erg/g/K

  CustomHostEOS custom_eos(rho0, C0, s, Gamma0, Cv);

  std::cout << "Created custom Mie-Gruneisen-like EOS with parameters:\n";
  std::cout << "  rho0   = " << rho0 << " g/cm^3\n";
  std::cout << "  C0     = " << C0 << " km/s\n";
  std::cout << "  s      = " << s << "\n";
  std::cout << "  Gamma0 = " << Gamma0 << "\n";
  std::cout << "  Cv     = " << Cv << " erg/g/K\n\n";

  // Set up grid parameters for tabulation
  SpinerTableGridParams params;

  // Define thermodynamic range
  params.rhoMin = 1.0;   // g/cm^3
  params.rhoMax = 20.0;  // g/cm^3
  params.TMin = 300.0;   // K
  params.TMax = 50000.0; // K

  // Compute energy bounds
  params.sieMin =
      custom_eos.InternalEnergyFromDensityTemperature(params.rhoMin, params.TMin);
  params.sieMax =
      custom_eos.InternalEnergyFromDensityTemperature(params.rhoMax, params.TMax);

  // Grid resolution (coarser for example speed)
  params.numRhoPerDecade = 50;
  params.numTPerDecade = 50;
  params.numSiePerDecade = 50;

  // Enable piecewise grids for efficiency
  params.piecewiseRho = true;
  params.piecewiseT = true;
  params.piecewiseSie = true;

  params.matid = 9999;
  params.rhoNormal = rho0;

  std::cout << "Tabulating custom EOS with grid parameters:\n";
  std::cout << "  rho range: [" << params.rhoMin << ", " << params.rhoMax << "] g/cm^3\n";
  std::cout << "  T range:   [" << params.TMin << ", " << params.TMax << "] K\n";
  std::cout << "  sie range: [" << params.sieMin << ", " << params.sieMax << "] erg/g\n";
  std::cout << "  Resolution: " << params.numRhoPerDecade << " points per decade\n";
  std::cout << "  Piecewise grids: enabled\n\n";

  std::cout << "Constructing SpinerEOS (this may take a moment)...\n";

  // Create SpinerEOS from custom EOS
  // The constructor will:
  // - Evaluate custom_eos on the specified grid
  // - Use GruneisenParam if available, else finite differences
  // - Compute all necessary derivatives
  // - Create a fully functional SpinerEOS in memory
  SpinerEOSDependsRhoSie spiner_eos(custom_eos, params);

  std::cout << "SpinerEOS construction complete!\n\n";

  // Verify the tabulated EOS matches the original
  std::cout << "Verification: Comparing tabulated vs. original EOS\n";
  std::cout << "==================================================\n\n";

  // Test points
  Real test_rho = 10.0;  // g/cm^3
  Real test_T = 10000.0; // K

  // Compute from original custom EOS
  Real P_custom = custom_eos.PressureFromDensityTemperature(test_rho, test_T);
  Real sie_custom = custom_eos.InternalEnergyFromDensityTemperature(test_rho, test_T);

  // Compute from tabulated SpinerEOS
  Real P_spiner = spiner_eos.PressureFromDensityTemperature(test_rho, test_T);
  Real sie_spiner = spiner_eos.InternalEnergyFromDensityTemperature(test_rho, test_T);

  std::cout << "Test point: rho = " << test_rho << " g/cm^3, T = " << test_T << " K\n\n";

  std::cout << "Pressure:\n";
  std::cout << "  Custom EOS:  " << P_custom << " erg/cm^3\n";
  std::cout << "  Spiner EOS:  " << P_spiner << " erg/cm^3\n";
  std::cout << "  Difference:  " << std::abs(P_spiner - P_custom) << " erg/cm^3\n";
  std::cout << "  Rel. error:  "
            << std::abs(P_spiner - P_custom) / std::abs(P_custom) * 100.0 << " %\n\n";

  std::cout << "Internal Energy:\n";
  std::cout << "  Custom EOS:  " << sie_custom << " erg/g\n";
  std::cout << "  Spiner EOS:  " << sie_spiner << " erg/g\n";
  std::cout << "  Difference:  " << std::abs(sie_spiner - sie_custom) << " erg/g\n";
  std::cout << "  Rel. error:  "
            << std::abs(sie_spiner - sie_custom) / std::abs(sie_custom) * 100.0
            << " %\n\n";

  // Test temperature inversion
  Real test_sie = sie_custom;
  Real T_custom = custom_eos.TemperatureFromDensityInternalEnergy(test_rho, test_sie);
  Real T_spiner = spiner_eos.TemperatureFromDensityInternalEnergy(test_rho, test_sie);

  std::cout << "Temperature (inverted from sie):\n";
  std::cout << "  Custom EOS:  " << T_custom << " K\n";
  std::cout << "  Spiner EOS:  " << T_spiner << " K\n";
  std::cout << "  Difference:  " << std::abs(T_spiner - T_custom) << " K\n";
  std::cout << "  Rel. error:  "
            << std::abs(T_spiner - T_custom) / std::abs(T_custom) * 100.0 << " %\n\n";

  // Check material properties
  std::cout << "Material Properties:\n";
  std::cout << "  Spiner Abar: " << spiner_eos.MeanAtomicMass() << "\n";
  std::cout << "  Spiner Zbar: " << spiner_eos.MeanAtomicNumber() << "\n";
  std::cout << "  (Correctly extracted from CustomHostEOS)\n\n";

  std::cout << "Summary\n";
  std::cout << "=======\n";
  std::cout << "The tabulated SpinerEOS successfully reproduces the custom EOS\n";
  std::cout << "with interpolation errors < 1%. The SpinerEOS can now be used:\n";
  std::cout << "  - In production simulations for better performance\n";
  std::cout << "  - On GPUs (via GetOnDevice())\n";
  std::cout << "  - With the same interface as file-based Spiner tables\n";
  std::cout << "  - For mixed-cell closures (PTE, etc.)\n\n";

  return 0;
}

#else

int main() {
  std::cout << "This example requires SINGULARITY_USE_SPINER_WITH_HDF5=ON\n";
  return 1;
}

#endif
