// © 2026. Triad National Security, LLC. All rights reserved.
#include <cstdlib>
#include <iostream>
#include <string>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/eos/eos_sesame.hpp>

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "Usage: get_native_sesame_state SESAME_FILE MATID RHO TEMPERATURE\n";
    return 1;
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::ScopeGuard runtime(argc, argv);
#endif
  singularity::SesameEOS eos(argv[1], std::stoi(argv[2]));
  const Real rho = std::stod(argv[3]), temp = std::stod(argv[4]);
  std::cout << "e = " << eos.InternalEnergyFromDensityTemperature(rho, temp)
            << "\nP = " << eos.PressureFromDensityTemperature(rho, temp)
            << "\nCv = " << eos.SpecificHeatFromDensityTemperature(rho, temp)
            << "\nBs = " << eos.BulkModulusFromDensityTemperature(rho, temp)
            << "\nGamma = " << eos.GruneisenParamFromDensityTemperature(rho, temp)
            << '\n';
  eos.Finalize();
  return 0;
}
