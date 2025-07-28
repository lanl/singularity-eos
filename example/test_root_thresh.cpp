#include <iostream>
#include <string>
#include <cstdlib>
#include <chrono>

// Spiner EOS
#include <singularity-eos/eos/eos.hpp>
#include <singularity-eos/eos/eos_spiner_rho_temp.hpp>

//Root finder
#include <singularity-eos/base/root-finding-1d/root_finding.hpp>


using namespace singularity;
namespace RootFinding1D {
//extern bool test_root_thresh;

template <typename T>
PORTABLE_INLINE_FUNCTION Status regula_falsi(const T &f, const Real ytarget,
                                             const Real guess, Real a, Real b,
                                             const Real xtol, const Real ytol,
                                             Real &xroot,
                                             const RootCounts *counts,
                                             const bool &verbose);
}
int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " sp5_file.h5 matid\n";
    return 1;
  }

  std::string sp5_file = argv[1];
  int matid = std::atoi(argv[2]);
  double rho = 0.1;
  double sie = 1.0e12;

  singularity::SpinerEOSDependsRhoT eos(sp5_file, matid);
  singularity::TableStatus status;

  auto start = std::chrono::high_resolution_clock::now();
  double temperature = eos.TemperatureFromDensityInternalEnergy(rho, sie);
  auto end = std::chrono::high_resolution_clock::now();
  double elapsed = std::chrono::duration<double>(end - start).count();

  eos.Root_Thresh_Print(); //prints the tolerance
  RootFinding1D::test_root_thresh = true;
  std::cout << "Temperature = " << temperature << " K\n";
  std::cout << "Elapsed time = " << elapsed << " s\n";

  return 0;
}
