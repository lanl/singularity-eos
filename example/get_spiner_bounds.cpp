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

// C++ headers
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

// Needed to import the eos models
#include <singularity-eos/eos/eos.hpp>

int main(int argc, char *argv[]) {
  // In case we're doing a Kokkos backend
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
#endif
  { // begin kokkos scope
    if (argc < 3) {
      std::cerr << "Usage: " << argv[0] << " matid sp5_filename" << std::endl;
      // std::cerr << "Arguments given: " << argc << std::endl;
      std::exit(1);
    }

    const int matid = std::atoi(argv[1]);
    const std::string filename = argv[2];

    singularity::SpinerEOSDependsRhoT deprt(filename, matid);
    singularity::SpinerEOSDependsRhoSie depre(filename, matid);

    printf("rho bounds, depends(r, T): %.14e %.14e\n", deprt.MinimumDensity(),
           deprt.rhoMax());
    printf("rho bounds, depends(r, sie): %.14e %.14e\n", depre.MinimumDensity(),
           depre.rhoMax());
    printf("T bounds, depends(r, T): %.14e %.14e\n", deprt.MinimumTemperature(),
           deprt.TMax());
    printf("T bounds, depends(r, sie): %.14e %.14e\n", depre.MinimumTemperature(),
           depre.TMax());

    deprt.Finalize();
    depre.Finalize();
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif
}
