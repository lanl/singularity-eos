//------------------------------------------------------------------------------
// © 2021-2022. Triad National Security, LLC. All rights reserved.  This
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

#ifdef SPINER_USE_HDF
#ifdef SINGULARITY_USE_EOSPAC

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/eos/eos.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>

using duration = std::chrono::nanoseconds;
constexpr int NTIMES = 10;

class Bounds {
 public:
  Bounds() {}

  Bounds(Real min, Real max, int N, Real offset)
      : grid(Spiner::RegularGrid1D(min, max, N)), offset(offset) {}

  Bounds(Real min, Real max, int N) : offset(0) {
    constexpr Real epsilon = std::numeric_limits<float>::epsilon();
    const Real min_offset = 10 * std::abs(epsilon);
    if (min <= 0) offset = 1.1 * std::abs(min) + min_offset;
    min += offset;
    max += offset;
    min = std::log10(std::abs(min));
    max = std::log10(std::abs(max));
    grid = Spiner::RegularGrid1D(min, max, N);
  }

  PORTABLE_INLINE_FUNCTION Real log2lin(Real xl) const {
    return std::pow(10., xl) - offset;
  }
  PORTABLE_INLINE_FUNCTION Real i2lin(int i) const { return log2lin(grid.x(i)); }

  friend std::ostream &operator<<(std::ostream &os, const Bounds &b) {
    os << "Bounds: [" << b.grid.min() << ", " << b.grid.max() << "]"
       << " + " << b.offset << ", "
       << "[N,dx] = [" << b.grid.nPoints() << ", " << b.grid.dx() << "]"
       << "\n";
    return os;
  }

 public:
  Spiner::RegularGrid1D grid;
  Real offset;
};

int main(int argc, char *argv[]) {

#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize(argc, argv);
#endif
  {
    if (argc < 5) {
      std::cout << "Usage: " << argv[0] << " filename nfine matid {0,1}" << std::endl;
      std::exit(1);
    }

    const std::string filename = argv[1];
    const int nFine = std::atoi(argv[2]);
    const int matid = std::atoi(argv[3]);
    const bool cache = std::atoi(argv[4]);

    const int nFineRho = nFine + 1;
    const int nFineT = nFine - 1;

    std::cout << "Profile root finding" << std::endl;
    std::cout << "\tLoading EOS..." << std::endl;
    singularity::SpinerEOSDependsRhoT eos_host(filename, matid);
    auto eos_device = eos_host.GetOnDevice();

    std::cout << "\tallocating scratch memory..." << std::endl;
    std::unique_ptr<Spiner::DataBox, Spiner::DBDeleter> plambda_host(
        new Spiner::DataBox(nFineRho, nFineT, eos_host.nlambda()));
    auto lambda_host = *plambda_host;
    std::unique_ptr<Spiner::DataBox, Spiner::DBDeleter> plambda_device(
        new Spiner::DataBox(Spiner::AllocationTarget::Device, nFineRho, nFineT,
                            eos_host.nlambda()));
    auto lambda_device = *plambda_device;
    std::unique_ptr<Spiner::DataBox, Spiner::DBDeleter> presults_host(
        new Spiner::DataBox(nFineRho, nFineT));
    auto results_host = *presults_host;
    std::unique_ptr<Spiner::DataBox, Spiner::DBDeleter> presults_device(
        new Spiner::DataBox(Spiner::AllocationTarget::Device, nFineRho, nFineT));
    auto results_device = *presults_device;
    std::unique_ptr<Spiner::DataBox, Spiner::DBDeleter> presults_hm(
        new Spiner::DataBox(nFineRho, nFineT));
    auto results_hm = *presults_hm;

    std::cout << "Filling bounds arrays..." << std::endl;
    const int rhoMin = (1 + 1e-8) * eos_host.rhoMin();
    const int rhoMax = (1 - 1e-8) * eos_host.rhoMax();
    const int TMin = (1 + 1e-8) * eos_host.TMin();
    const int TMax = (1 - 1e-8) * eos_host.TMax();
    Bounds lRhoBounds(rhoMin, rhoMax, nFineRho);
    Bounds lTBounds(TMin, TMax, nFineT);
    Spiner::DataBox rho_host(nFineRho);
    Spiner::DataBox sie_host(nFineRho, nFineT);
    for (int iRho = 0; iRho < nFineRho; ++iRho) {
      const Real rho = lRhoBounds.i2lin(iRho);
      rho_host(iRho) = rho;
      for (int iT = 0; iT < nFineT; ++iT) {
        const Real T = lTBounds.i2lin(iT);
        sie_host(iRho, iT) = eos_host.InternalEnergyFromDensityTemperature(rho, T);
      }
    }
    auto rho_device = rho_host.getOnDevice();
    auto sie_device = sie_host.getOnDevice();

    std::cout << "\tProfiling on host..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    for (int n = 0; n < NTIMES; ++n) {
      for (int j = 0; j < nFineRho; ++j) {
        Real rho = rho_host(j);
        for (int i = 0; i < nFineT; ++i) {
          Real sie = sie_host(j, i);
          Real *lambda = cache ? &(lambda_host(j, i, 0)) : nullptr;
          results_host(j, i) =
              eos_host.TemperatureFromDensityInternalEnergy(rho, sie, lambda);
        }
      }
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration_host = std::chrono::duration_cast<duration>(stop - start);

    std::cout << "\tProfiling on device..." << std::endl;
#ifdef PORTABILITY_STRATEGY_KOKKOS
    Kokkos::fence();
#endif
    start = std::chrono::high_resolution_clock::now();
    for (int n = 0; n < NTIMES; ++n) {
      portableFor(
          "sie from rho T", 0, nFineRho, 0, nFineT,
          PORTABLE_LAMBDA(const int j, const int i) {
            Real rho = rho_device(j);
            Real sie = sie_device(j, i);
            Real *lambda = cache ? &(lambda_device(j, i, 0)) : nullptr;
            results_device(j, i) =
                eos_device.TemperatureFromDensityInternalEnergy(rho, sie, lambda);
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
#endif
    }
    stop = std::chrono::high_resolution_clock::now();
    auto duration_device = std::chrono::duration_cast<duration>(stop - start);

    std::cout << "\tComputing error on host..." << std::endl;
    Real error_host = 0;
    for (int j = 0; j < nFineRho; ++j) {
      for (int i = 0; i < nFineT; ++i) {
        Real Troot = results_host(j, i);
        Real Ttrue = lTBounds.i2lin(i);
        Real diff = 2 * (Troot - Ttrue) / (std::abs(Troot) + std::abs(Ttrue) + 1e-10);
        error_host += diff * diff;
      }
    }
    error_host = std::sqrt(error_host) / (nFineRho * nFineT);

    std::cout << "\tCopying device results to host..." << std::endl;
#ifdef PORTABILITY_STRATEGY_KOKKOS
    {
      using HS = Kokkos::HostSpace;
      using DMS = Kokkos::DefaultExecutionSpace::memory_space;
      using memUnmanaged = Kokkos::MemoryUnmanaged;
      using HostView_t = Kokkos::View<Real *, HS, memUnmanaged>;
      using DeviceView_t = Kokkos::View<Real *, memUnmanaged>;
      DeviceView_t devView(results_device.data(), results_device.size());
      HostView_t hostView(results_hm.data(), results_hm.size());
      Kokkos::deep_copy(hostView, devView);
    }
#endif // PORTABILITY_STRATEGY_KOKKOS

    std::cout << "\tComparing host to device results.." << std::endl;
    Real diff_hd = 0;
    for (int j = 0; j < nFineRho; ++j) {
      for (int i = 0; i < nFineT; ++i) {
        Real host = results_host(j, i);
        Real device = results_hm(j, i);
        Real diff = 2 * (host - device) / (std::abs(host) + std::abs(device) + 1e-10);
        diff_hd += diff * diff;
      }
    }
    diff_hd = std::sqrt(diff_hd) / (nFineRho * nFineT);

    std::cout << "\t\tRoot finding:\n"
              << "\t\tits: counts\n";
    for (int i = 0; i < eos_host.counts.nBins(); i++) {
      std::cout << "\t\t\t" << i << ": "
                << 100. * eos_host.counts[i] / eos_host.counts.total() << "\n";
    }
    printf("Results:\n"
           "\terror on host = %.14e\n"
           "\tdiff between host/device = %.14e\n"
           "\ttime/point host = %.14e ns\n"
           "\ttime/point device = %.14e ns\n",
           error_host, diff_hd,
           duration_host.count() / static_cast<Real>(nFineRho * nFineT * NTIMES),
           duration_device.count() / static_cast<Real>(nFineRho * nFineT * NTIMES));

    free(rho_host);
    free(rho_device);
    free(sie_host);
    free(sie_device);
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif

  return 0;
}

#endif // SINGULARITY_USE_EOSPAC
#endif // SPINER_USE_HDF
