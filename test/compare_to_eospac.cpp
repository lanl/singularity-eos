//------------------------------------------------------------------------------
// © 2021-2023. Triad National Security, LLC. All rights reserved.  This
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

/*
  Comparison for SpinerEOS vs. EOSPAC
  Authors: Jonah Miller, Chad Meyers, Josh Dolence, Daniel Holiday
 */

// TODO there's a memory leak/double free somewhere

#ifdef SPINER_USE_HDF
#ifdef SINGULARITY_USE_EOSPAC

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/sp5/singularity_eos_sp5.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/sp5.hpp>
#include <spiner/spiner_types.hpp>

#include <eospac-wrapper/eospac_wrapper.hpp>

#include <ports-of-call/portability.hpp>

#include <singularity-eos/eos/eos.hpp>

using namespace singularity;
using namespace EospacWrapper;

using DataBox = Spiner::DataBox<Real>;
using RegularGrid1D = Spiner::RegularGrid1D<Real>;

using duration = std::chrono::microseconds;
constexpr char diffFileName[] = "diffs.sp5";
constexpr int NTIMES = 10;

class Bounds {

 public:
  Bounds() {}

  Bounds(Real min, Real max, int N, Real offset)
      : grid(RegularGrid1D(min, max, N)), offset(offset) {}

  Bounds(Real min, Real max, int N) : offset(0) {
    constexpr Real epsilon = std::numeric_limits<float>::epsilon();
    const Real min_offset = 10 * std::abs(epsilon);
    if (min <= 0) offset = 1.1 * std::abs(min) + min_offset;
    min += offset;
    max += offset;
    min = std::log10(std::abs(min));
    max = std::log10(std::abs(max));
    grid = RegularGrid1D(min, max, N);
  }

  PORTABLE_INLINE_FUNCTION Real log2lin(Real xl) const { return pow(10., xl) - offset; }
  PORTABLE_INLINE_FUNCTION Real i2lin(int i) const { return log2lin(grid.x(i)); }

  friend std::ostream &operator<<(std::ostream &os, const Bounds &b) {
    os << "Bounds: [" << b.grid.min() << ", " << b.grid.max() << "]"
       << " + " << b.offset << ", "
       << "[N,dx] = [" << b.grid.nPoints() << ", " << b.grid.dx() << "]"
       << "\n";
    return os;
  }

 public:
  RegularGrid1D grid;
  Real offset;
};

int main(int argc, char *argv[]) {

  herr_t status = H5_SUCCESS;

  std::vector<int> matids;

#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize();
#endif
  {

    if (argc < 2) {
      std::cerr << "Usage: " << argv[0]
                << " [-n nFine] [-f filename] matid1 matid2 ... matidN" << std::endl;
      std::exit(1);
    }

#ifdef PORTABILITY_STRATEGY_KOKKOS
    using Kokkos::SpaceAccessibility;
    using KDMS = Kokkos::DefaultExecutionSpace::memory_space;
    using HS = Kokkos::HostSpace;
    constexpr const bool is_host = SpaceAccessibility<HS, KDMS>::accessible;
    constexpr const int def_fine = is_host ? 512 : 2048;
#else
    constexpr const int def_fine = 512;
#endif // PORTABILITY_STRATEGY_KOKKOS

    int nFine = def_fine;
    std::string filename = "materials.sp5";
    for (int i = 1; i < argc; ++i) {
      if (std::strcmp(argv[i], "-f") == 0 && i < argc - 1) {
        filename = argv[i + 1];
        i++;
      } else if (std::strcmp(argv[i], "-n") == 0 && i < argc - 1) {
        nFine = std::atoi(argv[i + 1]);
        i++;
      } else {
        matids.push_back(std::atoi(argv[i]));
      }
    }

    int nFineRho = nFine + 1;
    int nFineT = nFine - 1;
    std::string sp5name = filename;

    std::cout << "Comparison to EOSPAC" << std::endl;

    if (status != H5_SUCCESS) {
      std::cerr << "H5 error. Aborting" << std::endl;
      std::exit(1);
    }

    hid_t file = H5Fcreate(diffFileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    constexpr int NT = 3;
    EOS_INTEGER nXYPairs = (nFineRho) * (nFineT);
    DataBox xVals(nFineRho, nFineT);
    DataBox yVals(nFineRho, nFineT);
    DataBox vars(nFineRho, nFineT);
    DataBox dx(nFineRho, nFineT);
    DataBox dy(nFineRho, nFineT);
    DataBox pressEOSPAC(nFineRho, nFineT);
    DataBox pressDiff_h(nFineRho, nFineT);
    DataBox pressDiff_d(nFineRho, nFineT);

#ifdef PORTABILITY_STRATEGY_KOKKOS
    using RView = Kokkos::View<Real *>;
    RView pressSpinerDView("pressSpinerDevice", nFineRho * nFineT);
    typename RView::HostMirror pressSpinerHView =
        Kokkos::create_mirror_view(pressSpinerDView);
    DataBox pressSpiner_d(pressSpinerDView.data(), nFineRho, nFineT);
    DataBox pressSpiner_hm(pressSpinerHView.data(), nFineRho, nFineT);
#else
    DataBox pressSpiner_d(nFineRho, nFineT);
    DataBox pressSpiner_hm = pressSpiner_d.slice(2, 0, nFineRho);
#endif

    DataBox pressSpiner_h(nFineRho, nFineT);
    DataBox tempSpiner_h(nFineRho, nFineT);
    DataBox tempSpinerE_h(nFineRho, nFineT);
    DataBox tempEOSPAC(nFineRho, nFineT);
    DataBox tempDiff_h(nFineRho, nFineT);
    DataBox tempDiff_d(nFineRho, nFineT);
    DataBox tempDiffE_h(nFineRho, nFineT);
    DataBox tempDiffE_d(nFineRho, nFineT);

#ifdef PORTABILITY_STRATEGY_KOKKOS
    RView tempSpinerDView("tempSpinerDevice", nFineRho * nFineT);
    RView tempSpinerDViewE("tempSpinerSieDevice", nFineRho * nFineT);
    RView::HostMirror tempSpinerHView = Kokkos::create_mirror_view(tempSpinerDView);
    RView::HostMirror tempSpinerHViewE = Kokkos::create_mirror_view(tempSpinerDViewE);
    DataBox tempSpiner_d(tempSpinerDView.data(), nFineRho, nFineT);
    DataBox tempSpiner_hm(tempSpinerHView.data(), nFineRho, nFineT);
    DataBox tempSpinerE_d(tempSpinerDViewE.data(), nFineRho, nFineT);
    DataBox tempSpinerE_hm(tempSpinerHViewE.data(), nFineRho, nFineT);
    RView rhos_v("rhos", nFineRho);
    RView Ts_v("Ts", nFineT);
    RView sies_v("sies", nFineT);
    DataBox rhos(rhos_v.data(), nFineRho);
    DataBox Ts(Ts_v.data(), nFineT);
    DataBox sies(sies_v.data(), nFineT);
#else
    DataBox tempSpiner_d(nFineRho, nFineT);
    DataBox tempSpiner_hm = tempSpiner_d.slice(2, 0, nFineRho);
    DataBox tempSpinerE_d(nFineRho, nFineT);
    DataBox tempSpinerE_hm = tempSpiner_h.slice(2, 0, nFineRho);
    DataBox rhos(nFineRho);
    DataBox Ts(nFineT);
    DataBox sies(nFineT);
#endif

    std::cout << "Comparing to EOSPAC for materials..." << std::endl;
    duration durationEospacTot = duration::zero();
    duration durationSpinerTot = duration::zero();
    duration durationSpinerDTot = duration::zero();
    duration durationEospac = duration::zero();
    duration durationSpiner = duration::zero();
    duration durationSpinerSie = duration::zero();
    duration durationSpinerDev = duration::zero();
    duration durationSpinerSieDev = duration::zero();
    auto start = std::chrono::high_resolution_clock::now();
    auto stop = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < matids.size(); i++) {
      durationEospac = duration::zero();
      durationSpiner = duration::zero();

      int matid = matids[i];
      std::cout << "\t..." << matid << std::endl;
      std::cout << "\t------------------------" << std::endl;

      SesameMetadata metadata;
      eosGetMetadata(matid, metadata, Verbosity::Debug);
      std::cout << metadata << std::endl;

      hid_t idGroup = H5Gcreate(file, std::to_string(matid).c_str(), H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT);

      std::cout << "\t\tloading EOSs\n"
                << "\t\t\tdepends rho T" << std::endl;
      SpinerEOSDependsRhoT eos_host(sp5name, matid);
      std::cout << "\t\t\tdepends rho sie" << std::endl;
      SpinerEOSDependsRhoSie eosE_host(sp5name, matid);
      std::cout << "\t\tmoving to device" << std::endl;
      EOS eos_device = eos_host.GetOnDevice();
      EOS eosE_device = eosE_host.GetOnDevice();

      Real *lambda_dp =
          (Real *)PORTABLE_MALLOC(sizeof(Real) * nFineRho * nFineT * eos_host.nlambda());
      Real *lambda_hp =
          (Real *)malloc(sizeof(Real) * nFineRho * nFineT * eos_host.nlambda());
      DataBox lambda_h(lambda_hp, nFineRho, nFineT, eos_host.nlambda());
      DataBox lambda_d(lambda_dp, nFineRho, nFineT, eos_host.nlambda());

      Real rhoMin = 1.1 * std::max(metadata.rhoMin, 1e-5);
      Real rhoMax = 0.9 * metadata.rhoMax;
      Real TMin = 1.1 * std::max(metadata.TMin, 1.0);
      Real TMax = 0.9 * metadata.TMax;
      Real sieMin = metadata.sieMin + 0.1 * std::abs(metadata.sieMin);
      Real sieMax = 0.9 * metadata.sieMax;

      Bounds lRhoBounds(rhoMin, rhoMax, nFineRho);
      Bounds lTBounds(TMin, TMax, nFineT);
      Bounds leBounds(sieMin, sieMax, nFineT);

      std::cout << "\t\trho bounds = [" << rhoMin << ", " << rhoMax << "]" << std::endl;
      std::cout << "\t\tT   bounds = [" << TMin << ", " << TMax << "]" << std::endl;
      std::cout << "\t\tsie bounds = [" << sieMin << ", " << sieMax << "]" << std::endl;

      EOS_INTEGER tableHandle[NT];
      EOS_INTEGER eospacPofRT, eospacTofRE;
      EOS_INTEGER tableType[NT] = {EOS_Pt_DT, EOS_T_DUt, EOS_Ut_DT};
      eosSafeLoad(NT, matid, tableType, tableHandle, Verbosity::Debug);
      eospacPofRT = tableHandle[0];
      eospacTofRE = tableHandle[1];

      std::cout << "\t\tGenerating interpolation points for rho-T tables" << std::endl;
      start = std::chrono::high_resolution_clock::now();
      for (int j = 0; j < nFineRho; j++) {
        Real rho = lRhoBounds.i2lin(j);
        for (int i = 0; i < nFineT; i++) {
          Real T = lTBounds.i2lin(i);
          xVals(j, i) = rho;
          yVals(j, i) = T;
        }
      }
      stop = std::chrono::high_resolution_clock::now();
      durationEospac += std::chrono::duration_cast<duration>(stop - start);

      portableFor(
          "densities", 0, nFineRho,
          PORTABLE_LAMBDA(const int &j) { rhos(j) = lRhoBounds.i2lin(j); });
      portableFor(
          "temperatures", 0, nFineT, PORTABLE_LAMBDA(const int &i) {
            Ts(i) = lTBounds.i2lin(i);
            sies(i) = leBounds.i2lin(i);
          });

      std::cout << "\t\tInterpolate pressure of rho T" << std::endl;
      {
        pressSpiner_h.setRange(0, lTBounds.grid);
        pressSpiner_h.setRange(1, lRhoBounds.grid);
        pressSpiner_d.setRange(0, lTBounds.grid);
        pressSpiner_d.setRange(1, lRhoBounds.grid);
        pressSpiner_hm.setRange(0, lTBounds.grid);
        pressSpiner_hm.setRange(1, lRhoBounds.grid);
        pressEOSPAC.copyMetadata(pressSpiner_h);
        pressDiff_h.copyMetadata(pressSpiner_h);
        pressDiff_d.copyMetadata(pressSpiner_h);

        std::cout << "\t\t...eospac..." << std::endl;
        start = std::chrono::high_resolution_clock::now();
        for (int n = 0; n < NTIMES; n++) {
          for (int i = 0; i < nXYPairs; i++)
            vars(i) = dx(i) = dy(i) = 0;
          eosSafeInterpolate(&eospacPofRT, nXYPairs, xVals.data(), yVals.data(),
                             vars.data(), dx.data(), dy.data(), "PofRT",
                             Verbosity::Quiet);
        }
        stop = std::chrono::high_resolution_clock::now();
        durationEospac += std::chrono::duration_cast<duration>(stop - start);

        std::cout << "\t\t...spiner on host..." << std::endl;
        start = std::chrono::high_resolution_clock::now();
        for (int n = 0; n < NTIMES; n++) {
#pragma omp simd
          for (int j = 0; j < nFineRho; j++) {
            const Real rho = xVals(j, 0);
            for (int i = 0; i < nFineT; i++) {
              const Real T = yVals(0, i);
              pressSpiner_h(j, i) = eos_host.PressureFromDensityTemperature(rho, T);
            }
          }
        }
        stop = std::chrono::high_resolution_clock::now();
        durationSpiner += std::chrono::duration_cast<duration>(stop - start);

        std::cout << "\t\t...spiner on device..." << std::endl;
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
#endif
        start = std::chrono::high_resolution_clock::now();
        for (int n = 0; n < NTIMES; ++n) {
          portableFor(
              "pressure from density and temperature", 0, nFineRho, 0, nFineT,
              PORTABLE_LAMBDA(const int &j, const int &i) {
                Real rho = rhos(j);
                Real T = Ts(i);
                pressSpiner_d(j, i) = eos_device.PressureFromDensityTemperature(rho, T);
              });
        }
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
#endif
        stop = std::chrono::high_resolution_clock::now();
        durationSpinerDev += std::chrono::duration_cast<duration>(stop - start);

        Real diffL2{0.0}, L2{0.0};
        std::cout << "\t\t...comparing on host..." << std::endl;
        for (int j = 0; j < nFineRho; j++) {
          for (int i = 0; i < nFineT; i++) {
            const Real p_true = pressureFromSesame(vars(j, i));
            pressEOSPAC(j, i) = p_true;
            const Real diff = pressSpiner_h(j, i) - p_true;
            const Real mean_p = 0.5 * diff + p_true;
            pressDiff_h(j, i) = diff;
            const Real reldiff = diff / (1e-10 + std::abs(mean_p));
            diffL2 += reldiff * reldiff;
            L2 += 1;
          }
        }
        diffL2 /= L2;
        diffL2 = std::sqrt(diffL2);
        std::cout << "\t\t...L2 difference = " << diffL2 << std::endl;

        hid_t pressGroup =
            H5Gcreate(idGroup, "pressuresRhoT", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status += pressSpiner_h.saveHDF(pressGroup, "spiner_h");
        status += pressEOSPAC.saveHDF(pressGroup, "EOSPAC");
        status += pressDiff_h.saveHDF(pressGroup, "diff_h");

#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::deep_copy(pressSpinerHView, pressSpinerDView);
#endif
        std::cout << "\t\t...comparing device results..." << std::endl;
        diffL2 = 0;
        L2 = 0;
        for (int j = 0; j < nFineRho; j++) {
          for (int i = 0; i < nFineT; i++) {
            // pressSpiner_d has been copied into pressSpiner_h
            const Real p_true = pressEOSPAC(j, i);
            const Real diff = pressSpiner_hm(j, i) - p_true;
            const Real mean_p = 0.5 * diff + p_true;
            pressDiff_d(j, i) = diff;
            const Real reldiff = diff / (1e-10 + std::abs(mean_p));
            diffL2 += reldiff * reldiff;
            L2 += 1;
          }
        }
        diffL2 /= L2;
        diffL2 = std::sqrt(diffL2);
        std::cout << "\t\t...L2 difference = " << diffL2 << std::endl;

        status += pressSpiner_hm.saveHDF(pressGroup, "spiner_d");
        status += pressDiff_d.saveHDF(pressGroup, "diff_d");
        status += H5Gclose(pressGroup);

        std::cout << "\t\t...EOSPAC time/point (microseconds) = "
                  << durationEospac.count() /
                         static_cast<Real>((nFineRho) * (nFineT)*NTIMES)
                  << std::endl;
        std::cout << "\t\t...spiner host time/point (microseconds) = "
                  << durationSpiner.count() /
                         static_cast<Real>((nFineRho) * (nFineT)*NTIMES)
                  << std::endl;
        std::cout << "\t\t...spiner device time/point (microseconds) = "
                  << durationSpinerDev.count() /
                         static_cast<Real>((nFineRho) * (nFineT)*NTIMES)
                  << std::endl;
        std::cout << "\t\t...saving to file..." << std::endl;
      }
      durationEospacTot += durationEospac;
      durationSpinerTot += durationSpiner;
      durationSpinerDTot += durationSpinerDev;
      durationEospac = duration::zero();
      durationSpiner = duration::zero();
      durationSpinerDev = duration::zero();

      std::cout << "\t\tGenerating interpolation points for rho-sie tables" << std::endl;
      start = std::chrono::high_resolution_clock::now();
      for (int j = 0; j < nFineRho; j++) {
        Real rho = lRhoBounds.i2lin(j);
        for (int i = 0; i < nFineT; i++) {
          Real sie = leBounds.i2lin(i);
          xVals(j, i) = densityToSesame(rho);
          yVals(j, i) = sieToSesame(sie);
        }
      }
      stop = std::chrono::high_resolution_clock::now();
      durationEospac += std::chrono::duration_cast<duration>(stop - start);
      std::cout << "\t\tRoot find temperature of rho sie" << std::endl;
      {
        tempSpiner_h.setRange(0, leBounds.grid);
        tempSpiner_h.setRange(1, lRhoBounds.grid);
        tempSpiner_d.setRange(0, leBounds.grid);
        tempSpiner_d.setRange(1, lRhoBounds.grid);
        tempSpiner_hm.setRange(0, leBounds.grid);
        tempSpiner_hm.setRange(1, lRhoBounds.grid);
        tempSpinerE_h.setRange(0, leBounds.grid);
        tempSpinerE_h.setRange(1, lRhoBounds.grid);
        tempSpinerE_d.setRange(0, leBounds.grid);
        tempSpinerE_d.setRange(1, lRhoBounds.grid);
        tempSpinerE_hm.setRange(0, leBounds.grid);
        tempSpinerE_hm.setRange(1, lRhoBounds.grid);
        tempEOSPAC.copyMetadata(tempSpiner_h);
        tempDiff_h.copyMetadata(tempSpiner_h);
        tempDiff_d.copyMetadata(tempSpiner_h);
        tempDiffE_h.copyMetadata(tempSpiner_h);
        tempDiffE_d.copyMetadata(tempSpiner_d);

        std::cout << "\t\t...eospac..." << std::endl;
        start = std::chrono::high_resolution_clock::now();
        for (int n = 0; n < NTIMES; n++) {
          eosSafeInterpolate(&eospacTofRE, nXYPairs, xVals.data(), yVals.data(),
                             vars.data(), dx.data(), dy.data(), "TofRE",
                             Verbosity::Quiet);
        }
        stop = std::chrono::high_resolution_clock::now();
        durationEospac += std::chrono::duration_cast<duration>(stop - start);

        std::cout << "\t\t...spiner(rho,T) on host..." << std::endl;
        start = std::chrono::high_resolution_clock::now();
        for (int n = 0; n < NTIMES; n++) {
          for (int j = 0; j < nFineRho; j++) {
            Real rho = lRhoBounds.i2lin(j);
            DataBox lambda_hj = lambda_h.slice(j);
            for (int i = 0; i < nFineT; i++) {
              Real sie = leBounds.i2lin(i);
              DataBox lambda = lambda_hj.slice(i);
              tempSpiner_h(j, i) =
                  eos_host.TemperatureFromDensityInternalEnergy(rho, sie, lambda.data());
            }
          }
        }
        stop = std::chrono::high_resolution_clock::now();
        durationSpiner += std::chrono::duration_cast<duration>(stop - start);

        std::cout << "\t\t...spiner(rho,T) on device..." << std::endl;
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
#endif
        start = std::chrono::high_resolution_clock::now();
        for (int n = 0; n < NTIMES; ++n) {
          portableFor(
              "T(rho,e)", 0, nFineRho, 0, nFineT,
              PORTABLE_LAMBDA(const int &j, const int &i) {
                Real rho = rhos(j);
                Real sie = sies(i);
                DataBox lambda = lambda_d.slice(j).slice(i);
                tempSpiner_d(j, i) = eos_device.TemperatureFromDensityInternalEnergy(
                    rho, sie, lambda.data());
              });
        }
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
#endif
        stop = std::chrono::high_resolution_clock::now();
        durationSpinerDev += std::chrono::duration_cast<duration>(stop - start);

        std::cout << "\t\t...spiner(rho,sie) on host..." << std::endl;
        start = std::chrono::high_resolution_clock::now();
        for (int n = 0; n < NTIMES; n++) {
          for (int j = 0; j < nFineRho; j++) {
            Real rho = lRhoBounds.i2lin(j);
            for (int i = 0; i < nFineT; i++) {
              Real sie = leBounds.i2lin(i);
              tempSpinerE_h(j, i) =
                  eosE_host.TemperatureFromDensityInternalEnergy(rho, sie);
            }
          }
        }
        stop = std::chrono::high_resolution_clock::now();
        durationSpinerSie += std::chrono::duration_cast<duration>(stop - start);

        std::cout << "\t\t...spiner(rho,sie) on device..." << std::endl;
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
#endif
        start = std::chrono::high_resolution_clock::now();
        for (int n = 0; n < NTIMES; ++n) {
          portableFor(
              "T(rho,e)", 0, nFineRho, 0, nFineT,
              PORTABLE_LAMBDA(const int &j, const int &i) {
                Real rho = rhos(j);
                Real sie = sies(i);
                tempSpinerE_d(j, i) =
                    eosE_device.TemperatureFromDensityInternalEnergy(rho, sie);
              });
        }
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
#endif
        stop = std::chrono::high_resolution_clock::now();
        durationSpinerSieDev += std::chrono::duration_cast<duration>(stop - start);

        std::cout << "\t\t...comparing host data..." << std::endl;
        Real diffL2, diffL2E, diffL2_d, diffL2E_d;
        Real L2{0.0};
        diffL2 = 0;
        for (int j = 0; j < nFineRho; j++) {
          for (int i = 0; i < nFineT; i++) {
            const Real t_true = temperatureFromSesame(vars(j, i));
            tempEOSPAC(j, i) = t_true;
            Real diff = tempSpiner_h(j, i) - t_true;
            tempDiff_h(j, i) = diff;
            const Real mean_t = 0.5 * diff + t_true;
            if (std::isnan(diff)) {
              std::cout << "NAN! " << j << ", " << i << ", " << tempSpiner_h(j, i) << ", "
                        << tempEOSPAC(j, i) << ", " << mean_t << ", " << diff
                        << std::endl;
              exit(1);
            }
#ifdef singularity_normalize_cold_point
            if (tempSpiner_h(j, i) <= 2. * TMin && tempEOSPAC(j, i) < TMin) {
              tempSpiner_h(j, i) = 0;
              tempEOSPAC(j, i) = 0;
              diff = 0;
            }
#endif
            Real reldiff = diff / (1e-10 + std::abs(mean_t));
            diffL2 += reldiff * reldiff;
            L2 += 1;
          }
        }
        diffL2 /= L2;
        diffL2 = std::sqrt(diffL2);
        diffL2E = 0;
        L2 = 0;
        for (int j = 0; j < nFineRho; j++) {
          for (int i = 0; i < nFineT; i++) {
            const Real t_true = tempEOSPAC(j, i);
            const Real diff = tempSpinerE_h(j, i) - t_true;
            tempDiffE_h(j, i) = diff;
            const Real mean_t = 0.5 * diff + t_true;
            const Real reldiff = diff / (1e-10 + std::abs(mean_t));
            diffL2E += reldiff * reldiff;
            L2 += 1;
          }
        }
        diffL2E /= L2;
        diffL2E = std::sqrt(diffL2E);
        std::cout << "\t\tRoot finding:\n"
                  << "\t\tits: counts\n";
        for (int i = 0; i < eos_host.counts.nBins(); i++) {
          std::cout << "\t\t\t" << i << ": "
                    << 100. * eos_host.counts[i] / eos_host.counts.total() << "\n";
        }
        hid_t pressGroup = H5Gcreate(idGroup, "temperaturesRhoSie", H5P_DEFAULT,
                                     H5P_DEFAULT, H5P_DEFAULT);
        status += tempSpiner_h.saveHDF(pressGroup, "spiner_h");
        status += tempEOSPAC.saveHDF(pressGroup, "EOSPAC");
        status += tempDiff_h.saveHDF(pressGroup, "diff_h");
        status += tempSpinerE_h.saveHDF(pressGroup, "spinerE_h");
        status += tempDiffE_h.saveHDF(pressGroup, "diffE_h");

#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::deep_copy(tempSpinerHView, tempSpinerDView);
        Kokkos::deep_copy(tempSpinerHViewE, tempSpinerDViewE);
#endif
        // tempSpiner_hm now contains copy of device data
        diffL2_d = 0;
        L2 = 0;
        for (int j = 0; j < nFineRho; j++) {
          for (int i = 0; i < nFineT; i++) {
            const Real t_true = tempEOSPAC(j, i);
            Real diff = tempSpiner_hm(j, i) - t_true;
            tempDiff_d(j, i) = diff;
            const Real mean_t = 0.5 * diff + t_true;
#ifdef singularity_normalize_cold_point
            if (tempSpiner_hm(j, i) < 2 * TMin && tempEOSPAC(j, i) < TMin) {
              tempSpiner_hm(j, i) = 0;
              tempEOSPAC(j, i) = 0;
              diff = 0;
            }
#endif
            Real reldiff = diff / (1e-10 + std::abs(mean_t));
            diffL2_d += reldiff * reldiff;
            L2 += 1;
          }
        }
        diffL2_d /= L2;
        diffL2_d = std::sqrt(diffL2_d);
        diffL2E_d = 0;
        L2 = 0;
        for (int j = 0; j < nFineRho; j++) {
          for (int i = 0; i < nFineT; i++) {
            const Real t_true = tempEOSPAC(j, i);
            const Real diff = tempSpinerE_hm(j, i) - t_true;
            tempDiffE_d(j, i) = diff;
            const Real mean_t = 0.5 * diff + t_true;
            Real reldiff = diff / (1e-10 + std::abs(mean_t));
            diffL2E_d += reldiff * reldiff;
            L2 += 1;
          }
        }
        diffL2E_d /= L2;
        diffL2E_d = std::sqrt(diffL2E_d);
        status += tempSpiner_hm.saveHDF(pressGroup, "spiner_d");
        status += tempDiff_d.saveHDF(pressGroup, "diff_d");
        status += tempSpinerE_hm.saveHDF(pressGroup, "spinerE_d");
        status += tempDiffE_d.saveHDF(pressGroup, "diffE_d");
        status += H5Gclose(pressGroup);

        std::cout << "\t------------------------" << std::endl;
        std::cout << "\t\t...L2 difference for spiner(rho,T) host = " << diffL2
                  << std::endl;
        std::cout << "\t\t...L2 difference for spiner(rho,sie) host = " << diffL2E
                  << std::endl;
        std::cout << "\t\t...L2 difference for spiner(rho,T) device = " << diffL2_d
                  << std::endl;
        std::cout << "\t\t...L2 difference for spiner(rho,sie) device = " << diffL2E_d
                  << std::endl;
        std::cout << "\t\t...EOSPAC time/point (microseconds) = "
                  << durationEospac.count() /
                         static_cast<Real>((nFineRho) * (nFineT)*NTIMES)
                  << std::endl;
        std::cout << "\t\t...spiner(rho,T) time/point (microseconds) host = "
                  << durationSpiner.count() /
                         static_cast<Real>((nFineRho) * (nFineT)*NTIMES)
                  << std::endl;
        std::cout << "\t\t...spiner(rho,sie) time/point (microseconds) host = "
                  << durationSpinerSie.count() /
                         static_cast<Real>((nFineRho) * (nFineT)*NTIMES)
                  << std::endl;
        std::cout << "\t\t...spiner(rho,T) time/point (microseconds) device = "
                  << durationSpinerDev.count() /
                         static_cast<Real>((nFineRho) * (nFineT)*NTIMES)
                  << std::endl;
        std::cout << "\t\t...spiner(rho,sie) time/point (microseconds) device = "
                  << durationSpinerSieDev.count() /
                         static_cast<Real>((nFineRho) * (nFineT)*NTIMES)
                  << std::endl;
      }
      std::cout << "\t------------------------" << std::endl;

      status += H5Gclose(idGroup);

      PORTABLE_FREE(lambda_dp);
      free(lambda_hp);
      eos_device.Finalize();
      eosE_device.Finalize();
    }

    status += H5Fclose(file);
    if (status != H5_SUCCESS) {
      std::cerr << "WARNING: problem with HDf5" << std::endl;
    }
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif

  return 0;
}

#endif /// SINGULARITY_USE_EOSPAC
#endif // SPINER_USE_HDF
