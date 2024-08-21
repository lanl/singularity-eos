//------------------------------------------------------------------------------
// © 2021-2024. Triad National Security, LLC. All rights reserved.  This
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

#ifndef SINGULARITY_EOS_BASE_TABLE_UTILS_HPP_
#define SINGULARITY_EOS_BASE_TABLE_UTILS_HPP_
#ifdef SINGULARITY_USE_SPINER

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/fast-math/logs.hpp>

#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>

namespace singularity {
namespace table_utils {
// For logarithmic interpolation, quantities may be negative.
// If they are, use offset to ensure negative values make sense.
template <int NGRIDS = 3>
class Bounds {
 public:
  using RegularGrid1D = Spiner::RegularGrid1D<Real>;
  using Grid_t = Spiner::PiecewiseGrid1D<Real, NGRIDS>;
  // for tag dispatch in constructors
  using OneGrid = std::integral_constant<int, 1>;
  using TwoGrids = std::integral_constant<int, 2>;
  using ThreeGrids = std::integral_constant<int, 3>;

  Bounds() : offset(0), piecewise(false), linmin_(0), linmax_(0) {}

  Bounds(const Real min, const Real max, const int N, const Real offset)
      : grid(Grid_t(std::vector<RegularGrid1D>{RegularGrid1D(min, max, N)})),
        offset(offset), piecewise(false), linmin_(min), linmax_(max) {}
  Bounds(OneGrid, const Real min, const Real max, const int N, const Real offset)
      : Bounds(min, max, N, offset) {}

  Bounds(Real min, Real max, int N, const bool convertToLog = false,
         const Real shrinkRange = 0,
         Real anchor_point = std::numeric_limits<Real>::signaling_NaN())
      : offset(0), piecewise(true), linmin_(min), linmax_(max) {
    if (convertToLog) {
      convertBoundsToLog_(min, max, shrinkRange);
      if (!(std::isnan(anchor_point))) {
        anchor_point = singularity::FastMath::log10(std::abs(anchor_point));
      }
    }
    if (!(std::isnan(anchor_point))) {
      adjustForAnchor_(min, max, N, anchor_point);
    }
    grid = Grid_t(std::vector<RegularGrid1D>{RegularGrid1D(min, max, N)});
  }
  Bounds(OneGrid, Real min, Real max, int N, const bool convertToLog = false,
         const Real shrinkRange = 0,
         Real anchor_point = std::numeric_limits<Real>::signaling_NaN())
      : Bounds(min, max, N, convertToLog, shrinkRange, anchor_point) {}

  Bounds(TwoGrids, Real global_min, Real global_max, Real anchor_point, Real splitPoint,
         const Real ppd_fine, const Real ppd_factor, const bool convertToLog,
         const Real shrinkRange = 0)
      : offset(0), piecewise(true), linmin_(global_min), linmax_(global_max) {
    const Real ppd_coarse = (ppd_factor > 0) ? ppd_fine / ppd_factor : ppd_fine;

    if (convertToLog) {
      convertBoundsToLog_(global_min, global_max, shrinkRange);
      anchor_point = toLog_(anchor_point, offset);
      splitPoint = toLog_(splitPoint, offset);
    }

    checkInterval_(splitPoint, global_min, global_max, "Split point");
    checkInterval_(anchor_point, global_min, splitPoint, "Anchor point");

    // add a point just to make sure we have enough points after adjusting for anchor
    int N_fine = getNumPointsFromDensity(global_min, splitPoint, ppd_fine) + 1;
    adjustForAnchor_(global_min, splitPoint, N_fine, anchor_point);
    RegularGrid1D grid_lower(global_min, splitPoint, N_fine);

    const int N_upper = getNumPointsFromDensity(splitPoint, global_max, ppd_coarse);
    RegularGrid1D grid_upper(splitPoint, global_max, N_upper);

    grid = Grid_t(std::vector<RegularGrid1D>{grid_lower, grid_upper});
  }

  Bounds(ThreeGrids, Real global_min, Real global_max, Real anchor_point, Real fine_min,
         Real fine_max, const Real ppd_fine, const Real ppd_factor_lo,
         const Real ppd_factor_hi, const bool convertToLog, const Real shrinkRange = 0)
      : offset(0), piecewise(true), linmin_(global_min), linmax_(global_max) {

    if (convertToLog) {
      convertBoundsToLog_(global_min, global_max, shrinkRange);
      anchor_point = toLog_(anchor_point, offset);
      fine_min = toLog_(fine_min, offset);
      fine_max = toLog_(fine_max, offset);
    }
    checkInterval_(anchor_point, global_min, global_max, "Anchor point");

    grid = gridFromIntervals_(ThreeGrids(), global_min, global_max, anchor_point,
                              fine_min, fine_max, ppd_fine, ppd_factor_lo, ppd_factor_hi);
  }

  Bounds(ThreeGrids, Real global_min, Real global_max, Real anchor_point,
         Real log_fine_diameter, const Real ppd_fine, const Real ppd_factor_lo,
         const Real ppd_factor_hi, const bool convertToLog, const Real shrinkRange = 0)
      : offset(0), piecewise(true), linmin_(global_min), linmax_(global_max) {

    if (convertToLog) {
      convertBoundsToLog_(global_min, global_max, shrinkRange);
      anchor_point = toLog_(anchor_point, offset);
    }

    checkInterval_(anchor_point, global_min, global_max, "Anchor point");
    Real mid_min = anchor_point - 0.5 * log_fine_diameter;
    Real mid_max = anchor_point + 0.5 * log_fine_diameter;

    grid = gridFromIntervals_(ThreeGrids(), global_min, global_max, anchor_point, mid_min,
                              mid_max, ppd_fine, ppd_factor_lo, ppd_factor_hi);
  }

  inline Real log2lin(const Real xl) const {
    // JMM: Need to guard this with the linear bounds passed in. The
    // reason is that our fast math routines, while completely
    // invertible at the level of machine epsilon, do introduce error
    // at machine epsilon, which can bring us out of the interpolation
    // range of eospac.
    return std::min(linmax_,
                    std::max(linmin_, singularity::FastMath::pow10(xl) - offset));
  }
  inline Real i2lin(const int i) const { return log2lin(grid.x(i)); }

  friend std::ostream &operator<<(std::ostream &os, const Bounds &b) {
    os << "Bounds: [" << b.grid.min() << ", " << b.grid.max() << "]"
       << " + " << b.offset << "\n"
       << "\tN = " << b.grid.nPoints() << "\n";
    for (int ig = 0; ig < b.grid.nGrids(); ++ig) {
      os << "\t[ig,dx] = [" << ig << ", " << b.grid.dx(ig) << "]"
         << "\n";
    }
    return os;
  }

  // This uses real logs
  template <typename T>
  static int getNumPointsFromPPD(Real min, Real max, const T ppd) {
    constexpr Real epsilon = std::numeric_limits<Real>::epsilon();
    const Real min_offset = 10 * std::abs(epsilon);
    // 1.1 so that the y-intercept isn't identically zero
    // when min < 0.
    // min_offset to handle the case where min=0
    Real offset = 0;
    if (min <= 0) offset = 1.1 * std::abs(min) + min_offset;
    min += offset;
    max += offset;

    Real lmin = std::log10(min);
    Real lmax = std::log10(max);
    int N = getNumPointsFromDensity(lmin, lmax, ppd);
    return N;
  }
  template <typename T>
  static int getNumPointsFromDensity(const Real min, const Real max, const T density) {
    Real delta = max - min;
    int N = std::max(2, static_cast<int>(std::ceil(density * delta)));
    return N;
  }

 private:
  Real toLog_(Real val, const Real offset) {
    val += offset;
    return singularity::FastMath::log10(std::abs(val));
  }

  void convertBoundsToLog_(Real &min, Real &max, const Real shrinkRange = 0) {
    // Log scales can't handle negative numbers or exactly zero. To
    // deal with that, we offset.
    constexpr Real epsilon = std::numeric_limits<Real>::epsilon();
    const Real min_offset = 10 * std::abs(epsilon);
    // 1.1 so that the y-intercept isn't identically zero
    // when min < 0.
    // min_offset to handle the case where min=0
    if (min <= 0) offset = 1.1 * std::abs(min) + min_offset;

    min = toLog_(min, offset);
    max = toLog_(max, offset);

    Real delta = max - min;
    min += 0.5 * shrinkRange * delta;
    max -= 0.5 * shrinkRange * delta;
  }

  static void adjustForAnchor_(const Real min, Real &max, int &N,
                               const Real anchor_point) {
    if (min < anchor_point && anchor_point < max) {
      Real Nfrac = (anchor_point - min) / (max - min);
      PORTABLE_REQUIRE((0 < Nfrac && Nfrac < 1), "anchor in bounds");

      int Nanchor = static_cast<int>(std::ceil(N * Nfrac));
      if ((Nanchor < 2) || (Nanchor >= N)) return; // not possible to shift this safely

      Real dx = (anchor_point - min) / static_cast<Real>(Nanchor - 1);
      int Nmax_new = static_cast<int>((max - min) / dx);
      Real max_new = dx * (Nmax_new - 1) + min;
      PORTABLE_REQUIRE(max_new <= max, "must not exceed table bounds");

      N = Nmax_new;
      max = max_new;
    } else {
      PORTABLE_ALWAYS_WARN("Anchor point out of bounds. Ignoring it.");
    }
  }

  static void checkInterval_(Real &p, const Real min, const Real max,
                             const std::string &name) {
    if (p <= min) {
      PORTABLE_ALWAYS_WARN(name + " less than minimum. Adjusting.");
      Real eps = 0.1 * std::abs(min);
      p = min + eps;
    }
    if (p >= max) {
      PORTABLE_ALWAYS_WARN(name + " greater than maximum. Adjusting.");
      Real eps = 0.1 * std::abs(max);
      p = max - eps;
    }
  }

  static Grid_t gridFromIntervals_(ThreeGrids, Real global_min, Real global_max,
                                   Real anchor_point, Real mid_min, Real mid_max,
                                   const Real ppd_fine, const Real ppd_factor_lo,
                                   const Real ppd_factor_hi) {
    const Real ppd_lo = (ppd_factor_lo > 0) ? ppd_fine / ppd_factor_lo : ppd_fine;
    const Real ppd_hi = (ppd_factor_hi > 0) ? ppd_fine / ppd_factor_hi : ppd_fine;

    if (mid_min <= global_min) {
      PORTABLE_ALWAYS_WARN(
          "Table bounds refined minimum lower than global minimum. Adjusting.");
      Real delta = std::abs(anchor_point - global_min);
      mid_min = anchor_point - 0.9 * delta;
    }
    if (mid_max >= global_max) {
      PORTABLE_ALWAYS_WARN(
          "Table bounds refined maximum greater than global maximum. Adjusting.");
      Real delta = std::abs(global_max - anchor_point);
      mid_max = anchor_point + 0.9 * delta;
    }

    // add a point just to make sure we have enough points after adjusting for anchor
    int N_fine = getNumPointsFromDensity(mid_min, mid_max, ppd_fine) + 1;
    adjustForAnchor_(mid_min, mid_max, N_fine, anchor_point);
    RegularGrid1D grid_middle(mid_min, mid_max, N_fine);

    const int N_lower = getNumPointsFromDensity(global_min, mid_min, ppd_lo);
    RegularGrid1D grid_lower(global_min, mid_min, N_lower);

    const int N_upper = getNumPointsFromDensity(mid_max, global_max, ppd_hi);
    RegularGrid1D grid_upper(mid_max, global_max, N_upper);

    return Grid_t(std::vector<RegularGrid1D>{grid_lower, grid_middle, grid_upper});
  }

 public:
  Grid_t grid;
  Real offset = 0;
  bool piecewise = false;

 private:
  Real linmin_, linmax_;
};

// JMM: Making this a struct with static methods, rather than a
// namespace, saves a few "friend" declarations
template <typename EOS>
struct SpinerTricks {
  static auto GetOnDevice(EOS *peos_h) {
    // trivially copy all but dynamic memory
    EOS eos_d = *peos_h;
    auto pdbs_d = eos_d.GetDataBoxPointers_();
    auto pdbs_h = peos_h->GetDataBoxPointers_();
    int idb = 0;
    for (auto *pdb_d : pdbs_d) {
      auto *pdb_h = pdbs_h[idb++];
      *pdb_d = pdb_h->getOnDevice();
    }
    // set memory status
    eos_d.memoryStatus_ = DataStatus::OnDevice;
    return eos_d;
  }
  static void Finalize(EOS *peos) {
    if (peos->memoryStatus_ != DataStatus::UnManaged) {
      for (auto *pdb : peos->GetDataBoxPointers_()) {
        pdb->finalize();
      }
    }
    peos->memoryStatus_ = DataStatus::Deallocated;
  }
  static std::size_t DynamicMemorySizeInBytes(const EOS *peos) {
    std::size_t out = 0;
    for (const auto *pdb : peos->GetDataBoxPointers_()) {
      out += pdb->sizeBytes();
    }
    return out;
  }
  static std::size_t DumpDynamicMemory(char *dst, const EOS *peos) {
    std::size_t offst = 0;
    for (const auto *pdb : peos->GetDataBoxPointers_()) {
      std::size_t size = pdb->sizeBytes();
      memcpy(dst + offst, pdb->data(), size);
      offst += size;
    }
    return offst;
  }
  static std::size_t SetDynamicMemory(char *src, EOS *peos) {
    std::size_t offst = 0;
    for (auto *pdb : peos->GetDataBoxPointers_()) {
      offst += pdb->setPointer(src + offst);
    }
    peos->memoryStatus_ = DataStatus::UnManaged;
    return offst;
  }
};
} // namespace table_utils
} // namespace singularity

#endif // SINGULARITY_USE_SPINER
#endif // SINGULARITY_EOS_BASE_TABLE_UTILS_HPP_
