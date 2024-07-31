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

#ifndef SINGULARITY_EOS_BASE_TABLE_BOUNDS_HPP_
#define SINGULARITY_EOS_BASE_TABLE_BOUNDS_HPP_
#ifdef SINGULARITY_USE_SPINER

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/fast-math/logs.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>

namespace singularity {

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

  Bounds() {}

  Bounds(Real min, Real max, int N, Real offset)
      : grid(Grid_t(std::vector<RegularGrid1D>{RegularGrid1D(min, max, N)})),
        offset(offset), piecewise(false) {}
  Bounds(OneGrid, Real min, Real max, int N, Real offset) : Bounds(min, max, N, offset) {}

  Bounds(Real min, Real max, int N, bool convertToLog = false, Real shrinkRange = 0,
         Real anchor_point = std::numeric_limits<Real>::signaling_NaN())
      : offset(0), piecewise(true) {
    if (convertToLog) {
      convertBoundsToLog_(min, max, shrinkRange);
      if (!(std::isnan(anchor_point))) {
        anchor_point += offset;
        anchor_point = singularity::FastMath::log10(std::abs(anchor_point));
      }
    }
    if (!(std::isnan(anchor_point))) {
      adjustForAnchor_(min, max, N, anchor_point);
    }
    grid = Grid_t(std::vector<RegularGrid1D>{RegularGrid1D(min, max, N)});
  }
  Bounds(OneGrid, Real min, Real max, int N, bool convertToLog = false,
         Real shrinkRange = 0,
         Real anchor_point = std::numeric_limits<Real>::signaling_NaN())
      : Bounds(min, max, N, convertToLog, shrinkRange, anchor_point) {}

  Bounds(TwoGrids, Real global_min, Real global_max, Real anchor_point, Real splitPoint,
         Real ppd_fine, Real ppd_factor, Real shrinkRange = 0)
      : offset(0), piecewise(true) {
    const Real ppd_coarse = (ppd_factor > 0) ? ppd_fine / ppd_factor : ppd_fine;

    convertBoundsToLog_(global_min, global_max, shrinkRange);
    anchor_point += offset;
    anchor_point = singularity::FastMath::log10(std::abs(anchor_point));
    splitPoint += offset;
    splitPoint = singularity::FastMath::log10(std::abs(splitPoint));

    // add a point just to make sure we have enough points after adjusting for anchor
    int N_fine = getNumPointsFromPPD(global_min, splitPoint, ppd_fine) + 1;
    adjustForAnchor_(global_min, splitPoint, N_fine, anchor_point);
    RegularGrid1D grid_lower(global_min, splitPoint, N_fine);

    const int N_upper = getNumPointsFromPPD(splitPoint, global_max, ppd_coarse);
    RegularGrid1D grid_upper(splitPoint, global_max, N_upper);

    grid = Grid_t(std::vector<RegularGrid1D>{grid_lower, grid_upper});
  }

  Bounds(ThreeGrids, Real global_min, Real global_max, Real anchor_point,
         Real log_fine_diameter, Real ppd_fine, Real ppd_factor, Real shrinkRange = 0) {
    const Real ppd_coarse = (ppd_factor > 0) ? ppd_fine / ppd_factor : ppd_fine;

    convertBoundsToLog_(global_min, global_max, shrinkRange);
    anchor_point += offset;
    anchor_point = singularity::FastMath::log10(std::abs(anchor_point));

    Real mid_min = anchor_point - 0.5 * log_fine_diameter;
    Real mid_max = anchor_point + 0.5 * log_fine_diameter;

    // add a point just to make sure we have enough points after adjusting for anchor
    int N_fine = getNumPointsFromPPD(mid_min, mid_max, ppd_fine) + 1;
    adjustForAnchor_(mid_min, mid_max, N_fine, anchor_point);
    RegularGrid1D grid_middle(mid_min, mid_max, N_fine);

    const int N_lower = getNumPointsFromPPD(global_min, mid_min, ppd_coarse);
    RegularGrid1D grid_lower(global_min, mid_min, N_lower);

    const int N_upper = getNumPointsFromPPD(mid_max, global_max, ppd_coarse);
    RegularGrid1D grid_upper(mid_max, global_max, N_upper);

    grid = Grid_t(std::vector<RegularGrid1D>{grid_lower, grid_middle, grid_upper});
  }

  inline Real log2lin(Real xl) const { return singularity::FastMath::pow10(xl) - offset; }
  inline Real i2lin(int i) const { return log2lin(grid.x(i)); }

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
  static int getNumPointsFromPPD(Real min, Real max, T ppd) {
    constexpr Real epsilon = std::numeric_limits<float>::epsilon();
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
    Real ndecades = lmax - lmin;
    return std::max(2, static_cast<int>(std::ceil(ppd * ndecades)));
  }

 private:
  void convertBoundsToLog_(Real &min, Real &max, Real shrinkRange = 0) {
    // Log scales can't handle negative numbers or exactly zero. To
    // deal with that, we offset.
    constexpr Real epsilon = std::numeric_limits<float>::epsilon();
    const Real min_offset = 10 * std::abs(epsilon);
    // 1.1 so that the y-intercept isn't identically zero
    // when min < 0.
    // min_offset to handle the case where min=0
    if (min <= 0) offset = 1.1 * std::abs(min) + min_offset;

    min += offset;
    max += offset;

    min = singularity::FastMath::log10(std::abs(min));
    max = singularity::FastMath::log10(std::abs(max));
    Real delta = max - min;
    min += 0.5 * shrinkRange * delta;
    max -= 0.5 * shrinkRange * delta;
  }

  void adjustForAnchor_(Real min, Real &max, int &N, Real anchor_point) {
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
    }
  }

 public:
  Grid_t grid;
  Real offset = 0;
  bool piecewise = false;
};

} // namespace singularity

#endif // SINGULARITY_USE_SPINER
#endif // SINGULARITY_EOS_BASE_TABLE_BOUNDS_HPP_
