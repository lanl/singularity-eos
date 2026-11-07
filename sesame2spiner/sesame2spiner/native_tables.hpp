//------------------------------------------------------------------------------
// © 2026. Triad National Security, LLC. All rights reserved. This program was
// produced under U.S. Government contract 89233218CNA000001 for Los Alamos
// National Laboratory (LANL), operated by Triad National Security, LLC for the
// U.S. Department of Energy/National Nuclear Security Administration. The
// Government is granted a nonexclusive, paid-up, irrevocable worldwide license
// to reproduce, prepare derivative works, distribute, perform and display this
// material, and permit others to do so.
//------------------------------------------------------------------------------
#ifndef SESAME2SPINER_NATIVE_TABLES_HPP_
#define SESAME2SPINER_NATIVE_TABLES_HPP_

#include <string>
#include <utility>
#include <vector>

#include <eospac-wrapper/eospac_wrapper.hpp>
#include <ports-of-call/portability.hpp>
#include <singularity-utils/constants.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation_stencils.hpp>
#include <spiner/nonuniform_grid_1d.hpp>
#include <spiner/piecewise_grid_1d.hpp>

namespace sesame2spiner {

using NativeGrid = Spiner::NonUniformGrid1D<Real>;
using NativeTable =
    Spiner::DataBox<Real, NativeGrid, Spiner::InterpolationStencils::Rational>;
using TemperatureGuessGrid = Spiner::PiecewiseGrid1D<Real, 3>;
using TemperatureGuessTable =
    Spiner::DataBox<Real, TemperatureGuessGrid, Spiner::InterpolationStencils::Linear>;

struct TemperatureGuessOptions {
  // Zero selects twice as many intervals as the native rho-T energy table.
  // A supplied piecewise grid overrides the corresponding point count.
  int rho_points = 0;
  int sie_points = 0;
  TemperatureGuessGrid rho_grid, sie_grid;
};

struct NativeMetadata {
  int matid = 0;
  Real normalDensity = 0;
  Real meanAtomicMass = 1;
  Real meanAtomicNumber = 1;
};

struct NativeLoadOptions {
  singularity::TableSplit split = singularity::TableSplit::Total;
  bool load_entropy = true;
  EospacWrapper::Verbosity eospacWarn = EospacWrapper::Verbosity::Quiet;
  TemperatureGuessOptions temperature_guess;
};

// Host owner of the tables returned by the loader. Copies are prohibited because
// Spiner tables use shallow copies and explicit finalization. A consumer may copy
// the table handles and call release() to assume responsibility for finalizing
// their values AND grids. A moved-from/released bundle no longer owns its tables.
struct NativeMaterial {
  NativeMetadata metadata;
  NativeTable pressure, energy, entropy;
  // Required by SesameEOS. Linear in physical rho and specific internal energy.
  // -1 marks an invalid
  // sample; finite EOSPAC extrapolations are bounded to nonnegative temperature.
  TemperatureGuessTable temperature_guess;
  // Empty on success or when entropy loading was not requested.
  std::string entropy_error;
  NativeMaterial() = default;
  NativeMaterial(const NativeMaterial &) = delete;
  NativeMaterial &operator=(const NativeMaterial &) = delete;
  NativeMaterial &operator=(NativeMaterial &&) = delete;
  NativeMaterial(NativeMaterial &&other) noexcept
      : metadata(other.metadata), pressure(other.pressure), energy(other.energy),
        entropy(other.entropy), temperature_guess(other.temperature_guess),
        entropy_error(std::move(other.entropy_error)), owns_(other.owns_) {
    other.release();
  }
  ~NativeMaterial() {
    if (owns_) {
      for (auto *table : {&pressure, &energy, &entropy}) {
        table->finalizeGrids();
        if (table->data() != nullptr && table->ownsAllocatedMemory()) table->finalize();
      }
      temperature_guess.finalizeGrids();
      if (temperature_guess.data() != nullptr && temperature_guess.ownsAllocatedMemory())
        temperature_guess.finalize();
    }
  }
  void release() noexcept { owns_ = false; }

 private:
  bool owns_ = true;
};

// Values use Singularity units and temperature is the fastest index:
// values[irho * temperatures.size() + iT]. Axes remain in physical coordinates.
// The caller owns the returned values and grids and must explicitly finalize both.
// Construction/validation failures use PORTABLE_ALWAYS_THROW_OR_ABORT: serial
// builds throw; other portability backends abort, matching EOSPAC construction.
NativeTable makeNativeTable(const std::vector<Real> &densities,
                            const std::vector<Real> &temperatures,
                            const std::vector<Real> &values);

// A caller-supplied inverse approximation, in Singularity units and rho-major,
// energy-fast ordering. Returns owned values; finalize values and grids together.
// Values must be finite nonnegative temperatures or the invalid sentinel -1.
TemperatureGuessTable makeTemperatureGuessTable(const TemperatureGuessGrid &rho_grid,
                                                const TemperatureGuessGrid &sie_grid,
                                                const std::vector<Real> &temperatures);

// Read native nodal arrays, rather than sampling an EOSPAC interpolant. Filename
// selects the actual SESAME database, not a sesame2spiner parameter deck. EOSPAC
// is used only during this call; no EOSPAC handles survive in the result.
// Required failures follow the portable error policy described above. Entropy
// load failures are optional and leave an empty table with an entropy_error.
NativeMaterial loadNativeMaterial(const std::string &filename, int matid,
                                  const NativeLoadOptions &options = {});

} // namespace sesame2spiner
#endif // SESAME2SPINER_NATIVE_TABLES_HPP_
