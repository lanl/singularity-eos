//------------------------------------------------------------------------------
// © 2026. Triad National Security, LLC. All rights reserved. This program was
// produced under U.S. Government contract 89233218CNA000001 for Los Alamos
// National Laboratory (LANL), operated by Triad National Security, LLC for the
// U.S. Department of Energy/National Nuclear Security Administration. The
// Government is granted a nonexclusive, paid-up, irrevocable worldwide license
// to reproduce, prepare derivative works, distribute, perform and display this
// material, and permit others to do so.
//------------------------------------------------------------------------------
#include "native_tables.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <vector>

#include <eos_Interface.h>
#include <eospac-wrapper/eospac_wrapper.hpp>
#include <ports-of-call/portable_errors.hpp>

namespace sesame2spiner {
namespace {

// These calls and array conventions follow EOSPAC 6's user manual, chapters 5
// and 6 and appendix E: https://www.osti.gov/servlets/purl/1489917 .
class TableHandle {
 public:
  TableHandle(const std::string &filename, int matid, EOS_INTEGER type,
              EospacWrapper::Verbosity eospacWarn, bool required = true)
      : context_(filename + ": material " + std::to_string(matid) + ", table " +
                 std::to_string(type)),
        eospacWarn_(eospacWarn) {
    const EOS_INTEGER error = EospacWrapper::eosSafeLoad(
        filename, 1, matid, &type, &handle_, {context_}, eospacWarn_);
    if (!ok(error)) {
      destroy();
      if (required) fail("eosSafeLoad", error);
      diagnostic_ = context_ + ": eosSafeLoad: " + EospacWrapper::eosErrorString(error);
    }
    loaded_ = ok(error);
  }
  TableHandle(const TableHandle &) = delete;
  TableHandle &operator=(const TableHandle &) = delete;
  ~TableHandle() { destroy(); }

  bool loaded() const { return loaded_; }
  const std::string &diagnostic() const { return diagnostic_; }
  EOS_REAL scalar(EOS_INTEGER item) { return array(item, 1)[0]; }
  std::vector<EOS_REAL> array(EOS_INTEGER item, std::size_t size) {
    std::vector<EOS_REAL> result(size);
    // For array flags numInfoItems is ONE; InfoVals must hold the entire array.
    EOS_INTEGER count = 1, error = EOS_OK;
    EospacWrapper::eosSafeTableInfo(&handle_, count, &item, result.data(), error,
                                    eospacWarn_);
    if (!ok(error)) fail("eos_GetTableInfo item " + std::to_string(item), error);
    return result;
  }
  void interpolateTemperature(std::vector<EOS_REAL> &rhos,
                              std::vector<EOS_REAL> &energies,
                              std::vector<EOS_REAL> &temperatures) {
    std::vector<EOS_REAL> dx(rhos.size()), dy(rhos.size());
    EOS_INTEGER status = EOS_OK;
    EospacWrapper::eosSafeInterpolate(
        &handle_, static_cast<EOS_INTEGER>(rhos.size()), rhos.data(), energies.data(),
        temperatures.data(), dx.data(), dy.data(), context_.c_str(), eospacWarn_, status);
    // An extrapolated inverse is useful as an initial guess. Other failures
    // must not silently produce a table of uninitialized temperatures.
    if (!ok(status) &&
        !EospacWrapper::eosErrorCodesEqual(status, EOS_INTERP_EXTRAPOLATED))
      fail("eosSafeInterpolate", status);
  }

 private:
  EOS_INTEGER handle_ = -1;
  std::string context_;
  std::string diagnostic_;
  EospacWrapper::Verbosity eospacWarn_;
  bool loaded_ = false;
  static bool ok(EOS_INTEGER error) {
    return EospacWrapper::eosErrorCodesEqual(error, EOS_OK);
  }
  [[noreturn]] void fail(const std::string &operation, EOS_INTEGER error) const {
    PORTABLE_ALWAYS_THROW_OR_ABORT(context_ + ": " + operation + ": " +
                                   EospacWrapper::eosErrorString(error));
  }
  void destroy() noexcept {
    if (handle_ >= 0) {
      EospacWrapper::eosSafeDestroy(1, &handle_, eospacWarn_);
      handle_ = -1;
    }
  }
};

int dimension(EOS_REAL n) {
  if (!std::isfinite(n) || n < 4 ||
      n > static_cast<EOS_REAL>(std::numeric_limits<int>::max()) || std::floor(n) != n) {
    PORTABLE_ALWAYS_THROW_OR_ABORT(
        "Native rational tables require at least four points per axis");
  }
  return static_cast<int>(n);
}

void extract(TableHandle &handle, const Real output_scale, NativeTable &table) {
  const int nr = dimension(handle.scalar(EOS_NR));
  const int nt = dimension(handle.scalar(EOS_NT));
  if (nr > std::numeric_limits<int>::max() / nt / static_cast<int>(sizeof(Real))) {
    PORTABLE_ALWAYS_THROW_OR_ABORT(
        "Native table exceeds Spiner's supported storage size");
  }
  const auto r = handle.array(EOS_R_Array, nr);
  const auto t = handle.array(EOS_T_Array, nt);
  const auto f = handle.array(EOS_F_Array, static_cast<std::size_t>(nr) * nt);
  std::vector<Real> densities(nr), temperatures(nt);
  for (int ir = 0; ir < nr; ++ir)
    densities[ir] = EospacWrapper::densityFromSesame(r[ir]);
  for (int it = 0; it < nt; ++it)
    temperatures[it] = EospacWrapper::temperatureFromSesame(t[it]);
  std::vector<Real> values(f.size());
  // EOSPAC exposes F(NR,NT) in Fortran column-major order: density moves fastest.
  // Spiner's native rho-T representation makes temperature move fastest.
  for (int ir = 0; ir < nr; ++ir) {
    for (int it = 0; it < nt; ++it) {
      values[static_cast<std::size_t>(ir) * nt + it] =
          output_scale * f[static_cast<std::size_t>(it) * nr + ir];
    }
  }
  // Density and temperature are already in Singularity units; output conversion
  // is applied exactly once. No interpolation, inserted points, or smoothing.
  table = makeNativeTable(densities, temperatures, values);
}

EOS_INTEGER select(singularity::TableSplit split, EOS_INTEGER total, EOS_INTEGER electron,
                   EOS_INTEGER ion) {
  if (split == singularity::TableSplit::Total) return total;
  if (split == singularity::TableSplit::ElectronOnly) return electron;
  if (split == singularity::TableSplit::IonCold) return ion;
  PORTABLE_ALWAYS_THROW_OR_ABORT("Invalid native SESAME table split");
}

void checkGuessGrid(const TemperatureGuessGrid &grid, bool density) {
  if (!grid.isWellFormed() || grid.nPoints() < 2 ||
      grid.nPoints() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      !std::isfinite(grid.min()) || !std::isfinite(grid.max()) ||
      grid.max() <= grid.min() || (density && grid.min() < 0)) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Invalid temperature-guess grid");
  }
  for (int i = 0; i < static_cast<int>(grid.nPoints()); ++i) {
    // Piecewise grids store both copies of a shared endpoint. Duplicate join
    // coordinates are valid; each subgrid still needs a finite positive spacing.
    if (!std::isfinite(grid.x(i)) || !std::isfinite(grid.dx(i)) || grid.dx(i) <= 0 ||
        !std::isfinite(grid.idx(i)) || grid.idx(i) <= 0 ||
        (i > 0 && grid.x(i) < grid.x(i - 1)))
      PORTABLE_ALWAYS_THROW_OR_ABORT("Invalid temperature-guess grid coordinates");
  }
}

void allocateGuessTable(TemperatureGuessTable &table,
                        const TemperatureGuessGrid &rho_grid,
                        const TemperatureGuessGrid &sie_grid) {
  checkGuessGrid(rho_grid, true);
  checkGuessGrid(sie_grid, false);
  const auto nr = rho_grid.nPoints(), ne = sie_grid.nPoints();
  if (nr > std::numeric_limits<int>::max() / sizeof(Real) / ne)
    PORTABLE_ALWAYS_THROW_OR_ABORT(
        "Temperature-guess table exceeds supported storage size");
  table = TemperatureGuessTable(static_cast<int>(nr), static_cast<int>(ne));
  table.setRange(0, sie_grid);
  table.setRange(1, rho_grid);
}

TemperatureGuessGrid defaultGuessGrid(Real minimum, Real maximum, int requested,
                                      int native_points) {
  if (!std::isfinite(minimum) || !std::isfinite(maximum) || minimum >= maximum ||
      requested < 0 || native_points > std::numeric_limits<int>::max() / 2)
    PORTABLE_ALWAYS_THROW_OR_ABORT("Invalid temperature-guess bounds or point count");
  const int points = requested == 0 ? 2 * (native_points - 1) + 1 : requested;
  if (points < 2)
    PORTABLE_ALWAYS_THROW_OR_ABORT("Temperature-guess grids need at least two points");
  return TemperatureGuessGrid{{Spiner::RegularGrid1D<Real>(minimum, maximum, points)}};
}

void populateTemperatureGuess(const std::string &filename, int matid,
                              const NativeLoadOptions &options,
                              NativeMaterial &material) {
  const auto &configuration = options.temperature_guess;
  const auto &energy = material.energy;
  const auto rho_grid =
      configuration.rho_grid.nGrids() > 0
          ? configuration.rho_grid
          : defaultGuessGrid(energy.range(1).min(), energy.range(1).max(),
                             configuration.rho_points, energy.dim(2));
  const auto sie_grid = configuration.sie_grid.nGrids() > 0
                            ? configuration.sie_grid
                            : defaultGuessGrid(energy.min(), energy.max(),
                                               configuration.sie_points, energy.dim(1));
  allocateGuessTable(material.temperature_guess, rho_grid, sie_grid);
  TableHandle inverse(filename, matid,
                      select(options.split, EOS_T_DUt, EOS_T_DUe, EOS_T_DUic),
                      options.eospacWarn);

  // Keep scratch bounded even for a finer table. As in eosDataOfRhoSie, batch
  // (rho,e) pairs in SESAME units and convert returned temperatures exactly once.
  constexpr int batch_size = 4096;
  const int ne = static_cast<int>(sie_grid.nPoints());
  auto &table = material.temperature_guess;
  int usable = 0;
  for (int start = 0; start < table.size();) {
    const int count = std::min(batch_size, table.size() - start);
    std::vector<EOS_REAL> rhos(count), energies(count),
        temperatures(count, std::numeric_limits<EOS_REAL>::quiet_NaN());
    for (int k = 0; k < count; ++k) {
      rhos[k] = EospacWrapper::densityToSesame(rho_grid.x((start + k) / ne));
      energies[k] = EospacWrapper::sieToSesame(sie_grid.x((start + k) % ne));
    }
    inverse.interpolateTemperature(rhos, energies, temperatures);
    for (int k = 0; k < count; ++k) {
      const Real temp = EospacWrapper::temperatureFromSesame(temperatures[k]);
      table(start + k) = std::isfinite(temp) ? std::max(temp, Real(0)) : Real(-1);
      usable += std::isfinite(temp) ? 1 : 0;
    }
    start += count;
  }
  if (usable == 0)
    PORTABLE_ALWAYS_THROW_OR_ABORT("EOSPAC produced no usable temperature guesses");
}

} // namespace

NativeTable makeNativeTable(const std::vector<Real> &densities,
                            const std::vector<Real> &temperatures,
                            const std::vector<Real> &values) {
  const std::size_t max_values = std::numeric_limits<int>::max() / sizeof(Real);
  for (const auto *axis : {&densities, &temperatures}) {
    if (axis->size() < 4 || axis->size() > max_values) {
      PORTABLE_ALWAYS_THROW_OR_ABORT(
          "Native rational tables require at least four points per axis");
    }
    for (std::size_t i = 0; i < axis->size(); ++i) {
      if (!std::isfinite((*axis)[i]) || (*axis)[i] < 0 ||
          (i > 0 && (*axis)[i] <= (*axis)[i - 1])) {
        PORTABLE_ALWAYS_THROW_OR_ABORT(
            "Native SESAME axes must be finite, nonnegative and increasing");
      }
    }
  }
  if (densities.size() > max_values / temperatures.size() ||
      values.size() != densities.size() * temperatures.size()) {
    PORTABLE_ALWAYS_THROW_OR_ABORT("Invalid native SESAME table dimensions");
  }
  for (Real value : values) {
    if (!std::isfinite(value))
      PORTABLE_ALWAYS_THROW_OR_ABORT("Nonfinite native SESAME table value");
  }
  // Use the bundle's cleanup on allocation/validation failure; transfer the table
  // only after both grids have been built successfully.
  NativeMaterial owner;
  owner.pressure = NativeTable(static_cast<int>(densities.size()),
                               static_cast<int>(temperatures.size()));
  owner.pressure.setRange(0, NativeGrid(temperatures));
  owner.pressure.setRange(1, NativeGrid(densities));
  for (std::size_t i = 0; i < values.size(); ++i)
    owner.pressure(i) = values[i];
  NativeTable result = owner.pressure;
  owner.release();
  return result;
}

TemperatureGuessTable makeTemperatureGuessTable(const TemperatureGuessGrid &rho_grid,
                                                const TemperatureGuessGrid &sie_grid,
                                                const std::vector<Real> &temperatures) {
  NativeMaterial owner;
  allocateGuessTable(owner.temperature_guess, rho_grid, sie_grid);
  if (temperatures.size() != static_cast<std::size_t>(owner.temperature_guess.size()))
    PORTABLE_ALWAYS_THROW_OR_ABORT("Invalid temperature-guess table dimensions");
  for (std::size_t i = 0; i < temperatures.size(); ++i) {
    const Real value = temperatures[i];
    if (!std::isfinite(value) || (value < 0 && value != -1))
      PORTABLE_ALWAYS_THROW_OR_ABORT("Invalid temperature-guess sample");
    owner.temperature_guess(i) = value;
  }
  TemperatureGuessTable result = owner.temperature_guess;
  owner.release();
  return result;
}

NativeMaterial loadNativeMaterial(const std::string &filename, int matid,
                                  const NativeLoadOptions &options) {
  try {
    if (filename.empty())
      PORTABLE_ALWAYS_THROW_OR_ABORT("An explicit SESAME filename is required");
    NativeMaterial result;
    result.metadata.matid = matid;
    TableHandle p(filename, matid,
                  select(options.split, EOS_Pt_DT, EOS_Pe_DT, EOS_Pic_DT),
                  options.eospacWarn);
    TableHandle e(filename, matid,
                  select(options.split, EOS_Ut_DT, EOS_Ue_DT, EOS_Uic_DT),
                  options.eospacWarn);
    extract(p, EospacWrapper::pressureFromSesame(1), result.pressure);
    extract(e, EospacWrapper::sieFromSesame(1), result.energy);

    // Use the material's information table, selected from the SAME source file.
    TableHandle info(filename, matid, EOS_Info, options.eospacWarn);
    result.metadata.normalDensity =
        EospacWrapper::densityFromSesame(info.scalar(EOS_Normal_Density));
    result.metadata.meanAtomicMass = info.scalar(EOS_Mean_Atomic_Mass);
    result.metadata.meanAtomicNumber = info.scalar(EOS_Mean_Atomic_Num);
    if (!std::isfinite(result.metadata.normalDensity) ||
        result.metadata.normalDensity <= 0 ||
        !std::isfinite(result.metadata.meanAtomicMass) ||
        result.metadata.meanAtomicMass <= 0 ||
        !std::isfinite(result.metadata.meanAtomicNumber) ||
        result.metadata.meanAtomicNumber <= 0) {
      PORTABLE_ALWAYS_THROW_OR_ABORT("Invalid native SESAME material metadata");
    }
    if (options.load_entropy) {
      TableHandle s(filename, matid,
                    select(options.split, EOS_St_DT, EOS_Se_DT, EOS_Sic_DT),
                    options.eospacWarn, false);
      // Like EOSPAC's EOS class, an unavailable entropy table does not prevent
      // using pressure and energy. Successfully loaded entropy must be valid.
      if (s.loaded())
        extract(s, EospacWrapper::entropyFromSesame(1), result.entropy);
      else
        result.entropy_error = s.diagnostic();
    }
    populateTemperatureGuess(filename, matid, options, result);
    return result;
  } catch (const std::exception &error) {
    PORTABLE_ALWAYS_THROW_OR_ABORT(filename + ": material " + std::to_string(matid) +
                                   ": " + error.what());
  }
}

} // namespace sesame2spiner
