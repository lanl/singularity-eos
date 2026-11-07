//------------------------------------------------------------------------------
// © 2026. Triad National Security, LLC. All rights reserved. This program was
// produced under U.S. Government contract 89233218CNA000001 for Los Alamos
// National Laboratory (LANL), operated by Triad National Security, LLC for the
// U.S. Department of Energy/National Nuclear Security Administration. The
// Government is granted a nonexclusive, paid-up, irrevocable worldwide license
// to reproduce, prepare derivative works, distribute, perform and display this
// material, and permit others to do so.
//------------------------------------------------------------------------------
#ifndef SINGULARITY_EOS_EOS_SESAME_HPP_
#define SINGULARITY_EOS_EOS_SESAME_HPP_
#ifdef SINGULARITY_USE_SESAME2SPINER

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <sesame2spiner/native_tables.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/eos/eos_base.hpp>
#include <singularity-utils/indexable_types.hpp>
#include <singularity-utils/robust_utils.hpp>
#include <singularity-utils/root-finding-1d/root_finding.hpp>
#include <singularity-utils/variadic_utils.hpp>

namespace singularity {
using namespace eos_base;

// SESAME's native physical rho-T grids and Spiner's Kerley rational interpolant.
// EOSPAC is called only by the host loader. Queries use native interpolants;
// converged inversions validate the forward energy residual. An auxiliary
// temperature table supplies a Newton seed and the failed-solve fallback.
// Derivatives come from the native interpolant; optional RootStatus reports
// whether an inverse query converged or returned its guess, as in SpinerEOS.
// Copies follow the usual EOS shallow-copy / explicit-Finalize convention.
class SesameEOS : public EosBase<SesameEOS> {
 public:
  using DataBox = sesame2spiner::NativeTable;
  using TemperatureGuessTable = sesame2spiner::TemperatureGuessTable;
  using Options = sesame2spiner::NativeLoadOptions;
  using InterpResult = DataBox::InterpResult2D;
  SG_ADD_BASE_CLASS_USINGS(SesameEOS)
  SG_ADD_DEFAULT_MEAN_ATOMIC_FUNCTIONS(AZbar_)
  using EosBase<SesameEOS>::DensityEnergyFromPressureTemperature;

  // The generic EosBase pointer/scratch wrappers ignore Transform. Implement
  // this overload explicitly so modifiers can scale arguments and outputs.
#define SG_SESAME_TRANSFORMED(NAME)                                                      \
  template <typename Space, typename LambdaIndexer,                                      \
            typename EnableIfSpace =                                                     \
                std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>               \
  inline void NAME(const Space &space, const Real *x, const Real *y, Real *out,          \
                   Real * /*scratch*/, const int num, LambdaIndexer &&lambdas,           \
                   Transform &&transform = Transform()) const {                          \
    const Real sx = transform.x.is_set() ? transform.x.get() : Real(1);                  \
    const Real sy = transform.y.is_set() ? transform.y.get() : Real(1);                  \
    const Real sf = transform.f.is_set() ? transform.f.get() : Real(1);                  \
    const SesameEOS eos = *this;                                                         \
    portableFor(                                                                         \
        "SesameEOS::" #NAME, space, 0, num, PORTABLE_LAMBDA(const int i) {               \
          out[i] = sf * eos.NAME(sx * x[i], sy * y[i], lambdas[i]);                      \
        });                                                                              \
  }                                                                                      \
  template <typename LambdaIndexer>                                                      \
  inline void NAME(const Real *x, const Real *y, Real *out, Real *scratch,               \
                   const int num, LambdaIndexer &&lambdas,                               \
                   Transform &&transform = Transform()) const {                          \
    NAME(PortsOfCall::Exec::Device(), x, y, out, scratch, num,                           \
         std::forward<LambdaIndexer>(lambdas), std::forward<Transform>(transform));      \
  }
  SG_SESAME_TRANSFORMED(TemperatureFromDensityInternalEnergy)
  SG_SESAME_TRANSFORMED(InternalEnergyFromDensityTemperature)
  SG_SESAME_TRANSFORMED(PressureFromDensityTemperature)
  SG_SESAME_TRANSFORMED(PressureFromDensityInternalEnergy)
  SG_SESAME_TRANSFORMED(EntropyFromDensityTemperature)
  SG_SESAME_TRANSFORMED(EntropyFromDensityInternalEnergy)
  SG_SESAME_TRANSFORMED(SpecificHeatFromDensityTemperature)
  SG_SESAME_TRANSFORMED(SpecificHeatFromDensityInternalEnergy)
  SG_SESAME_TRANSFORMED(BulkModulusFromDensityTemperature)
  SG_SESAME_TRANSFORMED(BulkModulusFromDensityInternalEnergy)
  SG_SESAME_TRANSFORMED(GruneisenParamFromDensityTemperature)
  SG_SESAME_TRANSFORMED(GruneisenParamFromDensityInternalEnergy)
  SG_SESAME_TRANSFORMED(GibbsFreeEnergyFromDensityTemperature)
  SG_SESAME_TRANSFORMED(GibbsFreeEnergyFromDensityInternalEnergy)
#undef SG_SESAME_TRANSFORMED

  template <
      typename Space, typename LambdaIndexer,
      typename EnableIfSpace = std::enable_if_t<!variadic_utils::has_int_index_v<Space>>>
  inline void MinInternalEnergyFromDensity(const Space &space, const Real *rho, Real *sie,
                                           Real * /*scratch*/, const int num,
                                           LambdaIndexer &&lambdas,
                                           Transform &&transform = Transform()) const {
    const Real sx = transform.x.is_set() ? transform.x.get() : Real(1);
    const Real sf = transform.f.is_set() ? transform.f.get() : Real(1);
    const SesameEOS eos = *this;
    portableFor(
        "SesameEOS::MinInternalEnergyFromDensity", space, 0, num,
        PORTABLE_LAMBDA(const int i) {
          sie[i] = sf * eos.MinInternalEnergyFromDensity(sx * rho[i], lambdas[i]);
        });
  }
  template <typename LambdaIndexer>
  inline void MinInternalEnergyFromDensity(const Real *rho, Real *sie, Real *scratch,
                                           const int num, LambdaIndexer &&lambdas,
                                           Transform &&transform = Transform()) const {
    MinInternalEnergyFromDensity(PortsOfCall::Exec::Device(), rho, sie, scratch, num,
                                 std::forward<LambdaIndexer>(lambdas),
                                 std::forward<Transform>(transform));
  }

  PORTABLE_INLINE_FUNCTION SesameEOS() = default;
  inline SesameEOS(const std::string &filename, int matid, const Options &options = {})
      : SesameEOS(sesame2spiner::loadNativeMaterial(filename, matid, options)) {}
  inline SesameEOS(const std::string &filename, int matid, TableSplit split)
      : SesameEOS(filename, matid, Options{split, true}) {}

  // Also permits construction from validated in-memory native tables. Ownership
  // is transferred from the host bundle; no copying of table data is required.
  explicit inline SesameEOS(sesame2spiner::NativeMaterial &&material) {
    for (const auto *table : {&material.pressure, &material.energy}) {
      PORTABLE_ALWAYS_REQUIRE(
          table->rank() == 2 && table->dim(1) >= 4 && table->dim(2) >= 4,
          "SesameEOS requires native two-dimensional rational tables");
      PORTABLE_ALWAYS_REQUIRE(table->range(0).isWellFormed() &&
                                  table->range(1).isWellFormed(),
                              "SesameEOS requires valid native axes");
    }
    if (material.entropy.data() != nullptr) {
      PORTABLE_ALWAYS_REQUIRE(
          material.entropy.rank() == 2 && material.entropy.dim(1) >= 4 &&
              material.entropy.dim(2) >= 4 && material.entropy.range(0).isWellFormed() &&
              material.entropy.range(1).isWellFormed(),
          "SesameEOS requires a valid entropy table");
    }
    const auto &guess = material.temperature_guess;
    PORTABLE_ALWAYS_REQUIRE(
        guess.data() != nullptr && guess.rank() == 2 && guess.dim(1) >= 2 &&
            guess.dim(2) >= 2 && guess.range(0).isWellFormed() &&
            guess.range(1).isWellFormed(),
        "SesameEOS requires an auxiliary temperature table with valid grids");
    const auto &p = material.pressure;
    const auto &e = material.energy;
    rho_min_ = std::max(p.range(1).min(), e.range(1).min());
    rho_max_ = std::min(p.range(1).max(), e.range(1).max());
    temp_min_ = std::max(p.range(0).min(), e.range(0).min());
    temp_max_ = std::min(p.range(0).max(), e.range(0).max());
    sie_min_ = e.min();
    sie_max_ = e.max();
    press_min_ = p.min();
    rho_ref_ = material.metadata.normalDensity;
    matid_ = material.metadata.matid;
    AZbar_ = MeanAtomicProperties(material.metadata.meanAtomicMass,
                                  material.metadata.meanAtomicNumber);
    CheckBounds_();
    AZbar_.CheckParams();
    PORTABLE_ALWAYS_REQUIRE(std::isfinite(rho_ref_) && rho_ref_ > 0,
                            "SesameEOS requires a positive reference density");
    tables_[Pressure] = material.pressure;
    tables_[Energy] = material.energy;
    tables_[Entropy] = material.entropy;
    temperature_guess_ = material.temperature_guess;
    entropy_available_ = material.entropy.data() != nullptr;
    owns_ = true;
    initialized_ = true;
    material.release();
  }

  PORTABLE_INLINE_FUNCTION void CheckParams() const {
    PORTABLE_ALWAYS_REQUIRE(initialized_, "SesameEOS is not initialized");
    PORTABLE_ALWAYS_REQUIRE(TemperatureGuessAvailable(),
                            "SesameEOS requires its temperature table");
    CheckBounds_();
    AZbar_.CheckParams();
  }
  static std::string EosType() { return "SesameEOS"; }
  static std::string EosPyType() { return EosType(); }
  static constexpr unsigned long PreferredInput() {
    return thermalqs::density | thermalqs::temperature;
  }
  PORTABLE_INLINE_FUNCTION int matid() const { return matid_; }
  PORTABLE_INLINE_FUNCTION bool EntropyAvailable() const { return entropy_available_; }
  PORTABLE_INLINE_FUNCTION bool TemperatureGuessAvailable() const {
    return temperature_guess_.data() != nullptr;
  }
  PORTABLE_INLINE_FUNCTION Real MinimumDensity() const { return rho_min_; }
  PORTABLE_INLINE_FUNCTION Real MaximumDensity() const { return rho_max_; }
  PORTABLE_INLINE_FUNCTION Real MinimumTemperature() const { return temp_min_; }
  PORTABLE_INLINE_FUNCTION Real MaximumTemperature() const { return temp_max_; }
  PORTABLE_INLINE_FUNCTION Real MinimumPressure() const { return press_min_; }
  PORTABLE_INLINE_FUNCTION Real MinimumInternalEnergy() const { return sie_min_; }
  PORTABLE_INLINE_FUNCTION Real MaximumInternalEnergy() const { return sie_max_; }
  PORTABLE_INLINE_FUNCTION void PrintParams() const {
    printf("SesameEOS material %d: rho=[%g,%g], T=[%g,%g]\n", matid_, rho_min_, rho_max_,
           temp_min_, temp_max_);
    AZbar_.PrintParams();
  }

  // Expose value AND physical-coordinate gradients without a finite difference.
  PORTABLE_INLINE_FUNCTION InterpResult PressureWithDerivatives(Real rho,
                                                                Real temp) const {
    return Evaluate_(Pressure, rho, temp);
  }
  PORTABLE_INLINE_FUNCTION InterpResult EnergyWithDerivatives(Real rho, Real temp) const {
    return Evaluate_(Energy, rho, temp);
  }

  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityTemperature(
      Real rho, Real temp, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return Evaluate_(Energy, rho, temp).value;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityTemperature(
      Real rho, Real temp, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return Evaluate_(Pressure, rho, temp).value;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real TemperatureFromDensityInternalEnergy(
      Real rho, Real sie, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    PORTABLE_ALWAYS_REQUIRE(initialized_ && TemperatureGuessAvailable(),
                            "SesameEOS requires initialized temperature tables");
    PORTABLE_ALWAYS_REQUIRE(std::isfinite(rho) && rho >= 0 && std::isfinite(sie),
                            "SesameEOS: invalid temperature inversion inputs");
    // Like SpinerEOS, return the initial guess if the numerical solve fails.
    // If the seed cell itself is unusable, use the minimum table temperature.
    Real guess = temp_min_, temperature = guess;
    const bool converged =
        TryTemperatureFromDensityInternalEnergy_(rho, sie, temperature, guess);
    IndexerUtils::SafeSet<IndexableTypes::RootStatus>(
        lambda, static_cast<Real>(converged ? RootFinding1D::Status::SUCCESS
                                            : RootFinding1D::Status::FAIL));
    return converged ? temperature : guess;
  }

  // Failure leaves temperature unchanged. Use the shared Newton routine without
  // a cell scan, bracket search, or fallback to another inversion method.
  PORTABLE_INLINE_FUNCTION bool
  TryTemperatureFromDensityInternalEnergy(Real rho, Real sie, Real &temperature) const {
    Real guess = temp_min_;
    return TryTemperatureFromDensityInternalEnergy_(rho, sie, temperature, guess);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real PressureFromDensityInternalEnergy(
      Real rho, Real sie, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return PressureFromDensityTemperature(
        rho, TemperatureFromDensityInternalEnergy(rho, sie, lambda), lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real MinInternalEnergyFromDensity(
      Real rho, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    // Thermal-branch bound on this native table, without a separate cold table.
    return InternalEnergyFromDensityTemperature(rho, temp_min_, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityTemperature(
      Real rho, Real temp, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    if (!entropy_available_) EntropyIsNotEnabled("SesameEOS");
    return Evaluate_(Entropy, rho, temp).value;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real EntropyFromDensityInternalEnergy(
      Real rho, Real sie, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return EntropyFromDensityTemperature(
        rho, TemperatureFromDensityInternalEnergy(rho, sie, lambda), lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityTemperature(
      Real rho, Real temp, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return std::max(Evaluate_(Energy, rho, temp).d_dx1, Real(0));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real SpecificHeatFromDensityInternalEnergy(
      Real rho, Real sie, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return SpecificHeatFromDensityTemperature(
        rho, TemperatureFromDensityInternalEnergy(rho, sie, lambda), lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityTemperature(
      Real rho, Real temp, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return BulkModulus_(rho, Evaluate_(Pressure, rho, temp),
                        Evaluate_(Energy, rho, temp));
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real BulkModulusFromDensityInternalEnergy(
      Real rho, Real sie, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return BulkModulusFromDensityTemperature(
        rho, TemperatureFromDensityInternalEnergy(rho, sie, lambda), lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityTemperature(
      Real rho, Real temp, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    const auto e = Evaluate_(Energy, rho, temp);
    const auto p = Evaluate_(Pressure, rho, temp);
    return rho > 0 ? robust::ratio(robust::ratio(p.d_dx1, e.d_dx1), rho) : Real(0);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real GruneisenParamFromDensityInternalEnergy(
      Real rho, Real sie, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    return GruneisenParamFromDensityTemperature(
        rho, TemperatureFromDensityInternalEnergy(rho, sie, lambda), lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION Real InternalEnergyFromDensityPressure(
      Real rho, Real press, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    PORTABLE_ALWAYS_REQUIRE(initialized_ && std::isfinite(press),
                            "SesameEOS: invalid pressure inversion inputs");
    // Preserve an exact thermal-boundary state before the outer pressure solve:
    // the inner energy inversion intentionally accepts a finite residual tolerance.
    if (press == PressureFromDensityTemperature(rho, temp_min_, lambda))
      return InternalEnergyFromDensityTemperature(rho, temp_min_, lambda);
    // Shared pressure inversion uses the table-backed temperature solve through
    // PressureFromDensityInternalEnergy. Start from a known reference state.
    Real sie = InternalEnergyFromDensityTemperature(rho, ROOM_TEMPERATURE, lambda);
    EosBase<SesameEOS>::InternalEnergyFromDensityPressure(rho, press, sie, lambda);
    return sie;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void InternalEnergyFromDensityPressure(
      Real rho, Real press, Real &sie,
      Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    sie = InternalEnergyFromDensityPressure(rho, press, lambda);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  DensityEnergyFromPressureTemperature(Real press, Real temp, Indexer_t &&lambda,
                                       Real &rho, Real &sie) const {
    PORTABLE_ALWAYS_REQUIRE(initialized_ && std::isfinite(press) && std::isfinite(temp) &&
                                temp >= 0,
                            "SesameEOS: invalid pressure inversion inputs");
    // Preserve the output-only density convention by supplying a known seed.
    Real density = rho_ref_;
    EosBase<SesameEOS>::DensityEnergyFromPressureTemperature(press, temp, lambda, density,
                                                             sie);
    rho = density;
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  PTDerivativesFromPreferred(Real rho, Real sie, Real press, Real temp,
                             Indexer_t &&lambda, Real &dedP_T, Real &drdP_T, Real &dedT_P,
                             Real &drdT_P) const {
    const auto p = Evaluate_(Pressure, rho, temp);
    const auto e = Evaluate_(Energy, rho, temp);
    dedP_T = robust::ratio(e.d_dx2, p.d_dx2);
    drdP_T = robust::ratio(Real(1), p.d_dx2);
    dedT_P = e.d_dx1 - robust::ratio(e.d_dx2 * p.d_dx1, p.d_dx2);
    drdT_P = -robust::ratio(p.d_dx1, p.d_dx2);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  FillEos(Real &rho, Real &temp, Real &sie, Real &press, Real &cv, Real &bmod,
          unsigned long output, Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    constexpr unsigned long supported = thermalqs::density | thermalqs::temperature |
                                        thermalqs::specific_internal_energy |
                                        thermalqs::pressure | thermalqs::specific_heat |
                                        thermalqs::bulk_modulus;
    PORTABLE_ALWAYS_REQUIRE(output != 0 && (output & ~supported) == 0,
                            "SesameEOS: unsupported FillEos output mask");
    const auto input = ~output;
    if (output & thermalqs::density) {
      PORTABLE_ALWAYS_REQUIRE((input & thermalqs::temperature) &&
                                  (input & thermalqs::pressure),
                              "SesameEOS: density requires pressure and temperature");
      Real unused_energy = 0;
      DensityEnergyFromPressureTemperature(press, temp, lambda, rho, unused_energy);
    }
    if (output & thermalqs::temperature) {
      PORTABLE_ALWAYS_REQUIRE(input & thermalqs::density,
                              "SesameEOS: temperature requires density");
      if (input & thermalqs::specific_internal_energy) {
        temp = TemperatureFromDensityInternalEnergy(rho, sie, lambda);
      } else {
        PORTABLE_ALWAYS_REQUIRE(input & thermalqs::pressure,
                                "SesameEOS: temperature requires energy or pressure");
        const Real solved_energy = InternalEnergyFromDensityPressure(rho, press, lambda);
        temp = TemperatureFromDensityInternalEnergy(rho, solved_energy, lambda);
      }
    }
    InterpResult p{}, e{};
    if (output & (thermalqs::specific_internal_energy | thermalqs::specific_heat |
                  thermalqs::bulk_modulus))
      e = Evaluate_(Energy, rho, temp);
    if (output & (thermalqs::pressure | thermalqs::bulk_modulus))
      p = Evaluate_(Pressure, rho, temp);
    if (output & thermalqs::specific_internal_energy) sie = e.value;
    if (output & thermalqs::pressure) press = p.value;
    if (output & thermalqs::specific_heat) cv = std::max(e.d_dx1, Real(0));
    if (output & thermalqs::bulk_modulus) bmod = BulkModulus_(rho, p, e);
  }
  template <typename Indexer_t = Real *>
  PORTABLE_INLINE_FUNCTION void
  ValuesAtReferenceState(Real &rho, Real &temp, Real &sie, Real &press, Real &cv,
                         Real &bmod, Real &dpde, Real &dvdt,
                         Indexer_t &&lambda = static_cast<Real *>(nullptr)) const {
    rho = rho_ref_;
    temp = ROOM_TEMPERATURE;
    const auto p = Evaluate_(Pressure, rho, temp);
    const auto e = Evaluate_(Energy, rho, temp);
    sie = e.value;
    press = p.value;
    cv = std::max(e.d_dx1, Real(0));
    bmod = BulkModulus_(rho, p, e);
    dpde = robust::ratio(p.d_dx1, e.d_dx1);
    dvdt = robust::ratio(p.d_dx1, rho * rho * p.d_dx2);
  }

  PORTABLE_INLINE_FUNCTION Real MaximumPressureAtTemperature(Real temp) const {
    const auto &p = tables_[Pressure];
    Real maximum = -std::numeric_limits<Real>::max();
    for (int i = 0; i < p.dim(2); ++i) {
      maximum = std::max(maximum, Evaluate_(Pressure, p.range(1).x(i), temp).value);
    }
    return maximum;
  }
  PORTABLE_INLINE_FUNCTION Real RhoPmin(Real temp) const {
    const auto &p = tables_[Pressure];
    Real rho = rho_min_, minimum = std::numeric_limits<Real>::max();
    for (int i = 0; i < p.dim(2); ++i) {
      const Real r = p.range(1).x(i);
      const Real pressure = Evaluate_(Pressure, r, temp).value;
      if (pressure < minimum) {
        minimum = pressure;
        rho = r;
      }
    }
    return rho;
  }

  inline SesameEOS GetOnDevice() {
    CheckParams();
    SesameEOS result = *this;
    for (int i = 0; i < NumTables; ++i)
      result.tables_[i] = tables_[i].getOnDevice();
    result.temperature_guess_ = temperature_guess_.getOnDevice();
    result.owns_ = true;
    return result;
  }
  inline void Finalize() {
    for (auto &table : tables_) {
      if (owns_) {
        table.finalizeGrids();
        if (table.data() != nullptr && table.ownsAllocatedMemory()) table.finalize();
      }
      table = DataBox{};
    }
    if (owns_) {
      temperature_guess_.finalizeGrids();
      if (temperature_guess_.data() != nullptr &&
          temperature_guess_.ownsAllocatedMemory())
        temperature_guess_.finalize();
    }
    temperature_guess_ = TemperatureGuessTable{};
    owns_ = false;
    initialized_ = false;
  }
  inline std::size_t DynamicMemorySizeInBytes() const {
    std::size_t size = 0;
    for (const auto &table : tables_)
      size += table.dynamicMemorySizeInBytes();
    return size + temperature_guess_.dynamicMemorySizeInBytes();
  }
  inline std::size_t DumpDynamicMemory(char *dst) {
    std::size_t offset = 0;
    for (const auto &table : tables_)
      offset += table.dumpDynamicMemory(dst + offset);
    return offset + temperature_guess_.dumpDynamicMemory(dst + offset);
  }
  inline std::size_t
  SetDynamicMemory(char *src, const SharedMemSettings &settings = DEFAULT_SHMEM_STNGS) {
    char *data = settings.data != nullptr ? settings.data : src;
    std::size_t offset = 0;
    for (auto &table : tables_)
      offset += table.setPointer(data + offset);
    offset += temperature_guess_.setPointer(data + offset);
    owns_ = false;
    return offset;
  }

 private:
  enum Field { Pressure = 0, Energy = 1, Entropy = 2, NumTables = 3 };
  DataBox tables_[NumTables];
  TemperatureGuessTable temperature_guess_;
  MeanAtomicProperties AZbar_;
  Real rho_min_ = 0, rho_max_ = 0, temp_min_ = 0, temp_max_ = 0;
  Real sie_min_ = 0, sie_max_ = 0, press_min_ = 0, rho_ref_ = 0;
  int matid_ = 0;
  bool entropy_available_ = false, owns_ = false, initialized_ = false;

  PORTABLE_INLINE_FUNCTION bool
  TryTemperatureFromDensityInternalEnergy_(Real rho, Real sie, Real &temperature,
                                           Real &initial_guess) const {
    if (!initialized_ || !TemperatureGuessAvailable() || !std::isfinite(rho) || rho < 0 ||
        !std::isfinite(sie))
      return false;
    const auto rg = temperature_guess_.range(1), eg = temperature_guess_.range(0);
    // Interpolate inside the grid, then linearly extend its boundary value. This
    // avoids converting a huge extrapolated coordinate into a grid index.
    const Real lookup_rho = std::min(std::max(rho, rg.min()), rg.max());
    const Real lookup_sie = std::min(std::max(sie, eg.min()), eg.max());
    const int ir = rg.index(lookup_rho), ie = eg.index(lookup_sie);
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        const Real value = temperature_guess_(ir + j, ie + i);
        if (!std::isfinite(value) || value < 0) return false;
      }
    }
    const auto guess = temperature_guess_.interpToRealWithGrads(lookup_rho, lookup_sie);
    Real temp =
        guess.value + guess.d_dx2 * (rho - lookup_rho) + guess.d_dx1 * (sie - lookup_sie);
    if (!std::isfinite(temp)) return false;
    temp = std::max(temp, Real(0));
    initial_guess = temp;
    const Real scale = std::max(Real(1), std::abs(sie));
    const Real relative_tolerance =
        std::max(Real(1e-12), Real(32) * std::numeric_limits<Real>::epsilon());
    const Real tolerance = relative_tolerance * scale;
    const Real invalid = std::numeric_limits<Real>::quiet_NaN();
    // The shared solver uses |ytol*ytarget|. A unit target and normalized energy
    // residual provide an absolute-plus-relative tolerance even for zero energy.
    auto energy = [&](Real candidate) {
      if (!std::isfinite(candidate) || candidate < 0)
        return std::make_tuple(invalid, invalid);
      const auto e = Evaluate_(Energy, rho, candidate);
      if (!std::isfinite(e.value)) return std::make_tuple(invalid, invalid);
      // Accept a converged state even if the local thermal branch is flat.
      if (std::abs(e.value - sie) <= tolerance) return std::make_tuple(Real(1), Real(0));
      const Real derivative = e.d_dx1 / scale;
      if (!std::isfinite(derivative) || derivative == 0)
        return std::make_tuple(invalid, invalid);
      return std::make_tuple(Real(1) + (e.value / scale - sie / scale), derivative);
    };
    Real root = temp;
    // Permit the same nonnegative-temperature extrapolation as forward queries;
    // these are representable bounds, not a sign-changing bracket or table bounds.
    const auto status = RootFinding1D::newton_raphson(energy, Real(1), temp, Real(0),
                                                      std::numeric_limits<Real>::max(),
                                                      relative_tolerance, root);
    if (status != RootFinding1D::Status::SUCCESS || !std::isfinite(root) || root < 0)
      return false;
    const Real residual = Evaluate_(Energy, rho, root).value - sie;
    if (!std::isfinite(residual) || std::abs(residual) > tolerance) return false;
    temperature = root;
    return true;
  }

  PORTABLE_INLINE_FUNCTION void CheckBounds_() const {
    PORTABLE_ALWAYS_REQUIRE(rho_min_ >= 0 && rho_max_ > rho_min_ && temp_min_ >= 0 &&
                                temp_max_ > temp_min_ && std::isfinite(sie_min_) &&
                                std::isfinite(sie_max_) && sie_max_ >= sie_min_,
                            "SesameEOS requires ordered finite native bounds");
  }
  PORTABLE_INLINE_FUNCTION InterpResult Evaluate_(Field field, Real rho,
                                                  Real temp) const {
    PORTABLE_ALWAYS_REQUIRE(initialized_, "SesameEOS is not initialized");
    PORTABLE_ALWAYS_REQUIRE(
        std::isfinite(rho) && rho >= 0 && std::isfinite(temp) && temp >= 0,
        "SesameEOS requires finite nonnegative density and temperature");
    return tables_[field].interpToRealWithGrads(rho, temp);
  }
  PORTABLE_INLINE_FUNCTION static Real BulkModulus_(Real rho, const InterpResult &p,
                                                    const InterpResult &e) {
    if (rho <= 0) return 0;
    const Real b = e.d_dx1 > 0 ? rho * p.d_dx2 +
                                     (p.d_dx1 / e.d_dx1) * (p.value / rho - rho * e.d_dx2)
                               : rho * p.d_dx2;
    return std::max(b, Real(0));
  }
};

} // namespace singularity
#endif // SINGULARITY_USE_SESAME2SPINER
#endif // SINGULARITY_EOS_EOS_SESAME_HPP_
