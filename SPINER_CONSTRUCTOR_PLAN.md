# Plan: SpinerEOS Constructor from Generic EOS

## Goal

Add a new constructor to `SpinerEOSDependsRhoSie` that can create Spiner tables from any EOS object (analytic or otherwise) that provides the necessary thermodynamic methods. This allows users to tabulate any equation of state into the high-performance Spiner format.

## Overview

The constructor will:
1. Take a source EOS object (template parameter for flexibility)
2. Take grid parameters (mimicking sesame2spiner options)
3. Evaluate the source EOS on a grid to populate all required Spiner DataBox tables
4. Handle offsets automatically for negative values (especially sie)
5. Compute derivatives (either from source EOS or via finite differences)
6. Create a fully-functional in-memory SpinerEOS object

**Note on persistence**: The constructed SpinerEOS exists in memory and can be used immediately. While the underlying DataBox objects have `saveHDF()` methods, a convenience method to save the entire SpinerEOS to an HDF5 file in the correct SP5 format would be a useful future enhancement (see Future Enhancements section).

## New Data Structures

### SpinerTableGridParams

Located in: `singularity-eos/eos/eos_spiner_rho_sie.hpp` (or possibly a separate header)

```cpp
struct SpinerTableGridParams {
  // Density bounds
  Real rhoMin, rhoMax;
  int numRho = -1;  // -1 means use numRhoPerDecade
  int numRhoPerDecade = 350;  // PPD_DEFAULT_RHO from sesame2spiner
  Real shrinklRhoBounds = 0.0;  // shrink log bounds by fraction

  // Temperature bounds
  Real TMin, TMax;
  int numT = -1;  // -1 means use numTPerDecade
  int numTPerDecade = 100;  // PPD_DEFAULT_T
  Real shrinklTBounds = 0.0;

  // SIE bounds (energy can be negative!)
  Real sieMin, sieMax;
  int numSie = -1;  // -1 means use numSiePerDecade
  int numSiePerDecade = 100;
  Real shrinkleBounds = 0.0;

  // Offset control (usually automatic, but allow override)
  // Set to -1 for auto-compute (default behavior)
  Real rhoOffset = -1.0;
  Real TOffset = -1.0;
  Real sieOffset = -1.0;

  // Enforce positive minimums (like sesame2spiner does for rho/T)
  // Set to <= 0 to disable enforcement
  Real strictlyPositiveMinRho = 1e-8;   // STRICTLY_POS_MIN_RHO
  Real strictlyPositiveMinT = 1e-2;     // STRICTLY_POS_MIN_T
  Real strictlyPositiveMinSie = -1.0;   // disabled for sie (can be negative)

  // Material properties
  int matid = 0;
  Real Abar = std::numeric_limits<Real>::signaling_NaN();
  Real Zbar = std::numeric_limits<Real>::signaling_NaN();
  Real rhoNormal = std::numeric_limits<Real>::signaling_NaN();

  // Piecewise grid options (advanced - follow sesame2spiner defaults)
  bool piecewiseRho = true;
  bool piecewiseT = true;
  bool piecewiseSie = true;
  Real rhoCoarseFactorLo = 3.0;   // COARSE_FACTOR_DEFAULT_RHO_LO
  Real rhoCoarseFactorHi = 5.0;   // COARSE_FACTOR_DEFAULT_RHO_HI
  Real TCoarseFactor = 1.5;       // COARSE_FACTOR_DEFAULT_T
  Real sieCoarseFactor = 1.5;
  Real rhoFineDiameterDecades = 1.5;  // RHO_FINE_DIAMETER_DEFAULT
  Real TSplitPoint = 1e4;             // T_SPLIT_POINT_DEFAULT

  // Optional: fine grid bounds override (advanced use)
  Real rhoFineMin = -1.0;  // -1 means use diameter
  Real rhoFineMax = -1.0;
};
```

## Constructor Signature

```cpp
template <template <class> class TransformerT>
template <typename EOS>
inline SpinerEOSDependsRhoSieTransformable<TransformerT>::SpinerEOSDependsRhoSieTransformable(
    const EOS& source_eos,
    const SpinerTableGridParams& params,
    bool reproducibility_mode = false);
```

**Note**: This constructor only supports `TableSplit::Total`. Support for `ElectronOnly` and `IonCold` splits requires additional physics that generic EOS may not provide.

## Source EOS Requirements

The source EOS must provide these methods (matching singularity EOS interface):

**Required:**
- `Real TemperatureFromDensityInternalEnergy(Real rho, Real sie, [lambda])`
- `Real InternalEnergyFromDensityTemperature(Real rho, Real T, [lambda])`
- `Real PressureFromDensityTemperature(Real rho, Real T, [lambda])`
- `Real PressureFromDensityInternalEnergy(Real rho, Real sie, [lambda])`

**Optional (will finite-difference if not available):**
All derivatives and bulk modulus are **optional**. The constructor will:
1. Use SFINAE to detect if methods exist on the source EOS
2. If available, call them directly for better accuracy
3. If not available, compute via finite differences

Optional methods:
- `Real BulkModulusFromDensityTemperature(Real rho, Real T, [lambda])`
- `Real BulkModulusFromDensityInternalEnergy(Real rho, Real sie, [lambda])`
- Derivatives: `dPdRho`, `dPdE`, `dTdRho`, `dTdE`, `dEdRho`, `dEdT`

**For material properties (if not in params):**
- `Real Abar()` or via params (if not provided, will use NaN)
- `Real Zbar()` or via params (if not provided, will use NaN)

## Implementation Steps

### 1. Parameter Processing

- Apply `strictlyPositiveMin*` enforcement to bounds
- Create `Bounds` objects for rho, T, and sie
  - Bounds class will auto-compute offsets if min ≤ 0
  - Apply shrinkRange parameters
  - Handle piecewise grid construction
- Override offsets if user specified (rhoOffset/TOffset/sieOffset ≥ 0)
- Determine reference state (rhoNormal, TNormal)

### 2. Grid Construction

Following sesame2spiner logic in `getMatBounds()`:

```cpp
// Apply strictly positive minimums
if (params.strictlyPositiveMinRho > 0 && params.rhoMin < params.strictlyPositiveMinRho) {
  params.rhoMin = params.strictlyPositiveMinRho;
}
// ... similar for T and sie

// Determine number of points
int numRho = params.numRho;
if (numRho <= 0) {
  numRho = Bounds::getNumPointsFromPPD(rhoMin, rhoMax, params.numRhoPerDecade);
}
// ... similar for T and sie

// Construct Bounds with piecewise grids if enabled
if (params.piecewiseRho) {
  lRhoBounds = Bounds(Bounds::ThreeGrids(), rhoMin, rhoMax,
                      rhoNormal, params.rhoFineDiameterDecades,
                      params.numRhoPerDecade, params.rhoCoarseFactorLo,
                      params.rhoCoarseFactorHi, true, params.shrinklRhoBounds);
} else {
  lRhoBounds = Bounds(rhoMin, rhoMax, numRho, true, params.shrinklRhoBounds);
}
// ... similar for T and sie
```

### 3. Table Population

For **dependsRhoT** tables (sie as function of rho, T):

```cpp
// Allocate tables
sie_.resize(lRhoBounds.grid, lTBounds.grid);
dependsRhoT_.P.resize(lRhoBounds.grid, lTBounds.grid);
// ... other fields

// Loop over grid
for (int j = 0; j < lRhoBounds.grid.nPoints(); j++) {
  Real lRho = lRhoBounds.grid.x(j);
  Real rho = singularity::FastMath::pow10(lRho) - lRhoBounds.offset;

  for (int i = 0; i < lTBounds.grid.nPoints(); i++) {
    Real lT = lTBounds.grid.x(i);
    Real T = singularity::FastMath::pow10(lT) - lTBounds.offset;

    // Evaluate source EOS
    sie_(j, i) = source_eos.InternalEnergyFromDensityTemperature(rho, T);
    dependsRhoT_.P(j, i) = source_eos.PressureFromDensityTemperature(rho, T);

    // Bulk modulus (or compute from derivatives)
    if constexpr (has_bulk_modulus_method<EOS>) {
      dependsRhoT_.bMod(j, i) = source_eos.BulkModulusFromDensityTemperature(rho, T);
    } else {
      // Compute from derivatives or finite difference
    }

    // ... populate other fields
  }
}
```

For **dependsRhoSie** tables (T as function of rho, sie):

```cpp
// Similar loop structure but over rho and sie grid
for (int j = 0; j < lRhoBounds.grid.nPoints(); j++) {
  for (int i = 0; i < leBounds.grid.nPoints(); i++) {
    Real rho = ...;
    Real sie = ...;

    T_(j, i) = source_eos.TemperatureFromDensityInternalEnergy(rho, sie);
    dependsRhoSie_.P(j, i) = source_eos.PressureFromDensityInternalEnergy(rho, sie);
    // ... other fields
  }
}
```

### 4. Derivative Computation

Strategy:
1. Check if source EOS has derivative methods (using SFINAE or `if constexpr`)
2. If available, call them directly
3. Otherwise, compute via finite differences

```cpp
// Example: dPdRho at constant T
Real dPdRho;
if constexpr (has_dPdRho_method<EOS>) {
  dPdRho = source_eos.dPdRhoFromDensityTemperature(rho, T);
} else {
  // Finite difference
  Real h = rho * 1e-6;  // relative step
  Real P_plus = source_eos.PressureFromDensityTemperature(rho + h, T);
  Real P_minus = source_eos.PressureFromDensityTemperature(rho - h, T);
  dPdRho = (P_plus - P_minus) / (2 * h);
}
```

### 5. Cold Curve Population

Compute cold curves (minimum sie at each rho):

```cpp
PCold_.resize(lRhoBounds.grid);
sieCold_.resize(lRhoBounds.grid);
// ...

for (int j = 0; j < lRhoBounds.grid.nPoints(); j++) {
  Real rho = ...;
  Real Tmin = ...; // minimum temperature

  sieCold_(j) = source_eos.InternalEnergyFromDensityTemperature(rho, Tmin);
  PCold_(j) = source_eos.PressureFromDensityTemperature(rho, Tmin);
  // ... derivatives
}
```

### 6. Bulk Modulus Fix-up

Call existing `calcBMod_()` method to ensure bulk modulus consistency (lines 509-510, 556-577).

### 7. Metadata and Reference State

```cpp
matid_ = params.matid;
split_ = split;
reproducible_ = reproducibility_mode;

// Material properties
AZbar_.Abar = std::isnan(params.Abar) ? source_eos.Abar() : params.Abar;
AZbar_.Zbar = std::isnan(params.Zbar) ? source_eos.Zbar() : params.Zbar;

// Reference state
rhoNormal_ = std::isnan(params.rhoNormal) ? /* compute sensible default */ : params.rhoNormal;
TNormal_ = ROOM_TEMPERATURE; // or from params
// ... compute sieNormal_, PNormal_, etc. at reference state

// Initialize transformer
InitializeTransformer();
```

### 8. Finalization

- Set memory status: `memoryStatus_ = DataStatus::OnHost`
- Compute PMin and rho_at_pmin curve
- Set up extrapolation data (PlRhoMax_, dPdRhoMax_)

## Helper Functions

### Method Detection (C++17 with C++20 migration path)

We'll use a clean detection idiom that's readable in C++17 and can be easily replaced with concepts in C++20.

**Strategy**: Keep all detection logic in a separate namespace that can be swapped out entirely when moving to C++20.

```cpp
// In a separate namespace for easy C++20 migration
namespace eos_builder {

// C++17 implementation using detection idiom
#if __cplusplus < 202002L

// Detection helper
template <typename...>
using void_t = void;

// Detect if EOS::BulkModulusFromDensityTemperature exists
template <typename EOS, typename = void>
struct has_bmod_rho_T : std::false_type {};

template <typename EOS>
struct has_bmod_rho_T<EOS, void_t<decltype(
    std::declval<EOS>().BulkModulusFromDensityTemperature(
        std::declval<Real>(), std::declval<Real>()))>> : std::true_type {};

// Similar for other methods...
template <typename EOS, typename = void>
struct has_bmod_rho_sie : std::false_type {};

template <typename EOS>
struct has_bmod_rho_sie<EOS, void_t<decltype(
    std::declval<EOS>().BulkModulusFromDensityInternalEnergy(
        std::declval<Real>(), std::declval<Real>()))>> : std::true_type {};

// Convenience constexpr bools
template <typename EOS>
inline constexpr bool has_bmod_rho_T_v = has_bmod_rho_T<EOS>::value;

template <typename EOS>
inline constexpr bool has_bmod_rho_sie_v = has_bmod_rho_sie<EOS>::value;

#else
// C++20 implementation using concepts (for future migration)
// This block can replace the above when moving to C++20

template <typename EOS>
concept has_bmod_rho_T = requires(EOS eos, Real rho, Real T) {
    { eos.BulkModulusFromDensityTemperature(rho, T) } -> std::same_as<Real>;
};

template <typename EOS>
concept has_bmod_rho_sie = requires(EOS eos, Real rho, Real sie) {
    { eos.BulkModulusFromDensityInternalEnergy(rho, sie) } -> std::same_as<Real>;
};

// For compatibility with C++17 code, provide _v helpers
template <typename EOS>
inline constexpr bool has_bmod_rho_T_v = has_bmod_rho_T<EOS>;

template <typename EOS>
inline constexpr bool has_bmod_rho_sie_v = has_bmod_rho_sie<EOS>;

#endif

} // namespace eos_builder
```

**Usage in constructor code**:

```cpp
// Clean, readable usage that works in both C++17 and C++20
if constexpr (eos_builder::has_bmod_rho_T_v<EOS>) {
  dependsRhoT_.bMod(j, i) = source_eos.BulkModulusFromDensityTemperature(rho, T);
} else {
  // Finite difference or compute from derivatives
  dependsRhoT_.bMod(j, i) = computeBulkModulus(source_eos, rho, T);
}
```

This approach:
- ✅ Readable in C++17 (just `if constexpr` with a clear name)
- ✅ Easy C++20 migration (replace one preprocessor block)
- ✅ All detection logic isolated in `eos_builder` namespace
- ✅ Consistent interface (`_v` helpers) across both standards

### Finite Difference Helper

```cpp
template <typename Func>
Real finiteDifference(Func f, Real x, Real fx) {
  Real h = std::max(std::abs(x) * 1e-6, 1e-12);
  Real f_plus = f(x + h);
  Real f_minus = f(x - h);
  return (f_plus - f_minus) / (2 * h);
}
```

## Testing Strategy

### Primary Test: IdealGas

**Focus on IdealGas EOS as the primary test case:**

1. **Basic construction test**:
   - Create IdealGas EOS with known parameters (Cv, gm1)
   - Build Spiner table from it with default grid parameters
   - Compare table values to direct IdealGas evaluation at grid points
   - Verify derivatives match analytic values

2. **Interpolation accuracy test**:
   - Sample random points within table bounds
   - Compare SpinerEOS interpolated values to IdealGas direct evaluation
   - Check that interpolation error is within acceptable tolerance

3. **Grid parameter variations**:
   - Test with different numRhoPerDecade, numTPerDecade, numSiePerDecade
   - Test piecewise vs uniform grids
   - Test shrinkRange parameters

4. **Negative energy handling** (if applicable):
   - Use shifted IdealGas if energy can be negative in test range
   - Verify offset is computed and applied correctly

### Additional Unit Tests

1. **Finite difference vs analytic derivatives**:
   - IdealGas has analytic derivatives - verify we use them
   - Create a minimal test EOS without derivative methods
   - Verify finite differences are used and give reasonable results

2. **Edge cases**:
   - Very narrow grids (few points)
   - Very wide grids (many decades)
   - Minimum strictly positive enforcement

### Integration Tests (Future)

1. Round-trip test: IdealGas → Spiner → Interpolate → Compare to original
2. Test GetOnDevice() works with constructed tables
3. Test with other analytic EOS (Gruneisen, JWL) once IdealGas works

## File Locations

- **Header**: `singularity-eos/eos/eos_spiner_rho_sie.hpp`
  - Add `SpinerTableGridParams` struct
  - Add templated constructor
  - Add private helper methods

- **Test**: `test/test_eos_spiner_construction.cpp` (new file)
  - Unit tests for constructor
  - Test various grid parameters
  - Test with different source EOS types

## Potential Issues & Solutions

### Issue 1: Template Instantiation
- **Problem**: Template constructor might cause compilation issues
- **Solution**: Use SFINAE carefully, provide good error messages

### Issue 2: Memory Allocation
- **Problem**: Large tables could exhaust memory
- **Solution**: Add sanity checks on grid sizes, warn if table is huge

### Issue 3: Derivative Accuracy
- **Problem**: Finite differences might be inaccurate
- **Solution**: Use adaptive step size, or require derivatives

### Issue 4: Cold Curve Definition
- **Problem**: What is "cold" for an arbitrary EOS?
- **Solution**: Use minimum temperature from grid bounds

### Issue 5: TableSplit for Generic EOS
- **Problem**: Generic EOS might not support electron/ion split
- **Solution**: Only support `TableSplit::Total` for this constructor (documented in signature)
- **Status**: Decided - only Total split supported

## Future Enhancements

1. **HDF5 saving capability**: Add a `SaveToHDF5(filename)` method to write constructed table to file
   - Would mirror sesame2spiner's `saveMaterial()` function structure
   - Save all DataBox tables with proper SP5 format
   - Include material metadata and offsets as attributes

2. **TableSplit support**: Extend to ElectronOnly and IonCold splits
   - Requires source EOS to provide electron/ion contributions separately
   - More complex interface requirements

3. **Mass fractions**: Support for multiphase materials
   - Requires source EOS to provide mass fraction data

4. **Parallel construction**: OpenMP over grid points for large tables

5. **Adaptive grid refinement**: Automatically refine grid in regions with high curvature

6. **StellarCollapse support**: Constructor for StellarCollapse-style EOS with Ye dependence

## Design Decisions (Finalized)

1. ✅ **TableSplit**: Only `TableSplit::Total` supported initially
2. ✅ **HDF5 save**: Not in initial constructor, but future enhancement
3. ✅ **Derivatives**: Optional - use if available, otherwise finite-difference
4. ✅ **Test EOS**: Primary testing with IdealGas
