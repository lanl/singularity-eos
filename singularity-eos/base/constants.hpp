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

#ifndef SINGULARITY_EOS_BASE_CONSTANTS_HPP_
#define SINGULARITY_EOS_BASE_CONSTANTS_HPP_

namespace singularity {

namespace thermalqs {
constexpr unsigned long none = 0;
constexpr unsigned long density = (1 << 0);
constexpr unsigned long specific_internal_energy = (1 << 1);
constexpr unsigned long pressure = (1 << 2);
constexpr unsigned long temperature = (1 << 3);
constexpr unsigned long specific_heat = (1 << 4);
constexpr unsigned long bulk_modulus = (1 << 5);
constexpr unsigned long all_values = (1 << 6) - 1;
} // namespace thermalqs

constexpr size_t MAX_NUM_LAMBDAS = 3;
enum class DataStatus { Deallocated = 0, OnDevice = 1, OnHost = 2 };
enum class TableStatus { OnTable = 0, OffBottom = 1, OffTop = 2 };
constexpr Real ROOM_TEMPERATURE = 293; // K
constexpr Real ATMOSPHERIC_PRESSURE = 1e6;
constexpr Real EPS = 10.0 * std::numeric_limits<Real>::epsilon();

} // namespace singularity

#endif // SINGULARITY_EOS_BASE_CONSTANTS_HPP_
