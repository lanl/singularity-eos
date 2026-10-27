//======================================================================
// sesame2spiner tool for converting eospac to spiner
// Author: Jonah Miller (jonahm@lanl.gov)
// © 2021-2026. Triad National Security, LLC. All rights reserved.  This
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
//======================================================================

#ifndef _SESAME2SPINER_GENERATE_FILES_HPP_
#define _SESAME2SPINER_GENERATE_FILES_HPP_

#include <initializer_list>
#include <string>
#include <vector>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <eospac-wrapper/eospac_wrapper.hpp>
#include <singularity-utils/spiner_params.hpp>

#include "io_eospac.hpp"
#include "parser.hpp"

using namespace EospacWrapper;
using singularity::table_utils::SpinerTableGridParams;

namespace sesame2spiner {

constexpr int PPD_DEFAULT_RHO = 350;
constexpr int PPD_DEFAULT_T = 100;
constexpr Real STRICTLY_POS_MIN_RHO = 1e-8;
constexpr Real STRICTLY_POS_MIN_T = 1e-2;
constexpr Real COARSE_FACTOR_DEFAULT_RHO_LO = 3;
constexpr Real COARSE_FACTOR_DEFAULT_RHO_HI = 5;
constexpr Real COARSE_FACTOR_DEFAULT_T = 1.5;
constexpr Real RHO_FINE_DIAMETER_DEFAULT = 1.5;
constexpr Real T_SPLIT_POINT_DEFAULT = 1e4;

herr_t saveMaterial(hid_t loc, const SesameMetadata &metadata, const Bounds &lRhoBounds,
                    const Bounds &lTBounds, const Bounds &leBounds,
                    const std::string &name, const bool addSubtables,
                    Verbosity eospacWarn = Verbosity::Quiet);
inline herr_t saveMaterial(hid_t loc, const SesameMetadata &metadata,
                           const Bounds &lRhoBounds, const Bounds &lTBounds,
                           const Bounds &leBounds, const std::string &name,
                           Verbosity eospacWarn = Verbosity::Quiet) {
  return saveMaterial(loc, metadata, lRhoBounds, lTBounds, leBounds, name, false,
                      eospacWarn);
}

// Write the sp5 root attributes (singularity version and log type). Call this once
// on a newly created file before adding any materials to it. Exposed for callers
// that manage the hid_t themselves; the savename-based overloads below do it for you.
herr_t writeSP5RootAttributes(hid_t file);

// Check that an existing sp5 file's log type matches the one this build was
// compiled with. Appending materials to a file written with a different log type
// would produce a file whose materials disagree with its log_type attribute, which
// the SpinerEOS constructors cannot detect. Called by the hid_t overload below.
herr_t checkSP5RootAttributes(hid_t file);

// Add materials to an already-open sp5 file. This is the core implementation; the
// overloads below all funnel into it. Safe to call repeatedly on the same file to
// build it up one material at a time -- duplicate matid and name detection query
// the file itself, so they remain correct across calls.
herr_t saveAllMaterials(hid_t file, const std::vector<int> &matids,
                        const std::vector<Params> &params, bool printMetadata,
                        Verbosity eospacWarn);

// Save a list of materials, with per-material parameter overrides, to a new file.
herr_t saveAllMaterials(const std::string &savename, const std::vector<int> &matids,
                        const std::vector<Params> &params, bool printMetadata,
                        Verbosity eospacWarn);

// Save a list of materials to a new file using standard sesame2spiner defaults.
herr_t saveAllMaterials(const std::string &savename, const std::vector<int> &matids,
                        bool printMetadata, Verbosity eospacWarn);

// Save all materials described by a list of input files to a new file.
herr_t saveAllMaterials(const std::string &savename,
                        const std::vector<std::string> &filenames, bool printMetadata,
                        Verbosity eospacWarn);

// Disambiguates a braced list of string literals, e.g.
//   saveAllMaterials(savename, {"air.dat", "steel.dat"}, false, warn);
// Without this, such a call is ambiguous against the std::vector<int> overload
// above: std::vector<int> has an iterator-pair constructor, and a pair of
// const char* satisfies it (as a range of char, which converts to int). An
// exact-match std::initializer_list parameter outranks both vector candidates.
inline herr_t saveAllMaterials(const std::string &savename,
                               std::initializer_list<const char *> filenames,
                               bool printMetadata, Verbosity eospacWarn) {
  return saveAllMaterials(savename,
                          std::vector<std::string>(filenames.begin(), filenames.end()),
                          printMetadata, eospacWarn);
}

herr_t saveTablesRhoSie(hid_t loc, int matid, TableSplit split, const Bounds &lRhoBounds,
                        const Bounds &leBounds, Verbosity eospacWarn = Verbosity::Quiet);
herr_t saveTablesRhoT(hid_t loc, int matid, TableSplit split, const Bounds &lRhoBounds,
                      const Bounds &lTBounds, Verbosity eospacWarn = Verbosity::Quiet);

void getMatBounds(int matid, const SesameMetadata &metadata, const Params &params,
                  Bounds &lRhoBounds, Bounds &lTBounds, Bounds &leBounds);

// Convert string-based Params to structured SpinerTableGridParams
// This allows sesame2spiner to use the same grid construction logic as the EOS
// constructors
SpinerTableGridParams paramsToGridParams(int matid, const SesameMetadata &metadata,
                                         const Params &params);

bool checkValInMatBounds(int matid, const std::string &name, Real val, Real vmin,
                         Real vmax);

// Sanity check metadata returned by eosGetMetadata. eosGetMetadata discards the
// eosSafeLoad error code and leaves its buffers zero-initialized on failure, so a
// matid that is absent from the sesame file comes back with all-zero bounds rather
// than an error. Left unchecked, those bounds reach log() and produce NaN grids.
bool checkMetadataValid(int matid, const SesameMetadata &metadata);

int getNumPointsFromPPD(Real min, Real max, int ppd);
} // namespace sesame2spiner

#endif // _SESAME2SPINER_GENERATE_FILES_HPP_
