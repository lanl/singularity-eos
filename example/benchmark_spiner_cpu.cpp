//------------------------------------------------------------------------------
// © 2021-2025. Triad National Security, LLC. All rights reserved.  This
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
  This script compares the performance of EOSPAC SpinerEOSDependsRhoT and SpinerEOSDependsRhoSie in the following cases:

  1. P(rho, T) 
  2. e(rho, T) 
  3. P(rho, e(rho, T)) 
  4. T(rho, e(rho, T)) 

  The output of this script is as follows:
  1. For each material, model, and lookup type, a .csv of the grid values (to compare accuracy)
  2. A .csv of the lookup times for x-amount of trials (to compare performance)

  Example usage: ./Benchmark ./output_dir 100 100 ./tables/my_eos.sp5 2020 3720

  Authors: Erin O'Neil and Joshua Basabe
 */

 //C++ Headers
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <fstream>   
#include <iomanip> 


//Data headers
#include <hdf5.h>
#include <hdf5_hl.h>

//Spiner headers
#include <ports-of-call/portability.hpp>
#include <singularity-eos/base/sp5/singularity_eos_sp5.hpp>
#include <spiner/databox.hpp>
#include <spiner/interpolation.hpp>
#include <spiner/sp5.hpp>
#include <spiner/spiner_types.hpp>
#include <ports-of-call/portability.hpp>

//Get EOS Models
#include <singularity-eos/eos/eos.hpp>
#include "singularity-eos/eos/eos_spiner.hpp"
#include <eospac-wrapper/eospac_wrapper.hpp>

using namespace singularity;
using namespace EospacWrapper;
using namespace std::chrono;
#include <experimental/filesystem> //not sure why I had to do this
namespace fs = std::experimental::filesystem;
using DataBox = Spiner::DataBox<Real>;
using RegularGrid1D = Spiner::RegularGrid1D<Real>;
using Real = double;

// == Creates the Bounds of the Grid ==
class Bounds { //check for edge effects
 public:
  Bounds(Real min, Real max, int N) : offset(0.0) {
    constexpr Real epsilon = std::numeric_limits<float>::epsilon();
    const Real min_offset = 10 * std::abs(epsilon);
    if (min <= 1e-16) offset = 1.1 * std::abs(min) + min_offset; //changed this from (min <= 0) because small values will make the log blow up
    min += offset;
    max += offset;
    min = std::log10(min);
    max = std::log10(max);
    dx = (max - min) / (N - 1);
    x0 = min;
    size = N;
  }
    Real i2lin(int i) const {
    return std::pow(10.0, x0 + static_cast<Real>(i) * dx) - offset; 
  }
 private:
  Real x0, dx, offset;
  int size; 
};

// == The following function creates and fills in a csv of grid values ==
void write_matrix_csv(const std::string& filename, const std::vector<Real>& flat_matrix, size_t nRows, size_t nCols) {
    std::ofstream file(filename);
    file << std::setprecision(14);
    for (size_t i = 0; i < nRows; ++i) {
        for (size_t j = 0; j < nCols; ++j) {
            file << flat_matrix[i * nCols + j];
            if (j < nCols - 1) file << ",";
        }
        file << "\n";
    }
}

// == The following saves the rho and T values corresponding to the indices (for plotting purposes), as well as time trials ==
void write_vector_csv(const std::string& filename, const std::vector<Real>& vec) {
    std::ofstream file(filename);
    file << std::setprecision(14);
    for (const auto& val : vec)
        file << val << "\n";
}



// == Main loop ===
int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " path/to/output nRho nT sp5_file matid1 matid2 ..."; //Is this good usage? 
        return 1;
    }
    fs::path base_output_path(argv[1]);
    //if the folder does not exist, try to make one
    if (!fs::exists(base_output_path)) {
        if (!fs::create_directories(base_output_path)) {
            std::cerr << "Failed to create output directory: " << base_output_path << "\n";
            return 1; }}
    int nRho = std::stoi(argv[2]);  //test on batch sizes of 144, 269, 512, 2048 or round powers of ten
    int nT = std::stoi(argv[3]); 
    std::string sp5file = argv[4];

    std::vector<int> matids;
    for (int i = 5; i < argc; ++i)
        matids.push_back(std::atoi(argv[i]));

    // == Iterate through each material ==
    for (int matid : matids) {
        const int ntimes = 20; // Perform 20 time trials

        // == Get material metadata ==
        SesameMetadata meta;
        eosGetMetadata(matid, meta);

        // == Set up bounds ==
        Real rhoMin = 1.1 * std::max(meta.rhoMin, 1e-5);
        Real rhoMax = 0.9 * meta.rhoMax;
        Real TMin   = 1.1 * std::max(meta.TMin, 1.0);
        Real TMax   = 0.9 * meta.TMax;
        Bounds rhoBounds(rhoMin, rhoMax, nRho);
        Bounds TBounds(TMin, TMax, nT);

        std::vector<Real> rhos(nRho), temps(nT);
        for (int i = 0; i < nRho; ++i) rhos[i] = rhoBounds.i2lin(i);
        for (int j = 0; j < nT; ++j) temps[j] = TBounds.i2lin(j);

        // == This will save the values of rho and T corresponding to the indices chosen ==
        std::string rho_file = (base_output_path / ("rho_axis_" + std::to_string(matid) + "_nRho-" + std::to_string(nRho) + ".csv")).string();
        std::string temp_file = (base_output_path / ("temp_axis_" + std::to_string(matid) + "_nT-" + std::to_string(nT) + ".csv")).string();
        write_vector_csv(rho_file, rhos);
        write_vector_csv(temp_file, temps);

        // === Load EOS Models ===
        SpinerEOSDependsRhoT eos_rt(sp5file, matid);
        SpinerEOSDependsRhoSie eos_rs(sp5file, matid); //eventually make it SpinerEOSDependsRhoSie<NullTransfom>, for example
	    EOSPAC eos_ref(matid); 

        //These vectors will store the compute time for each model for the 20 trials
        std::vector<double> time_sie_eos_ref_list, time_sie_rt_list, time_sie_rs_list;
        std::vector<double> time_P_eos_ref_T_list, time_P_rt_T_list, time_P_rs_T_list;
        std::vector<double> time_P_eos_ref_sie_list, time_P_rt_sie_list, time_P_rs_sie_list;
        std::vector<double> time_T_back_eos_ref_list, time_T_back_rt_list, time_T_back_rs_list;

        //These matrices will be filled with the EOS values (overwriting for each trail)
        //Then at the end of the program, we export the matrix to a .csv to later use for accuracy tests
        size_t grid_size = nRho * nT;

        // Internal energy
        std::vector<Real> sie_eos_ref(grid_size, 0.0);
        std::vector<Real> sie_rt(grid_size, 0.0);
        std::vector<Real> sie_rs(grid_size, 0.0);

        // Pressure from (rho, T)
        std::vector<Real> P_eos_ref_T(grid_size, 0.0);
        std::vector<Real> P_rt_T(grid_size, 0.0);
        std::vector<Real> P_rs_T(grid_size, 0.0);

        // Pressure from (rho, sie)
        std::vector<Real> P_eos_ref_sie(grid_size, 0.0);
        std::vector<Real> P_rt_sie(grid_size, 0.0);
        std::vector<Real> P_rs_sie(grid_size, 0.0);

        // Temperature from (rho, sie)
        std::vector<Real> T_back_eos_ref(grid_size, 0.0);
        std::vector<Real> T_back_rt(grid_size, 0.0);
        std::vector<Real> T_back_rs(grid_size, 0.0);


        // === Benchmark Loop ===
        for (int rep = 0; rep < ntimes; ++rep) {
                // Clear matrices
                std::fill(sie_eos_ref.begin(), sie_eos_ref.end(), 0.0);
                std::fill(sie_rt.begin(), sie_rt.end(), 0.0);
                std::fill(sie_rs.begin(), sie_rs.end(), 0.0);

                std::fill(P_eos_ref_T.begin(), P_eos_ref_T.end(), 0.0);
                std::fill(P_rt_T.begin(), P_rt_T.end(), 0.0);
                std::fill(P_rs_T.begin(), P_rs_T.end(), 0.0);

                std::fill(P_eos_ref_sie.begin(), P_eos_ref_sie.end(), 0.0);
                std::fill(P_rt_sie.begin(), P_rt_sie.end(), 0.0);
                std::fill(P_rs_sie.begin(), P_rs_sie.end(), 0.0);

                std::fill(T_back_eos_ref.begin(), T_back_eos_ref.end(), 0.0);
                std::fill(T_back_rt.begin(), T_back_rt.end(), 0.0);
                std::fill(T_back_rs.begin(), T_back_rs.end(), 0.0);



                // == EOSPAC (_eos) model ==
                auto t0 = high_resolution_clock::now(); //e(rho, T) 
                for (int i = 0; i < nRho; ++i) { //replace these with a single for loop that loops through (i,j)
                    for (int j = 0; j < nT; ++j) { 
                        sie_eos_ref[i * nT + j] = eos_ref.InternalEnergyFromDensityTemperature(rhos[i], temps[j]); } }
                auto t1 = high_resolution_clock::now();
                std::chrono::duration<double, std::micro> elapsed = t1 - t0; //define "elapsed"
                time_sie_eos_ref_list.push_back(elapsed.count());

                t0 = high_resolution_clock::now(); // P(rho, T)
                for (int i = 0; i < nRho; ++i) {
                    for (int j = 0; j < nT; ++j) {
                        P_eos_ref_T[i * nT + j] = eos_ref.PressureFromDensityTemperature(rhos[i], temps[j]);}}
                t1 = high_resolution_clock::now();
                elapsed = t1 - t0;
                time_P_eos_ref_T_list.push_back(elapsed.count());

                t0 = high_resolution_clock::now(); // P(rho, sie)
                for (int i = 0; i < nRho; ++i) {
                    for (int j = 0; j < nT; ++j) { //is this line supposed to still be nT, it uses sie_eos_ref(i,j)
                        P_eos_ref_sie[i * nT + j] = eos_ref.PressureFromDensityInternalEnergy(rhos[i], sie_eos_ref[i * nT + j]);}}
                t1 = high_resolution_clock::now();
                elapsed = t1 - t0;
                time_P_eos_ref_sie_list.push_back(elapsed.count());

                t0 = high_resolution_clock::now(); //T(rho, sie)
                for (int i = 0; i < nRho; ++i) {
                    for (int j = 0; j < nT; ++j) { //same with this one
                        T_back_eos_ref[i * nT + j] = eos_ref.TemperatureFromDensityInternalEnergy(rhos[i], sie_eos_ref[i * nT + j]);}}
                t1 = high_resolution_clock::now();
                elapsed = t1 - t0;
                time_T_back_eos_ref_list.push_back(elapsed.count());



                // == SpinerEOSDependsRhoT (_rt) model ==
                t0 = high_resolution_clock::now(); //e(rho, T) 
                for (int i = 0; i < nRho; ++i) {
                    for (int j = 0; j < nT; ++j) {
                        sie_rt[i * nT + j] = eos_rt.InternalEnergyFromDensityTemperature(rhos[i], temps[j]);}}
                t1 = high_resolution_clock::now();
                elapsed = t1 - t0;
                time_sie_rt_list.push_back(elapsed.count());

	            t0 = high_resolution_clock::now(); // P(rho, T)
                for (int i = 0; i < nRho; ++i) {
                    for (int j = 0; j < nT; ++j) {
                        P_rt_T[i * nT + j] = eos_rt.PressureFromDensityTemperature(rhos[i], temps[j]);}}
                t1 = high_resolution_clock::now();
                elapsed = t1 - t0;
                time_P_rt_T_list.push_back(elapsed.count());

                t0 = high_resolution_clock::now(); // P(rho, sie)
                for (int i = 0; i < nRho; ++i) {
                    for (int j = 0; j < nT; ++j) { //same with this one
                        P_rt_sie[i * nT + j] = eos_rt.PressureFromDensityInternalEnergy(rhos[i], sie_rt[i * nT + j]);}}
                t1 = high_resolution_clock::now();
                elapsed = t1 - t0;
                time_P_rt_sie_list.push_back(elapsed.count());

                t0 = high_resolution_clock::now(); //T(rho, sie)
                for (int i = 0; i < nRho; ++i) {
                    for (int j = 0; j < nT; ++j) { //same with this one
                        T_back_rt[i * nT + j] = eos_rt.TemperatureFromDensityInternalEnergy(rhos[i], sie_rt[i * nT + j]);}}
                t1 = high_resolution_clock::now();
                elapsed = t1 - t0;
                time_T_back_rt_list.push_back(elapsed.count());



                // == SpinerEOSDependsRhoSie (_rs) model ==
                t0 = high_resolution_clock::now(); //e(rho, T) 
                for (int i = 0; i < nRho; ++i) {
                    for (int j = 0; j < nT; ++j) {
                        sie_rs[i * nT + j] = eos_rs.InternalEnergyFromDensityTemperature(rhos[i], temps[j]);}}
                t1 = high_resolution_clock::now();
                elapsed = t1 - t0;
                time_sie_rs_list.push_back(elapsed.count());

		        t0 = high_resolution_clock::now(); // P(rho, T)
                for (int i = 0; i < nRho; ++i) {
                    for (int j = 0; j < nT; ++j) {
                        P_rs_T[i * nT + j] = eos_rs.PressureFromDensityTemperature(rhos[i], temps[j]);}}
                t1 = high_resolution_clock::now();
                elapsed = t1 - t0;
                time_P_rs_T_list.push_back(elapsed.count());

                t0 = high_resolution_clock::now(); // P(rho, sie)
                for (int i = 0; i < nRho; ++i) {
                    for (int j = 0; j < nT; ++j) { //same with this one
                        P_rs_sie[i * nT + j] = eos_rs.PressureFromDensityInternalEnergy(rhos[i], sie_rs[i * nT + j]);}}
                t1 = high_resolution_clock::now();
                elapsed = t1 - t0;
                time_P_rs_sie_list.push_back(elapsed.count());

                t0 = high_resolution_clock::now(); //T(rho, sie)
                for (int i = 0; i < nRho; ++i) {
                    for (int j = 0; j < nT; ++j) { //same with this one
                        T_back_rs[i * nT + j] = eos_rs.TemperatureFromDensityInternalEnergy(rhos[i], sie_rs[i * nT + j]);}}
                t1 = high_resolution_clock::now();
                elapsed = t1 - t0;
                time_T_back_rs_list.push_back(elapsed.count());

        } //end for loop for # of trials

    //Only save the filled in grid values on the last loop (only one grid exported to compare accuracies)
    std::string filename;

    //internal energy grids
    filename = "sie_eos_ref_" + std::to_string(matid) + "_nRho-" + std::to_string(nRho) + "_nT-" + std::to_string(nT) + ".csv";
    write_matrix_csv((base_output_path / filename).string(), sie_eos_ref, nRho, nT);

    filename = "sie_rt_" + std::to_string(matid) + "_nRho-" + std::to_string(nRho) + "_nT-" + std::to_string(nT) + ".csv";
    write_matrix_csv((base_output_path / filename).string(), sie_rt, nRho, nT);

    filename = "sie_rs_" + std::to_string(matid) + "_nRho-" + std::to_string(nRho) + "_nT-" + std::to_string(nT) + ".csv";
    write_matrix_csv((base_output_path / filename).string(), sie_rs, nRho, nT);



    //pressure from density & temperature grids
    filename = "P_eos_ref_T_" + std::to_string(matid) + "_nRho-" + std::to_string(nRho) + "_nT-" + std::to_string(nT) + ".csv";
    write_matrix_csv((base_output_path / filename).string(), P_eos_ref_T, nRho, nT);

    filename = "P_rt_T_" + std::to_string(matid) + "_nRho-" + std::to_string(nRho) + "_nT-" + std::to_string(nT) + ".csv";
    write_matrix_csv((base_output_path / filename).string(), P_rt_T, nRho, nT);

    filename = "P_rs_T_" + std::to_string(matid) + "_nRho-" + std::to_string(nRho) + "_nT-" + std::to_string(nT) + ".csv";
    write_matrix_csv((base_output_path / filename).string(), P_rs_T, nRho, nT);



    //pressure from density & internal energy grids
    filename = "P_eos_ref_sie_" + std::to_string(matid) + "_nRho-" + std::to_string(nRho) + "_nT-" + std::to_string(nT) + ".csv";
    write_matrix_csv((base_output_path / filename).string(), P_eos_ref_sie, nRho, nT);

    filename = "P_rt_sie_" + std::to_string(matid) + "_nRho-" + std::to_string(nRho) + "_nT-" + std::to_string(nT) + ".csv";
    write_matrix_csv((base_output_path / filename).string(), P_rt_sie, nRho, nT);

    filename = "P_rs_sie_" + std::to_string(matid) + "_nRho-" + std::to_string(nRho) + "_nT-" + std::to_string(nT) + ".csv";
    write_matrix_csv((base_output_path / filename).string(), P_rs_sie, nRho, nT);



    //temperature grids
    filename = "T_back_eos_ref_" + std::to_string(matid) + "_nRho-" + std::to_string(nRho) + "_nT-" + std::to_string(nT) + ".csv";
    write_matrix_csv((base_output_path / filename).string(), T_back_eos_ref, nRho, nT);

    filename = "T_back_rt_" + std::to_string(matid) + "_nRho-" + std::to_string(nRho) + "_nT-" + std::to_string(nT) + ".csv";
    write_matrix_csv((base_output_path / filename).string(), T_back_rt, nRho, nT);

    filename = "T_back_rs_" + std::to_string(matid) + "_nRho-" + std::to_string(nRho) + "_nT-" + std::to_string(nT) + ".csv";
    write_matrix_csv((base_output_path / filename).string(), T_back_rs, nRho, nT);


    //after creating this .csv, should I average the last 15, can you explain again why the first 5 may be bad
    std::string timing_file = (base_output_path / ("timing_" + std::to_string(matid) + "_nRho-" + std::to_string(nRho) + "_nT-" + std::to_string(nT) + ".csv")).string();
    std::ofstream timing_out(timing_file);

    timing_out << "sie_eos_ref,sie_rt,sie_rs,P_eos_ref_T,P_rt_T,P_rs_T,P_eos_ref_sie,P_rt_sie,P_rs_sie,T_back_eos_ref,T_back_rt,T_back_rs\n";

    for (size_t i = 0; i < time_sie_eos_ref_list.size(); ++i) {
        timing_out << time_sie_eos_ref_list[i] << ","
                << time_sie_rt_list[i] << ","
                << time_sie_rs_list[i] << ","
                << time_P_eos_ref_T_list[i] << ","
                << time_P_rt_T_list[i] << ","
                << time_P_rs_T_list[i] << ","
                << time_P_eos_ref_sie_list[i] << ","
                << time_P_rt_sie_list[i] << ","
                << time_P_rs_sie_list[i] << ","
                << time_T_back_eos_ref_list[i] << ","
                << time_T_back_rt_list[i] << ","
                << time_T_back_rs_list[i] << "\n";
    }

    std::cout << "Benchmark complete for material " << std::to_string(matid) << "\n";
    }
    return 0;
}