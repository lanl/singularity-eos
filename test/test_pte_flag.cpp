//------------------------------------------------------------------------------
// © 2021. Triad National Security, LLC. All rights reserved.  This
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

#include <iostream>
//#include <vector>
#include <chrono>
#include <memory>
#include <stdlib.h>
#include <time.h>

#include <ports-of-call/portability.hpp>
#include <singularity-eos/closure/mixed_cell_models.hpp>
#include <singularity-eos/eos/bedroom_door.hpp>
#include <singularity-eos/eos/eos.hpp>

static constexpr const int NMAT{3};
static constexpr const int NTRIAL{3};

template <typename T>
using TempArray = std::unique_ptr<T[]>;
using RealArray = TempArray<Real>;
using IntArray = TempArray<int>;

using namespace singularity;

void set_state(RealArray &rho, RealArray &vfrac, RealArray &sie, BetterEOS *eos) {

  rho[0] = 8.93;
  rho[1] = 1.89;
  rho[2] = 2.5;

  Real vsum = 0.;
  for (int i = 0; i < NMAT; i++) {
    sie[i] = eos[i].InternalEnergyFromDensityTemperature(rho[i], 600.0);
    vfrac[i] = rand() / (1.0 * RAND_MAX);
    vsum += vfrac[i];
  }

  for (int i = 0; i < NMAT; i++)
    vfrac[i] *= 1.0 / vsum;
}

void compare_values(const int n, const RealArray &a, const RealArray &b, const char *s) {
  for (int i{0}; i < n; ++i) {
    Real diff = fabs(b[i] - a[i]);
    const Real divide = a[i] == 0.0 ? 1.0 : fabs(a[i]);
    printf("%s[%i]: %e %e %e\n", s, i, a[i], b[i], diff / divide);
  }
  return;
}

void print_values(const int n, const RealArray &a, const char *s) {
  for (int i{0}; i < n; ++i) {
    printf("%s[%i]: %e\n", s, i, a[i]);
  }
  return;
}

int main(int argc, char *argv[]) {

  Real *lambda[NMAT];
  BetterEOS *eos;

  init_bd_eos(NMAT, eos);
  eos[0] = Gruneisen(394000.0, 1.489, 0.0, 0.0, 2.02, 0.47, 8.93, 297.0, 1.0e6, 0.383e7);
  eos[1] = DavisReactants(1.890, 4.115e10, 1.0e6, 297.0, 1.8, 4.6, 0.34, 0.56, 0.0,
                          0.4265, 0.001074e10);
  eos[2] =
      DavisProducts(0.798311, 0.58, 1.35, 2.66182, 0.75419, 3.2e10, 0.001072e10, 0.0);

  TempArray<Real> rho(new Real[NMAT]);
  auto vfrac = std::make_unique<Real[]>(NMAT);
  auto sie = std::make_unique<Real[]>(NMAT);
  auto temp = std::make_unique<Real[]>(NMAT);
  auto press = std::make_unique<Real[]>(NMAT);

  auto vfrac_flag = std::make_unique<Real[]>(NMAT);
  auto sie_flag = std::make_unique<Real[]>(NMAT);
  auto mass_flag = std::make_unique<Real[]>(NMAT);
  auto mats = std::make_unique<int[]>(NMAT);

  srand(time(NULL));
  for (int i = 0; i < NMAT; i++)
    lambda[i] = nullptr;

  std::chrono::duration<double> sum_time = std::chrono::duration<double>::zero();
  int nsuccess = 0;
  for (int n = 0; n < NTRIAL; n++) {
    set_state(rho, vfrac, sie, eos);
    Real sie_tot = 0.0;
    Real mass_tot = 0.0;
    for (int i = 0; i < NMAT; i++) {
      // std::cout << i << " " << rho[i] << " " << vfrac[i] << " " << sie[i] << std::endl;
      mass_flag[i] = rho[i] * vfrac[i];
      mass_tot += rho[i] * vfrac[i];
      sie_tot += rho[i] * vfrac[i] * sie[i];
      sie_flag[i] = sie[i];
      vfrac_flag[i] = vfrac[i];
      mats[i] = i;
    }
    sie_tot /= mass_tot;
    printf("Input set #%i:\nSIE: %e V: 1.0\n", n + 1, sie_tot);
    print_values(NMAT, sie, "sie");
    print_values(NMAT, vfrac, "vfrac");
    print_values(NMAT, rho, "rho");
    print_values(NMAT, mass_flag, "mass");
    // std::cout << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    bool pte_base_success{
        PTE(NMAT, eos, 1.0, sie_tot, rho, vfrac, sie, temp, press, lambda)};
    bool pte_flag_success{pte_closure_flag_offset(NMAT, eos, 1.0, sie_tot, mats.get(),
                                                  mass_flag.get(), vfrac_flag.get(),
                                                  sie_flag.get(), lambda)};
    if (pte_base_success && pte_flag_success) {
      printf("both solvers succeeded\n");
    } else if (!pte_base_success && !pte_flag_success) {
      printf("both solvers failed\n");
    } else if (!pte_base_success && pte_flag_success) {
      printf("base failed and flag succeeded\n");
    } else if (pte_base_success && !pte_flag_success) {
      printf("base succeeded and flag failed\n");
    }
    printf("          PTE          flag         |PTE - flag|/|PTE| \n");
    compare_values(NMAT, vfrac, vfrac_flag, "vfrac");
    compare_values(NMAT, sie, sie_flag, "sie");
    auto stop = std::chrono::high_resolution_clock::now();
    sum_time += stop - start;
    // nsuccess += (success ? 1 : 0);
  }

  // std::cout << "PTE solver returned in " << n << " iterations:" << std::endl;
  // for (int i = 0 ; i < NMAT; i++) {
  //  std::cout << i << " " << rho[i] << " " << vfrac[i] << " " << sie[i] << " " << " " <<
  //  rho[i]*vfrac[i]/rho_tot << " " << press[i] << " " << temp[i] << std::endl;
  //}
  // std::cout << "Finished " << NTRIAL << " solves in " << sum_time.count() << " seconds"
  // << std::endl; std::cout << "Solves/second = " << NTRIAL/sum_time.count() <<
  // std::endl; std::cout << "Success: " << nsuccess << "   Failure: " << NTRIAL-nsuccess
  // << std::endl;

  finalize_bd_eos(NMAT, eos, 1);
  return 0;
}
