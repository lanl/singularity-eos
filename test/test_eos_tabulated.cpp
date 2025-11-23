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

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream> // debug

#include <ports-of-call/portability.hpp>
#include <ports-of-call/portable_arrays.hpp>
#include <ports-of-call/portable_errors.hpp>
#include <singularity-eos/base/constants.hpp>
#include <singularity-eos/base/variadic_utils.hpp>
#include <singularity-eos/eos/eos.hpp>

#ifndef CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_FAST_COMPILE
#include <catch2/catch_test_macros.hpp>
#endif

#include <test/eos_unit_test_helpers.hpp>

#ifdef SPINER_USE_HDF
#include <singularity-eos/base/spiner_table_utils.hpp>
using singularity::SpinerEOSDependsRhoSie;
using singularity::SpinerEOSDependsRhoT;
#endif
using singularity::UnitSystem;

#ifdef SINGULARITY_USE_EOSPAC
using singularity::EOSPAC;
#endif

namespace thermalqs = singularity::thermalqs;
using singularity::variadic_utils::np;

const std::string eosName = "../materials.sp5";
const std::string airName = "air";
const std::string steelName = "stainless steel 347";

#ifdef SPINER_USE_HDF
#ifdef SINGULARITY_TEST_SESAME
constexpr int steelID = 4272;
constexpr int airID = 5030;
constexpr int DTID = 5267;
constexpr int gID = 2700;
constexpr int titaniumID = 2961;
constexpr Real ev2k = 1.160451812e4;
#endif // SINGULARITY_TEST_SESAME
#endif // SPINER_USE_HDF

#ifdef SPINER_USE_HDF
#ifdef SINGULARITY_TEST_SESAME
#ifdef SINGULARITY_USE_EOSPAC
using EOS = singularity::Variant<SpinerEOSDependsRhoSie, SpinerEOSDependsRhoT, EOSPAC>;

SCENARIO("SpinerEOS depends on Rho and T", "[SpinerEOS][DependsRhoT][EOSPAC]") {

  GIVEN("SpinerEOS and EOSPAC EOS for steel can be initialized with matid") {
    EOS steelEOS_host_polymorphic = SpinerEOSDependsRhoT(eosName, steelID);
    EOS steelEOS = steelEOS_host_polymorphic.GetOnDevice();
    auto steelEOS_host = steelEOS_host_polymorphic.get<SpinerEOSDependsRhoT>();

    EOS eospac = EOSPAC(steelID);

    THEN("The correct metadata is read in") {
      REQUIRE(steelEOS_host.matid() == steelID);

      AND_THEN("We can get a reference density and temperature") {
        Real rho, T, sie, P, cv, bmod, dpde, dvdt;
        Real rho_pac, T_pac, sie_pac, P_pac, cv_pac, bmod_pac, dpde_pac, dvdt_pac;
        steelEOS_host.ValuesAtReferenceState(rho, T, sie, P, cv, bmod, dpde, dvdt);
        eospac.ValuesAtReferenceState(rho_pac, T_pac, sie_pac, P_pac, cv_pac, bmod_pac,
                                      dpde_pac, dvdt_pac);
        REQUIRE(isClose(rho, rho_pac));
        REQUIRE(isClose(T, T_pac));
      }

      AND_THEN("We get the correct mean atomic mass and number") {
        REQUIRE(steelEOS_host.MeanAtomicMass() == eospac.MeanAtomicMass());
        REQUIRE(steelEOS_host.MeanAtomicNumber() == eospac.MeanAtomicNumber());
      }

      AND_THEN("We get the correct density and temperature bounds") {
        REQUIRE(std::abs(steelEOS_host.MinimumTemperature() - 1e-2) <= 1e-12);
        REQUIRE(std::abs(steelEOS_host.MinimumDensity() - 1e-8) <= 1e-12);
      }

      // TODO: this needs to be a much more rigorous test
      AND_THEN("Quantities can be read from density and temperature") {
        const Real sie_pac = eospac.InternalEnergyFromDensityTemperature(1e0, 1e6);

        int nw_ie{0};
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::fence();
        using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
        Kokkos::View<int, atomic_view> n_wrong_ie("wrong_ie");
#else
        PortableMDArray<int> n_wrong_ie(&nw_ie, 1);
#endif
        portableFor(
            "calc ie's steel", 0, 100, PORTABLE_LAMBDA(const int &i) {
              const double ie = steelEOS.InternalEnergyFromDensityTemperature(1e0, 1e6);
              if (!isClose(ie, sie_pac)) n_wrong_ie() += 1;
            });
#ifdef PORTABILITY_STRATEGY_KOKKOS
        Kokkos::deep_copy(nw_ie, n_wrong_ie);
#endif
        REQUIRE(nw_ie == 0);
      }

      AND_THEN("rho(P,T) correct for P=1atm, T=freezing") {
        Real T = 273;  // Kelvin
        Real P = 1e6;  // barye
        Real rho, sie; // output
        Real rho_pac, sie_pac;
        std::vector<Real> lambda(steelEOS_host_polymorphic.nlambda());
        steelEOS_host_polymorphic.DensityEnergyFromPressureTemperature(
            P, T, lambda.data(), rho, sie);
        eospac.DensityEnergyFromPressureTemperature(P, T, np<Real>(), rho_pac, sie_pac);
        REQUIRE(isClose(rho, rho_pac));
      }
    }
    // Failing to call finalize leads to a memory leak,
    // but otherwise behaviour is as expected.
    steelEOS_host_polymorphic.Finalize(); // host and device must be
                                          // finalized separately.
    steelEOS.Finalize();                  // cleans up memory on device.
    eospac.Finalize();
  }

  GIVEN("SpinerEOS and EOSPAC for air can be initialized with matid") {
    SpinerEOSDependsRhoT airEOS_host = SpinerEOSDependsRhoT(eosName, airID);
    EOS airEOS = airEOS_host.GetOnDevice();
    EOS eospac = EOSPAC(airID);
    constexpr Real rho = 1e0;
    constexpr Real sie = 2.5e-4;
    THEN("We can get a reference state") {
      Real rho, T, sie, P, cv, bmod, dpde, dvdt;
      Real rho_pac, T_pac, sie_pac, P_pac, cv_pac, bmod_pac, dpde_pac, dvdt_pac;
      airEOS_host.ValuesAtReferenceState(rho, T, sie, P, cv, bmod, dpde, dvdt);
      eospac.ValuesAtReferenceState(rho_pac, T_pac, sie_pac, P_pac, cv_pac, bmod_pac,
                                    dpde_pac, dvdt_pac);
      REQUIRE(isClose(rho, rho_pac));
      REQUIRE(isClose(T, T_pac));
    }
    THEN("Quantities of rho and sie look sane on both host and device") {
      const Real gm1_host = airEOS_host.GruneisenParamFromDensityInternalEnergy(rho, sie);
      const Real T_host = airEOS_host.TemperatureFromDensityInternalEnergy(rho, sie);
      const Real cv_host = airEOS_host.SpecificHeatFromDensityInternalEnergy(rho, sie);

      int nw_gm1{0}, nw_cv{0};
#ifdef PORTABILITY_STRATEGY_KOKKOS
      using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
      Kokkos::fence();
      Kokkos::View<int, atomic_view> n_wrong_gm1("wrong_gm1");
      Kokkos::View<int, atomic_view> n_wrong_cv("wrong_cv");
#else
      PortableMDArray<int> n_wrong_gm1(&nw_gm1, 1);
      PortableMDArray<int> n_wrong_cv(&nw_cv, 1);
#endif
      portableFor(
          "calc gm1 and cv", 0, 100, PORTABLE_LAMBDA(const int &i) {
            const Real gm1 = airEOS.GruneisenParamFromDensityInternalEnergy(rho, sie);
            const Real cv = airEOS.SpecificHeatFromDensityInternalEnergy(rho, sie);
            if (!isClose(gm1, gm1_host)) n_wrong_gm1() += 1;
            if (!isClose(cv, cv_host)) n_wrong_cv() += 1;
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::deep_copy(nw_gm1, n_wrong_gm1);
      Kokkos::deep_copy(nw_cv, n_wrong_cv);
#endif
      REQUIRE(nw_gm1 == 0);
      REQUIRE(nw_cv == 0);
    }
    AND_THEN("P(rho, sie) correct for extrapolation regime") {
      Real rho = 1;
      Real sie = 2.43e16;
      Real P_pac = eospac.PressureFromDensityInternalEnergy(rho, sie);
      Real P_spi = airEOS_host.PressureFromDensityInternalEnergy(rho, sie);
      REQUIRE(isClose(P_pac, P_spi));
    }
    airEOS_host.Finalize();
    airEOS.Finalize();
    eospac.Finalize();
  }

  GIVEN("EOS initialized with matid") {
    SpinerEOSDependsRhoT eos_spiner = SpinerEOSDependsRhoT(eosName, DTID);
    EOS eos_eospac = EOSPAC(DTID);
    Real P = 1e8;
    Real rho = 1.28e-3;
    THEN("Inversion for T(rho,P) works on host") {
      Real T, sie, cv, bmod;
      Real T_pac, sie_pac, cv_pac, bmod_pac;
      const unsigned long output =
          (thermalqs::temperature | thermalqs::specific_internal_energy |
           thermalqs::specific_heat | thermalqs::bulk_modulus);
      eos_spiner.FillEos(rho, T, sie, P, cv, bmod, output);
      eos_eospac.FillEos(rho, T_pac, sie_pac, P, cv_pac, bmod_pac, output);
      REQUIRE(isClose(T, T_pac));
      REQUIRE(isClose(sie, sie_pac));
      REQUIRE(isClose(cv, cv_pac));
    }
    eos_spiner.Finalize();
    eos_eospac.Finalize();
  }

  GIVEN("EOS initialized with matid") {
    SpinerEOSDependsRhoT eos_spiner = SpinerEOSDependsRhoT(eosName, gID);
    EOS eos_eospac = EOSPAC(gID);
    THEN("PT lookup works on the host") {
      Real P = 1e6;          // cgs
      Real T = 0.025 / ev2k; // K
      Real rho, sie;
      Real rho_pac, sie_pac;
      std::vector<Real> lambda(eos_spiner.nlambda());
      eos_spiner.DensityEnergyFromPressureTemperature(P, T, lambda.data(), rho, sie);
      eos_eospac.DensityEnergyFromPressureTemperature(P, T, lambda.data(), rho_pac,
                                                      sie_pac);
      REQUIRE(isClose(rho, rho_pac));
    }
    eos_spiner.Finalize();
    eos_eospac.Finalize();
  }
}

SCENARIO("SpinerEOS depends on rho and sie", "[SpinerEOS][DependsRhoSie]") {

  GIVEN("A titanium EOS in density energy space") {
    using EOS = Variant<SpinerEOSDependsRhoSie, UnitSystem<SpinerEOSDependsRhoSie>>;
    constexpr Real time_unit = 1e6; // These are the numbers you MULTIPLY BY to get cgs
    constexpr Real mass_unit = 1;
    constexpr Real length_unit = 1;
    constexpr Real temp_unit = 1;
    constexpr Real press_unit = (mass_unit * length_unit) / (time_unit * time_unit);
    auto tag = singularity::eos_units_init::LengthTimeUnitsInit();

    EOS eos = SpinerEOSDependsRhoSie(eosName, titaniumID);
    THEN("We can modify it with a unit system with milliseconds") {
      EOS eos_mod =
          eos.Modify<UnitSystem>(tag, time_unit, mass_unit, length_unit, temp_unit);
      Real rho_normal_cgs = 4.447;
      Real T = 298;
      AND_THEN("Pressure from density and temperature produces the right value for both "
               "EOSs") {
        Real P_cgs = eos.PressureFromDensityTemperature(rho_normal_cgs, T);
        Real P_cgmu = eos_mod.PressureFromDensityTemperature(rho_normal_cgs, T);
        REQUIRE(isClose(P_cgs, P_cgmu * press_unit, 1e-10));
      }
      AND_THEN("We can request density energy from pressure temperature") {
        Real P_cgs = 1e6;
        Real P_cgmu = P_cgs / press_unit;
        Real rho_cgs = 4.447;
        Real rho_cgmu = rho_cgs;
        Real sie_cgs, sie_cgmu;
        eos.DensityEnergyFromPressureTemperature(P_cgs, T, static_cast<Real *>(nullptr),
                                                 rho_cgs, sie_cgs);
        eos_mod.DensityEnergyFromPressureTemperature(
            P_cgmu, T, static_cast<Real *>(nullptr), rho_cgmu, sie_cgmu);
        REQUIRE(isClose(rho_cgs, rho_cgmu, 1e-10));
      }
    }
  }

  GIVEN("SpinerEOSes for steel can be initialised with matid") {
    SpinerEOSDependsRhoSie steelEOS_host(eosName, steelID);
    EOS steelEOS = steelEOS_host.GetOnDevice();
    THEN("The correct metadata is read in") {
      REQUIRE(steelEOS_host.matid() == steelID);

      AND_THEN("We get the correct mean atomic mass and number") {
        REQUIRE(isClose(steelEOS_host.MeanAtomicMass(), 55.37));
        REQUIRE(isClose(steelEOS_host.MeanAtomicNumber(), 25.80));
      }
      AND_THEN("We get the correct density and temperature bounds") {
        REQUIRE(std::abs(steelEOS_host.MinimumTemperature() - 1e-2) <= 1e-12);
        REQUIRE(std::abs(steelEOS_host.MinimumDensity() - 1e-8) <= 1e-12);
      }

      int nw_ie2{0}, nw_te2{0};
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::fence();
      using atomic_view = Kokkos::MemoryTraits<Kokkos::Atomic>;
      Kokkos::View<int, atomic_view> n_wrong_ie("wrong_ie");
      Kokkos::View<int, atomic_view> n_wrong_te("wrong_te");
#else
      PortableMDArray<int> n_wrong_ie(&nw_ie2, 1);
      PortableMDArray<int> n_wrong_te(&nw_te2, 1);
#endif
      portableFor(
          "calc ie's steel 2", 0, 100, PORTABLE_LAMBDA(const int &i) {
            const Real ie{steelEOS.InternalEnergyFromDensityTemperature(1e0, 1e6)};
            const Real te{steelEOS.TemperatureFromDensityInternalEnergy(1e0, 1e12)};
            if (!isClose(ie, 4.96415e13)) n_wrong_ie() += 1;
            if (!isClose(te, 75594.6)) n_wrong_te() += 1;
          });
#ifdef PORTABILITY_STRATEGY_KOKKOS
      Kokkos::deep_copy(nw_ie2, n_wrong_ie);
      Kokkos::deep_copy(nw_te2, n_wrong_te);
#endif
      REQUIRE(nw_ie2 == 0);
      REQUIRE(nw_te2 == 0);
    }
    // this can be removed with with reference counting or other tricks
    steelEOS_host.Finalize(); // cleans up host memory
    steelEOS.Finalize();      // cleans up device memory
  }
}

SCENARIO("SpinerEOS and EOSPAC Serialization",
         "[SpinerEOS][DependsRhoT][DependsRhoSie][EOSPAC][Serialization]") {
  GIVEN("Eoses initialized with matid") {
    SpinerEOSDependsRhoT rhoT_orig = SpinerEOSDependsRhoT(eosName, steelID);
    SpinerEOSDependsRhoSie rhoSie_orig = SpinerEOSDependsRhoSie(eosName, steelID);
    EOS eospac_orig = EOSPAC(steelID);
    // we want to stress test that we can serialize and deserialize
    // multiple EOSPAC objects
    EOS eospac_air = EOSPAC(airID);
    THEN("They report dynamic vs static memory correctly") {
      REQUIRE(rhoT_orig.AllDynamicMemoryIsShareable());
      REQUIRE(rhoSie_orig.AllDynamicMemoryIsShareable());
      REQUIRE(!eospac_orig.AllDynamicMemoryIsShareable());
      REQUIRE(eospac_orig.SerializedSizeInBytes() >
              eospac_orig.DynamicMemorySizeInBytes());
      REQUIRE(eospac_air.SerializedSizeInBytes() > eospac_air.DynamicMemorySizeInBytes());
    }
    WHEN("We serialize") {
      auto [rhoT_size, rhoT_data] = rhoT_orig.Serialize();
      REQUIRE(rhoT_size == rhoT_orig.SerializedSizeInBytes());

      auto [rhoSie_size, rhoSie_data] = rhoSie_orig.Serialize();
      REQUIRE(rhoSie_size == rhoSie_orig.SerializedSizeInBytes());

      auto [eospac_size, eospac_data] = eospac_orig.Serialize();
      REQUIRE(eospac_size == eospac_orig.SerializedSizeInBytes());

      auto [air_size, air_data] = eospac_air.Serialize();
      REQUIRE(air_size == eospac_air.SerializedSizeInBytes());

      const std::size_t rhoT_shared_size = rhoT_orig.DynamicMemorySizeInBytes();
      REQUIRE(rhoT_size > rhoT_shared_size);

      const std::size_t rhoSie_shared_size = rhoSie_orig.DynamicMemorySizeInBytes();
      REQUIRE(rhoSie_size > rhoSie_shared_size);

      const std::size_t eospac_shared_size = eospac_orig.DynamicMemorySizeInBytes();
      REQUIRE(eospac_size > eospac_shared_size);

      const std::size_t air_shared_size = eospac_air.DynamicMemorySizeInBytes();
      REQUIRE(air_size > air_shared_size);

      THEN("We can deserialize into shared memory") {
        using singularity::SharedMemSettings;
        using RhoTTricks = singularity::table_utils::SpinerTricks<SpinerEOSDependsRhoT>;
        using RhoSieTricks =
            singularity::table_utils::SpinerTricks<SpinerEOSDependsRhoSie>;

        char *rhoT_shared_data = (char *)malloc(rhoT_shared_size);
        char *rhoSie_shared_data = (char *)malloc(rhoSie_shared_size);
        char *eospac_shared_data = (char *)malloc(eospac_shared_size);
        char *air_shared_data = (char *)malloc(air_shared_size);

        SpinerEOSDependsRhoT eos_rhoT;
        std::size_t read_size_rhoT =
            eos_rhoT.DeSerialize(rhoT_data, SharedMemSettings(rhoT_shared_data, true));
        REQUIRE(read_size_rhoT == rhoT_size);
        REQUIRE(RhoTTricks::DataBoxesPointToDifferentMemory(rhoT_orig, eos_rhoT));

        SpinerEOSDependsRhoSie eos_rhoSie;
        std::size_t read_size_rhoSie = eos_rhoSie.DeSerialize(
            rhoSie_data, SharedMemSettings(rhoSie_shared_data, true));
        REQUIRE(read_size_rhoSie == rhoSie_size);
        REQUIRE(RhoSieTricks::DataBoxesPointToDifferentMemory(rhoSie_orig, eos_rhoSie));

        eospac_orig.Finalize();
        EOS eos_eospac = EOSPAC();
        std::size_t read_size_eospac = eos_eospac.DeSerialize(
            eospac_data, SharedMemSettings(eospac_shared_data, true));
        REQUIRE(read_size_eospac == eospac_size);

        eospac_air.Finalize();
        EOS eos_air_2 = EOSPAC();
        std::size_t read_size_air =
            eos_air_2.DeSerialize(air_data, SharedMemSettings(air_shared_data, true));
        REQUIRE(read_size_air == air_size);

        AND_THEN("EOS lookups work") {
          constexpr Real rho_trial = 1;
          constexpr Real sie_trial = 1e12;
          const Real P_eospac =
              eos_eospac.PressureFromDensityInternalEnergy(rho_trial, sie_trial);
          const Real P_spiner_orig =
              rhoT_orig.PressureFromDensityInternalEnergy(rho_trial, sie_trial);
          const Real P_spiner_rhoT =
              eos_rhoT.PressureFromDensityInternalEnergy(rho_trial, sie_trial);
          const Real P_spiner_rhoSie =
              eos_rhoSie.PressureFromDensityInternalEnergy(rho_trial, sie_trial);
          REQUIRE(isClose(P_eospac, P_spiner_orig));
          REQUIRE(isClose(P_eospac, P_spiner_rhoT));
          REQUIRE(isClose(P_eospac, P_spiner_rhoSie));

          // Check that MinInternalEnergyFromDensity serializes correctly for all EOS
          const Real ecold_eospac = eos_eospac.MinInternalEnergyFromDensity(rho_trial);
          const Real ecold_spiner_orig =
              rhoT_orig.MinInternalEnergyFromDensity(rho_trial);
          const Real ecold_spiner_rhoT = eos_rhoT.MinInternalEnergyFromDensity(rho_trial);
          const Real ecold_spiner_rhoSie =
              eos_rhoSie.MinInternalEnergyFromDensity(rho_trial);
          REQUIRE(isClose(ecold_eospac, ecold_spiner_orig));
          REQUIRE(isClose(ecold_eospac, ecold_spiner_rhoT));
          REQUIRE(isClose(ecold_eospac, ecold_spiner_rhoSie));
        }

        eos_eospac.Finalize();
        free(rhoT_shared_data);
        free(rhoSie_shared_data);
        free(eospac_shared_data);
      }
      free(rhoT_data);
      free(rhoSie_data);
      free(eospac_data);
    }

    rhoT_orig.Finalize();
    rhoSie_orig.Finalize();
  }
}

SCENARIO("Tabulated EOS with partial ionization",
         "[SpinerEOS][DependsRhoT][DependsRhoSie][EOSPAC][SplitTables]") {
  using singularity::TableSplit;
  GIVEN("Spiner EOSes for total, electron, and ion EOS") {
    SpinerEOSDependsRhoT rt_tot(eosName, steelID, TableSplit::Total);
    SpinerEOSDependsRhoT rt_e(eosName, steelID, TableSplit::ElectronOnly);
    SpinerEOSDependsRhoT rt_ic(eosName, steelID, TableSplit::IonCold);
    std::vector<SpinerEOSDependsRhoT> rt = {rt_tot, rt_e, rt_ic};

    SpinerEOSDependsRhoSie re_tot(eosName, steelID, TableSplit::Total);
    SpinerEOSDependsRhoSie re_e(eosName, steelID, TableSplit::ElectronOnly);
    SpinerEOSDependsRhoSie re_ic(eosName, steelID, TableSplit::IonCold);
    std::vector<SpinerEOSDependsRhoSie> re = {re_tot, re_e, re_ic};

    EOSPAC pac_tot(steelID, TableSplit::Total);
    EOSPAC pac_e(steelID, TableSplit::ElectronOnly);
    EOSPAC pac_ic(steelID, TableSplit::IonCold);
    std::vector<EOSPAC> pac = {pac_tot, pac_e, pac_ic};

    WHEN("We evaluate pressure from density and temperature") {
      constexpr Real rho = 1e1;
      constexpr Real T = 1e4;
      const std::size_t N = rt.size();
      std::vector<Real> P_rt(N), P_re(N), P_pac(N);
      for (std::size_t i = 0; i < N; ++i) {
        P_rt[i] = rt[i].PressureFromDensityTemperature(rho, T);
        P_re[i] = re[i].PressureFromDensityTemperature(rho, T);
        P_pac[i] = pac[i].PressureFromDensityTemperature(rho, T);
      }
      THEN("The subtables add to the total for DependsRhoT") {
        REQUIRE(isClose(P_rt[0], P_rt[1] + P_rt[2]));
      }
      THEN("The subtables add to the total for DependsRhoSie") {
        REQUIRE(isClose(P_re[0], P_re[1] + P_re[2]));
      }
      THEN("The subtables add to the total for EOSPAC") {
        REQUIRE(isClose(P_pac[0], P_pac[1] + P_pac[2]));
      }
      THEN("DependsRhoT agrees with EOSPAC") {
        for (std::size_t i = 0; i < N; ++i) {
          REQUIRE(isClose(P_rt[i], P_pac[i]));
        }
      }
      THEN("DependsRhoSie agrees with EOSPAC") {
        for (std::size_t i = 0; i < N; ++i) {
          REQUIRE(isClose(P_re[i], P_pac[i]));
        }
      }
    }

    for (auto &eos : rt) {
      eos.Finalize();
    }
    for (auto &eos : re) {
      eos.Finalize();
    }
    for (auto &eos : pac) {
      eos.Finalize();
    }
  }
}

#endif // SINGULARITY_USE_EOSPAC
#endif // SINGULARITY_TEST_SESAME
#endif // SPINER_USE_HDF
