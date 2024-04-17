#------------------------------------------------------------------------------
# © 2021-2023. Triad National Security, LLC. All rights reserved.  This
# program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by Triad
# National Security, LLC for the U.S.  Department of Energy/National
# Nuclear Security Administration. All rights in the program are
# reserved by Triad National Security, LLC, and the U.S. Department of
# Energy/National Nuclear Security Administration. The Government is
# granted for itself and others acting on its behalf a nonexclusive,
# paid-up, irrevocable worldwide license in this material to reproduce,
# prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.
#------------------------------------------------------------------------------
import math
import unittest
import singularity_eos
import numpy as np
from numpy.testing import assert_allclose

class EOSTestBase(object):
    def assertIsClose(self, a, b, eps=5e-2):
        rel = abs(b - a) / (abs(a + b) + 1e-20)
        if rel > eps:
            raise AssertionError("{} > {}".format(rel, eps))

    def assertEqualEOS(self, test_e, ref_e):
        """compare all individual member functions with 1 as inputs,
        this function is meant to catch mis-implementations of
        modifiers that can be initialized in such a way as to
        be equivalent of an unmodified eos. Best used with analytic eoss.
        """
        self.assertIsClose(test_e.TemperatureFromDensityInternalEnergy(1, 1),
            ref_e.TemperatureFromDensityInternalEnergy(1, 1), 1.e-15)
        self.assertIsClose(test_e.InternalEnergyFromDensityTemperature(1, 1),
            ref_e.InternalEnergyFromDensityTemperature(1, 1), 1.e-15)
        self.assertIsClose(test_e.PressureFromDensityInternalEnergy(1, 1),
            ref_e.PressureFromDensityInternalEnergy(1, 1), 1.e-15)
        self.assertIsClose(test_e.SpecificHeatFromDensityInternalEnergy(1, 1),
            ref_e.SpecificHeatFromDensityInternalEnergy(1, 1), 1.e-15)
        self.assertIsClose(test_e.BulkModulusFromDensityInternalEnergy(1, 1),
            ref_e.BulkModulusFromDensityInternalEnergy(1, 1), 1.e-15)
        self.assertIsClose(test_e.GruneisenParamFromDensityInternalEnergy(1, 1),
            ref_e.GruneisenParamFromDensityInternalEnergy(1, 1), 1.e-15)
        self.assertIsClose(test_e.PressureFromDensityTemperature(1, 1),
            ref_e.PressureFromDensityTemperature(1, 1), 1.e-15)
        self.assertIsClose(test_e.SpecificHeatFromDensityTemperature(1, 1),
            ref_e.SpecificHeatFromDensityTemperature(1, 1), 1.e-15)
        self.assertIsClose(test_e.BulkModulusFromDensityTemperature(1, 1),
            ref_e.BulkModulusFromDensityTemperature(1, 1), 1.e-15)
        self.assertIsClose(test_e.GruneisenParamFromDensityTemperature(1, 1),
            ref_e.GruneisenParamFromDensityTemperature(1, 1), 1.e-15)
        self.assertIsClose(test_e.MinimumDensity(), ref_e.MinimumDensity(), 1.e-15)
        self.assertIsClose(test_e.MinimumTemperature(), ref_e.MinimumTemperature(), 1.e-15)

class EOS(unittest.TestCase):
    def testConstants(self):
        from singularity_eos import thermalqs
        self.assertEqual(thermalqs.all_values, thermalqs.none |
                                               thermalqs.density |
                                               thermalqs.specific_internal_energy |
                                               thermalqs.pressure      |
                                               thermalqs.temperature   |
                                               thermalqs.specific_heat |
                                               thermalqs.bulk_modulus  |
                                               thermalqs.do_lambda)

    def testIdealGas(self):
        eos = singularity_eos.IdealGas(1,1)

    def testShiftedIdealGas(self):
        eos = singularity_eos.Shifted(singularity_eos.IdealGas(1,1),1)

    def testScaledIdealGas(self):
        eos = singularity_eos.Scaled(singularity_eos.IdealGas(1,1),1)

    def testScaledShiftedIdealGas(self):
        shifted = singularity_eos.Shifted(singularity_eos.IdealGas(1,1),1)
        eos = singularity_eos.Scaled(shifted,1)

    def testGruneisen(self):
        eos = singularity_eos.Gruneisen(1,1,1,1,1,1,1,1,1,1)

    def testJWL(self):
        eos = singularity_eos.JWL(1,1,1,1,1,1,1)

    def testDavisReactants(self):
        eos = singularity_eos.DavisReactants(1,1,1,1,1,1,1,1,1,1,1)

    def testDavisProducts(self):
        eos = singularity_eos.DavisProducts(1,1,1,1,1,1,1,1)

class Modifiers(unittest.TestCase, EOSTestBase):
    def setUp(self):
        "Parameters for a shifted and scaled ideal gas"
        self.Cv = 2.0
        self.gm1 = 0.5
        self.scale = 2.0
        self.shift = 0.1
        self.rho = 2.0
        self.sie = 0.5

    def testScaledShiftedIdealGas(self):
        from singularity_eos import IdealGas, Shifted, Scaled

        # We construct a shifted, scaled IdealGas by hand
        a = IdealGas(self.gm1, self.Cv)
        b = Shifted(a, self.shift)
        eos = Scaled(b, self.scale)

        # The shift and scale parameters pass through correctly"
        self.assertIsClose(eos.PressureFromDensityInternalEnergy(self.rho, self.sie), 0.3)

    def testBilinearRampScaledShiftedIdealGas(self):
        from singularity_eos import IdealGas, Shifted, Scaled, BilinearRamp

        # We construct a shifted, scaled IdealGas by hand
        a = IdealGas(self.gm1, self.Cv)
        b = Shifted(a, self.shift)
        eos = Scaled(b, self.scale)

        # We add a ramp
        r0 = 1
        a = 1
        b = 0
        c = 0
        eos_ramped = BilinearRamp(eos, r0, a, b, c)

    def testNonModifying(self):
        from singularity_eos import IdealGas, Shifted, Scaled
        ig = IdealGas(self.gm1, self.Cv)
        igsh = Scaled(IdealGas(self.gm1, self.Cv), 1.0)
        igsc = Shifted(IdealGas(self.gm1, self.Cv), 0.0)
        #int enabled[4] = {0, 0, 1, 0}
        #Real vals[6] = {0.0, 0.0, 1.e9, 1.0, 2.0, 1.0}
        #Real rho0 = 1.e6 / (gm1 * Cv * 293.0)
        #init_sg_IdealGas(0, &igra, gm1, Cv, enabled, vals)

        # The modified EOS should produce equivalent results
        self.assertEqualEOS(igsh, ig)
        self.assertEqualEOS(igsc, ig)
        #  compare_two_eoss(igra, ig)

    def testRelativisticIdealGas(self):
        """[RelativisticEOS][IdealGas]"""
        from singularity_eos import IdealGas, Relativistic
        Cv = 2.0
        gm1 = 0.5
        eos = Relativistic(IdealGas(gm1, Cv), 1.0)

        # The EOS has finite sound speeds"
        rho = 1e3
        sie = 1e3
        bmod = eos.BulkModulusFromDensityInternalEnergy(rho, sie)
        cs2 = bmod / rho
        self.assertLess(cs2, 1)

    def testUnitSystemIdealGas_ThermalUnits(self):
        """[UnitSystem][IdealGas][ThermalUnits]"""
        from singularity_eos import IdealGas, UnitSystem

        # Parameters for an ideal gas
        Cv = 2.0
        gm1 = 0.5

        # Units with a thermal unit system
        rho_unit = 1e1
        sie_unit = 1e-1
        temp_unit = 123
        eos = UnitSystem(IdealGas(gm1, Cv), rho_unit=rho_unit, sie_unit=sie_unit, temp_unit=temp_unit)

        # Units cancel out for an ideal gas
        rho = 1e3
        sie = 1e3
        P = eos.PressureFromDensityInternalEnergy(rho, sie)
        Ptrue = gm1 * rho * sie
        self.assertLess(abs(P - Ptrue) / Ptrue,1e-3)

    def testUnitSystemIdealGas_LengthTimeUnits(self):
        """[UnitSystem][IdealGas][LengthTimeUnits]"""
        from singularity_eos import IdealGas, UnitSystem
        from singularity_eos.eos_units import LengthTimeUnits

        # Parameters for an ideal gas
        Cv = 2.0
        gm1 = 0.5

        # Units with length and time units
        time_unit = 456
        length_unit = 1e2
        mass_unit = 1e6
        temp_unit = 789
        eos = UnitSystem(IdealGas(gm1, Cv), units=LengthTimeUnits, time_unit=time_unit, length_unit=length_unit, mass_unit=mass_unit, temp_unit=temp_unit)

        # Units cancel out for an ideal gas
        rho = 1e3
        sie = 1e3
        P = eos.PressureFromDensityInternalEnergy(rho, sie)
        Ptrue = gm1 * rho * sie
        self.assertLess(abs(P - Ptrue) / Ptrue, 1e-3)

    def testBilinearRamp(self):
        from singularity_eos import IdealGas, BilinearRamp, pAlpha2BilinearRampParams
        Cv = 2.0;
        gm1 = 0.5;
        # We construct a ramp from a p-alpha model
        Pe = 5.e7
        Pc = 1.e8
        alpha0 = 1.5
        T0 = 293.0
        rho0 = 1.e6 / (gm1 * Cv * T0)
        r0 = rho0 / alpha0
        r1 = Pc / (gm1 * Cv * T0)
        rmid = Pe / (gm1 * Cv * T0 * alpha0)
        # P(alpha0 * rmid)
        P_armid = alpha0 * gm1 * Cv * rmid * T0

        ig = IdealGas(gm1, Cv)
        param_r0, param_a, param_b, param_c = pAlpha2BilinearRampParams(ig, alpha0, Pe, Pc)
        igra = BilinearRamp(IdealGas(gm1, Cv), param_r0, param_a, param_b, param_c)

        # construct ramp params and evaluate directly for test
        a = r0 * Pe / (rmid - r0)
        b = r0 * (Pc - Pe) / (r1 - rmid)
        c = (Pc * rmid - Pe * r1) / (r0 * (Pc - Pe))
        # density in the middle of the first slope
        rho_t1 = 0.5 * (r0 + rmid)
        # density in the middle of the second slope
        rho_t2 = 0.5 * (rmid + r1)
        # P (rho_t1) note that r0 = rho0 / alpha0
        Prhot1 = a * (rho_t1 / r0 - 1.0)
        # P (rho_t2)
        Prhot2 = b * (rho_t2 / r0 - c)
        # bmod (rho_t1)
        bmodrt1 = rho_t1 * a / r0
        # bmod (rho_t2)
        bmodrt2 = rho_t2 * b / r0

        # Then P_eos(alpha_0*rmid, T0) = P_ramp(rmid,T0)
        self.assertIsClose(P_armid, igra.PressureFromDensityTemperature(rmid, T0), 1.e-12)

        # We obtain correct ramp behavior in P(rho) for rho <r0, [r0,rmid], [rmid,r1] and >r1
        self.assertIsClose(Prhot1, igra.PressureFromDensityTemperature(rho_t1, T0), 1.e-12)
        self.assertIsClose(Prhot2, igra.PressureFromDensityTemperature(rho_t2, T0), 1.e-12)

        # check pressure below and beyond ramp matches unmodified ideal gas
        self.assertIsClose(0.8 * r0 * gm1 * Cv * T0, igra.PressureFromDensityTemperature(0.8 * r0, T0), 1.e-12)
        self.assertIsClose(1.2 * r1 * gm1 * Cv * T0, igra.PressureFromDensityTemperature(1.2 * r1, T0), 1.e-12)

        # We obtain correct ramp behavior in bmod(rho) for rho <r0, [r0,rmid], [rmid,r1] and >r1
        # check bulk moduli on both pieces of ramp
        self.assertIsClose(bmodrt1, igra.BulkModulusFromDensityTemperature(rho_t1, T0), 1.e-12)
        self.assertIsClose(bmodrt2, igra.BulkModulusFromDensityTemperature(rho_t2, T0), 1.e-12)

        # check bulk modulus below and beyond ramp matches unmodified ideal gas
        self.assertIsClose(0.8 * r0 * gm1 * (gm1 + 1.0) * Cv * T0, igra.BulkModulusFromDensityTemperature(0.8 * r0, T0), 1.e-12)
        self.assertIsClose(1.2 * r1 * gm1 * (gm1 + 1.0) * Cv * T0, igra.BulkModulusFromDensityTemperature(1.2 * r1, T0), 1.e-12)

class VectorEOS_IdealGas_Given_Rho_Sie(unittest.TestCase):
    def setUp(self):
        Cv = 5.0
        gm1 = 0.4
        self.eos = singularity_eos.IdealGas(gm1, Cv)

        self.num = 3
        self.density = np.zeros(self.num)
        self.energy = np.zeros(self.num)

        # Populate the input arrays
        self.density[0] = 1.0
        self.density[1] = 2.0
        self.density[2] = 5.0
        self.energy[0] = 5.0
        self.energy[1] = 10.0
        self.energy[2] = 15.0

        # Gold standard values
        self.pressure_true = np.array((2.0, 8.0, 30.0))
        self.temperature_true = np.array((1., 2., 3.))
        self.bulkmodulus_true = np.array((2.8, 11.2, 42.))
        self.heatcapacity_true = np.array((Cv, Cv, Cv))
        self.gruneisen_true = np.array((gm1, gm1, gm1))

        # Create arrays for the outputs
        self.temperature = np.zeros(self.num)
        self.pressure = np.zeros(self.num)
        self.heatcapacity = np.zeros(self.num)
        self.bulkmodulus = np.zeros(self.num)
        self.gruneisen = np.zeros(self.num)

    def test_temp(self):
        """[Vector EOS][IdealGas][Energies and densities] A T(rho, e) lookup is performed"""
        self.eos.TemperatureFromDensityInternalEnergy(self.density, self.energy, self.temperature, self.num)
        assert_allclose(self.temperature, self.temperature_true, rtol=1e-12)

    def test_pressure(self):
        """[Vector EOS][IdealGas][Energies and densities] A P(rho, e) lookup is performed"""
        self.eos.PressureFromDensityInternalEnergy(self.density, self.energy, self.pressure, self.num)
        assert_allclose(self.pressure, self.pressure_true, rtol=1e-12)

    def test_cv(self):
        """[Vector EOS][IdealGas][Energies and densities] A C_v(rho, e) lookup is performed"""
        self.eos.SpecificHeatFromDensityInternalEnergy(self.density, self.energy, self.heatcapacity, self.num)
        assert_allclose(self.heatcapacity, self.heatcapacity_true, rtol=1e-12)

    def test_bmod(self):
        """[Vector EOS][IdealGas][Energies and densities] A B_S(rho, e) lookup is performed"""
        self.eos.BulkModulusFromDensityInternalEnergy(self.density, self.energy, self.bulkmodulus, self.num)
        assert_allclose(self.bulkmodulus, self.bulkmodulus_true, rtol=1e-12)

    def test_gamma(self):
        """[Vector EOS][IdealGas][Energies and densities] A Gamma(rho, e) lookup is performed"""
        self.eos.GruneisenParamFromDensityInternalEnergy(self.density, self.energy, self.gruneisen, self.num)
        assert_allclose(self.gruneisen, self.gruneisen_true, rtol=1e-12)


class VectorEOS_IdealGas_Given_Rho_Temp(unittest.TestCase):
    def setUp(self):
        Cv = 5.0
        gm1 = 0.4
        self.eos = singularity_eos.IdealGas(gm1, Cv)

        self.num = 3
        self.density = np.zeros(self.num)
        self.temperature = np.zeros(self.num)

        # Populate the input arrays
        self.density[0] = 1.0
        self.density[1] = 2.0
        self.density[2] = 5.0
        self.temperature[0] = 50.0
        self.temperature[1] = 100.0
        self.temperature[2] = 150.0

        # Gold standard values
        self.energy_true = np.array((250., 500., 750.))
        self.pressure_true = np.array((100., 400., 1500.))
        self.bulkmodulus_true = np.array((140., 560., 2100.))
        self.heatcapacity_true = np.array((Cv, Cv, Cv))
        self.gruneisen_true = np.array((gm1, gm1, gm1))

        # Create arrays for the outputs
        self.energy = np.zeros(self.num)
        self.pressure = np.zeros(self.num)
        self.heatcapacity = np.zeros(self.num)
        self.bulkmodulus = np.zeros(self.num)
        self.gruneisen = np.zeros(self.num)

    def test_energy(self):
        """[Vector EOS][IdealGas][Densities and temperatures] A e(rho, T) lookup is performed"""
        self.eos.InternalEnergyFromDensityTemperature(self.density, self.temperature, self.energy, self.num)
        assert_allclose(self.energy, self.energy_true, rtol=1e-12)

    def test_pressure(self):
        """[Vector EOS][IdealGas][Densities and temperatures] A P(rho, T) lookup is performed"""
        self.eos.PressureFromDensityTemperature(self.density, self.temperature, self.pressure, self.num)
        assert_allclose(self.pressure, self.pressure_true, rtol=1e-12)

    def test_cv(self):
        """[Vector EOS][IdealGas][Densities and temperatures] A C_v(rho, T) lookup is performed"""
        self.eos.SpecificHeatFromDensityTemperature(self.density, self.temperature, self.heatcapacity, self.num)
        assert_allclose(self.heatcapacity, self.heatcapacity_true, rtol=1e-12)

    def test_bmod(self):
        """[Vector EOS][IdealGas][Densities and temperatures] A B_S(rho, T) lookup is performed"""
        self.eos.BulkModulusFromDensityTemperature(self.density, self.temperature, self.bulkmodulus, self.num)
        assert_allclose(self.bulkmodulus, self.bulkmodulus_true, rtol=1e-12)

    def test_gamma(self):
        """[Vector EOS][IdealGas][Densities and temperatures] A Gamma(rho, T) lookup is performed"""
        self.eos.GruneisenParamFromDensityTemperature(self.density, self.temperature, self.gruneisen, self.num)
        assert_allclose(self.gruneisen, self.gruneisen_true, rtol=1e-12)


@unittest.skipIf('SpinerEOSDependsRhoT' not in dir(singularity_eos) or 'EOSPAC' not in dir(singularity_eos), "No Spiner or EOSPAC support")
class SpinerEOSdependsOnRhoT_Steel(unittest.TestCase, EOSTestBase):
    "[SpinerEOS],[DependsRhoT][EOSPAC]"

    def setUp(self):
        "SpinerEOS and EOSPAC EOS for steel can be initialized with matid"
        self.eosName = "../materials.sp5"
        self.steelID = 4272
        self.steelName = "stainless steel 347"

        self.steelEOS_host = singularity_eos.SpinerEOSDependsRhoT(self.eosName, self.steelID)
        self.eospac = singularity_eos.EOSPAC(self.steelID)

    def test_metadata(self):
        "The correct metadata is read in"
        self.assertEqual(self.steelEOS_host.matid, self.steelID)

    def test_reference(self):
        "We can get a reference density and temperature"
        state = self.steelEOS_host.ValuesAtReferenceState()
        state_pac = self.eospac.ValuesAtReferenceState()
        self.assertIsClose(state.density, state_pac.density)
        self.assertIsClose(state.temperature, state_pac.temperature)

    # TODO: this needs to be a much more rigorous test
    def test_quantities_from_rho_temp(self):
        "Quantities can be read from density and temperature"
        sie_pac = self.eospac.InternalEnergyFromDensityTemperature(1e0, 1e6)
        ie = self.steelEOS_host.InternalEnergyFromDensityTemperature(1e0, 1e6)
        self.assertIsClose(ie, sie_pac)

    def test_rho_of_P_T(self):
        "rho(P,T) correct for P=1atm, T=freezing"
        T = 273  # Kelvin
        P = 1e6  # barye
        lmbda = np.zeros(self.steelEOS_host.nlambda, dtype=np.double)
        rho, sie = self.steelEOS_host.DensityEnergyFromPressureTemperature(P, T, lmbda)
        rho_pac, sie_pac = self.eospac.DensityEnergyFromPressureTemperature(P, T)
        self.assertIsClose(rho, rho_pac)

    def tearDown(self):
        # Failing to call finalize leads to a memory leak,
        # but otherwise behaviour is as expected.
        # It's possible to this automatically clean up with
        # some form of reference counting. If this is a priority,
        # we can re-examine.
        self.steelEOS_host.Finalize()
        self.eospac.Finalize()


@unittest.skipIf('SpinerEOSDependsRhoT' not in dir(singularity_eos) or 'EOSPAC' not in dir(singularity_eos), "No Spiner or EOSPAC support")
class SpinerEOSdependsOnRhoT_Air(unittest.TestCase, EOSTestBase):
    "SpinerEOS and EOSPAC for air can be initialized with matid"

    def setUp(self):
        self.eosName = "../materials.sp5"
        self.airID = 5030
        self.airEOS_host = singularity_eos.SpinerEOSDependsRhoT(self.eosName, self.airID)
        self.eospac = singularity_eos.EOSPAC(self.airID)

    def test_reference(self):
        "We can get a reference state"
        state = self.airEOS_host.ValuesAtReferenceState()
        state_pac = self.eospac.ValuesAtReferenceState()
        self.assertIsClose(state.density, state_pac.density)
        self.assertIsClose(state.temperature, state_pac.temperature)

    def test_P_from_rho_sie(self):
        "P(rho, sie) correct for extrapolation regime"
        rho = 1
        sie = 2.43e16
        P_pac = self.eospac.PressureFromDensityInternalEnergy(rho, sie)
        P_spi = self.airEOS_host.PressureFromDensityInternalEnergy(rho, sie)
        self.assertIsClose(P_pac, P_spi)

    def tearDown(self):
        self.airEOS_host.Finalize()


@unittest.skipIf('SpinerEOSDependsRhoT' not in dir(singularity_eos) or 'EOSPAC' not in dir(singularity_eos), "No Spiner or EOSPAC support")
class EOS_init_with_matid2(unittest.TestCase, EOSTestBase):
    "EOS initialized with matid"

    def setUp(self):
        self.eosName = "../materials.sp5"
        self.DTID = 5267
        self.eos_spiner = singularity_eos.SpinerEOSDependsRhoT(self.eosName, self.DTID)
        self.eos_eospac = singularity_eos.EOSPAC(self.DTID)

    def test_inversion_for_T_rho_P(self):
        "Inversion for T(rho,P) works on host"
        from singularity_eos import thermalqs
        P = 1e8
        rho = 1.28e-3
        output = (thermalqs.temperature | thermalqs.specific_internal_energy | thermalqs.specific_heat | thermalqs.bulk_modulus)
        state = self.eos_spiner.FillEos(rho=rho, press=P, output=output)
        state_pac = self.eos_eospac.FillEos(rho=rho, press=P, output=output)
        self.assertIsClose(state.temperature, state_pac.temperature)
        self.assertIsClose(state.specific_internal_energy, state_pac.specific_internal_energy)
        self.assertIsClose(state.specific_heat, state_pac.specific_heat)

    def tearDown(self):
        self.eos_spiner.Finalize()
        self.eos_eospac.Finalize()


@unittest.skipIf('SpinerEOSDependsRhoT' not in dir(singularity_eos) or 'EOSPAC' not in dir(singularity_eos), "No Spiner or EOSPAC support")
class EOS_init_with_matid(unittest.TestCase, EOSTestBase):
    "EOS initialized with matid"

    def setUp(self):
        self.eosName = "../materials.sp5"
        self.gID = 2700
        self.eos_spiner = singularity_eos.SpinerEOSDependsRhoT(self.eosName, self.gID)
        self.eos_eospac = singularity_eos.EOSPAC(self.gID)

    def test_P_T_lookup(self):
        "PT lookup works on the host"
        ev2k = 1.160451812e4
        P = 1e6          # cgs
        T = 0.025 / ev2k # K
        lmbda = np.zeros(self.eos_spiner.nlambda, dtype=np.double)
        rho, sie = self.eos_spiner.DensityEnergyFromPressureTemperature(P, T, lmbda)
        rho_pac, sie_pac = self.eos_eospac.DensityEnergyFromPressureTemperature(P, T, lmbda)
        self.assertIsClose(rho, rho_pac)

    def tearDown(self):
        self.eos_spiner.Finalize()

@unittest.skipIf('SpinerEOSDependsRhoT' not in dir(singularity_eos) or 'EOSPAC' not in dir(singularity_eos), "No Spiner or EOSPAC support")
class VectorEOS_EOSPAC_Steel(unittest.TestCase, EOSTestBase):
    "[Vector EOS][EOSPAC]"

    def setUp(self):
        "SpinerEOS and EOSPAC EOS for steel can be initialized with matid"
        self.eosName = "../materials.sp5"
        self.steelID = 4272
        self.steelName = "stainless steel 347"

        self.steelEOS_host = singularity_eos.SpinerEOSDependsRhoT(self.eosName, self.steelID)
        self.eospac = singularity_eos.EOSPAC(self.steelID)

        self.num = 3

        self.density = np.zeros(self.num)
        self.temperature = np.zeros(self.num)

        # Populate the input arrays
        self.density[0] = 1e0
        self.density[1] = 1e0
        self.density[2] = 1e0
        self.temperature[0] = 1e6
        self.temperature[1] = 1e6
        self.temperature[2] = 1e6

    def test_energy(self):
        """[Vector EOS][Steel][Densities and temperatures] A e(rho, T) lookup is performed"""
        energy_eospac = np.zeros(self.num)
        energy_spiner = np.zeros(self.num)
        num_scratch = int(self.eospac.scratch_size("InternalEnergyFromDensityTemperature", self.num) / np.double(1).nbytes)
        scratch = np.zeros(num_scratch)
        self.eospac.InternalEnergyFromDensityTemperature(self.density, self.temperature, energy_eospac, scratch, self.num)
        self.steelEOS_host.InternalEnergyFromDensityTemperature(self.density, self.temperature, energy_spiner, self.num)

        for a, b in zip(energy_eospac, energy_eospac):
            self.assertIsClose(a, b)

    def test_pressure(self):
        """[Vector EOS][Steel][Densities and temperatures] A P(rho, T) lookup is performed"""
        pressure_eospac = np.zeros(self.num)
        pressure_spiner = np.zeros(self.num)
        num_scratch = int(self.eospac.scratch_size("PressureFromDensityTemperature", self.num) / np.double(1).nbytes)
        scratch = np.zeros(num_scratch)
        self.eospac.PressureFromDensityTemperature(self.density, self.temperature, pressure_eospac, scratch, self.num)
        self.steelEOS_host.PressureFromDensityTemperature(self.density, self.temperature, pressure_spiner, self.num)

        for a, b in zip(pressure_eospac, pressure_eospac):
            self.assertIsClose(a, b)

    def test_cv(self):
        """[Vector EOS][Steel][Densities and temperatures] A C_v(rho, T) lookup is performed"""
        heatcapacity_eospac = np.zeros(self.num)
        heatcapacity_spiner = np.zeros(self.num)
        num_scratch = int(self.eospac.scratch_size("SpecificHeatFromDensityTemperature", self.num) / np.double(1).nbytes)
        scratch = np.zeros(num_scratch)
        self.eospac.SpecificHeatFromDensityTemperature(self.density, self.temperature, heatcapacity_eospac, scratch, self.num)
        self.steelEOS_host.SpecificHeatFromDensityTemperature(self.density, self.temperature, heatcapacity_spiner, self.num)

        for a, b in zip(heatcapacity_eospac, heatcapacity_eospac):
            self.assertIsClose(a, b)

    def test_bmod(self):
        """[Vector EOS][Steel][Densities and temperatures] A B_S(rho, T) lookup is performed"""
        bulkmodulus_eospac = np.zeros(self.num)
        bulkmodulus_spiner = np.zeros(self.num)
        num_scratch = int(self.eospac.scratch_size("BulkModulusFromDensityTemperature", self.num) / np.double(1).nbytes)
        scratch = np.zeros(num_scratch)
        self.eospac.BulkModulusFromDensityTemperature(self.density, self.temperature, bulkmodulus_eospac, scratch, self.num)
        self.steelEOS_host.BulkModulusFromDensityTemperature(self.density, self.temperature, bulkmodulus_spiner, self.num)

        for a, b in zip(bulkmodulus_eospac, bulkmodulus_eospac):
            self.assertIsClose(a, b)

    def test_gamma(self):
        """[Vector EOS][Steel][Densities and temperatures] A Gamma(rho, T) lookup is performed"""
        gruneisen_eospac = np.zeros(self.num)
        gruneisen_spiner = np.zeros(self.num)
        num_scratch = int(self.eospac.scratch_size("GruneisenParamFromDensityTemperature", self.num) / np.double(1).nbytes)
        scratch = np.zeros(num_scratch)
        self.eospac.GruneisenParamFromDensityTemperature(self.density, self.temperature, gruneisen_eospac, scratch, self.num)
        self.steelEOS_host.GruneisenParamFromDensityTemperature(self.density, self.temperature, gruneisen_spiner, self.num)

        for a, b in zip(gruneisen_eospac, gruneisen_eospac):
            self.assertIsClose(a, b)

    def tearDown(self):
        # Failing to call finalize leads to a memory leak,
        # but otherwise behaviour is as expected.
        # It's possible to this automatically clean up with
        # some form of reference counting. If this is a priority,
        # we can re-examine.
        self.steelEOS_host.Finalize()
        self.eospac.Finalize()


class VectorEOS_Gruneisen_Given_Rho_Sie(unittest.TestCase):
    def setUp(self):
        """Parameters for a Gruneisen EOS"""

        # Unit conversions
        cm = 1.
        us = 1e-06
        Mbcc_per_g = 1e12

        # Gruneisen parameters for copper
        C0 = 0.394 * cm / us
        S1 = 1.489
        S2 = 0.
        S3 = 0.
        Gamma0 = 2.02
        b = 0.47
        rho0 = 8.93
        T0 = 298.
        P0 = 0.
        Cv = 0.383e-05 * Mbcc_per_g

        # Create the EOS
        self.eos = singularity_eos.Gruneisen(C0, S1, S2, S3, Gamma0, b, rho0, T0, P0, Cv)

        # Densities and energies
        self.num = 4
        self.density = np.zeros(self.num, dtype=np.double)
        self.energy = np.zeros(self.num, dtype=np.double)

        # Populate the input arrays
        self.density[0] = 8.0
        self.density[1] = 9.0
        self.density[2] = 9.5
        self.density[3] = 0.
        self.energy[0] = 1.e9
        self.energy[1] = 5.e8
        self.energy[2] = 1.e8
        self.energy[3] = 0.

        # Gold standard values for a subset of lookups
        self.pressure_true = np.array((-1.282094800000000e+11, 1.998504088912181e+10, 9.595823319513451e+10, P0 - C0 * C0 * rho0))
        self.bulkmodulus_true = np.array((9.990648504000005e+11, 1.460692677162573e+12, 1.851227213843747e+12, Gamma0 * (P0 - C0 * C0 * rho0)))
        self.temperature_true = np.array((5.590966057441253e+02, 4.285483028720627e+02, 3.241096605744125e+02, T0))
        self.gruneisen_true = np.array((Gamma0, 2.007944444444444e+00, 1.927000000000000e+00, Gamma0))

        self.temperature = np.zeros(self.num, dtype=np.double)
        self.pressure = np.zeros(self.num, dtype=np.double)
        self.bulkmodulus = np.zeros(self.num, dtype=np.double)
        self.gruneisen = np.zeros(self.num, dtype=np.double)

    def test_pressure(self):
        """[VectorEOS][GruneisenEOS][Densities and Energies] A P(rho, e) lookup is performed"""
        self.eos.PressureFromDensityInternalEnergy(self.density, self.energy, self.pressure, self.num)
        assert_allclose(self.pressure, self.pressure_true, rtol=1e-12)

    def test_bmod(self):
        """[VectorEOS][GruneisenEOS][Densities and Energies] A B_S(rho, e) lookup is performed"""
        self.eos.BulkModulusFromDensityInternalEnergy(self.density, self.energy, self.bulkmodulus, self.num)
        assert_allclose(self.bulkmodulus, self.bulkmodulus_true, rtol=1e-12)

    def test_temp(self):
        """[Vector EOS][GruneisenEOS][Densities and Energies] A T(rho, e) lookup is performed"""
        self.eos.TemperatureFromDensityInternalEnergy(self.density, self.energy, self.temperature, self.num)
        assert_allclose(self.temperature, self.temperature_true, rtol=1e-12)

    def test_gamma(self):
        """[Vector EOS][GruneisenEOS][Densities and Energies] A Gamma(rho, e) lookup is performed"""
        self.eos.GruneisenParamFromDensityInternalEnergy(self.density, self.energy, self.gruneisen, self.num)
        assert_allclose(self.gruneisen, self.gruneisen_true, rtol=1e-12)


class AluminumGruneisenEOS_SoundSpeedAndPressureComp(unittest.TestCase, EOSTestBase):
    "Aluminum Gruneisen EOS Sound Speed and Pressure Comparison"

    def setUp(self):
        """Parameters for a Gruneisen EOS"""

        # Unit conversions
        self.mm = 10.
        self.cm = 1.
        self.us = 1.e-06
        self.Mbar = 1.e12
        self.Mbcc_per_g = 1e12

        # Gruneisen parameters for copper
        self.C0 = 0.535 * self.cm / self.us
        self.S1 = 1.34
        self.S2 = 0.
        self.S3 = 0.
        self.Gamma0 = 1.97
        self.b = 0.
        self.rho0 = 2.714000
        self.T0 = 298.
        self.P0 = 1e-06 * self.Mbar
        self.Cv = 0.383e-05 * self.Mbcc_per_g

        # Create the EOS
        self.eos = singularity_eos.Gruneisen(self.C0, self.S1, self.S2, self.S3, self.Gamma0, self.b, self.rho0, self.T0, self.P0, self.Cv)

        self.density = 5.92418956756592            # g/cm^3
        self.energy = 792486007.804619             # erg/g
        self.true_pres = 2.620656373250729         # Mbar
        self.true_sound_speed = 1.5247992468363685 # cm/us

    def test_pressure(self):
        "[GruneisenEOS] A P(rho, e) lookup is performed"
        pres = self.eos.PressureFromDensityInternalEnergy(self.density, self.energy)
        pres = pres / self.Mbar
        self.assertIsClose(pres, self.true_pres, 1e-12)

    def test_bmod(self):
        "[GruneisenEOS] A B_S(rho, e) lookup is performed"
        bulk_modulus =  self.eos.BulkModulusFromDensityInternalEnergy(self.density, self.energy)
        sound_speed = math.sqrt(bulk_modulus / self.density) / (self.cm / self.us)
        self.assertIsClose(sound_speed, self.true_sound_speed, 1e-12)

    def test_bmod_approx(self):
        """A finite difference approximation is used for the bulk modulus"""
        # Bulk modulus approximation:
        #  B_S = rho * dPdr_e + P / rho * dPde_r
        drho = 1e-06 * self.density
        de = 1e-06 * self.energy
        P1 = self.eos.PressureFromDensityInternalEnergy(self.density, self.energy)
        P2 = self.eos.PressureFromDensityInternalEnergy(self.density + drho, self.energy)
        dPdr_e = (P2 - P1) / drho
        P2 = self.eos.PressureFromDensityInternalEnergy(self.density, self.energy + de)
        dPde_r = (P2 - P1) / de
        bmod_approx = self.density * dPdr_e + P1 / self.density * dPde_r
        # The finite difference solution should approximate the exact solution
        bulk_modulus = self.eos.BulkModulusFromDensityInternalEnergy(self.density, self.energy)
        ss_approx = math.sqrt(bmod_approx / self.density)
        sound_speed = math.sqrt(bulk_modulus / self.density)
        self.assertIsClose(sound_speed, ss_approx, 1e-5)

    def test_particle_vel_and_hugoniot(self):
        """A particle velocity and the same Hugoniot fit used in the EOS"""
        # Use Rankine-Hugoniot jump conditions to calculate a consistent point on the EOS
        up = 2 * self.cm / self.us # 20 km/s is a pretty strong shock
        Us = self.C0 + self.S1 * up
        e0 = 0
        # We have the density, energy, and pressure at this point") {
        density = Us * self.rho0 / (Us - up)
        true_pres = self.P0 + self.rho0 * Us * up
        energy =  e0 + 1. / 2. * Us * up * (1 - self.rho0 / density) + self.P0 * (1 / self.rho0 - 1 / density)
        # A P(rho, e) lookup is performed for the Hugoniot energy and density
        pres = self.eos.PressureFromDensityInternalEnergy(density, energy)
        self.assertIsClose(pres, true_pres, 1e-12)

if __name__ == "__main__":
    unittest.main()
