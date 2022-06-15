import unittest
import singularity_eos
import numpy as np
from numpy.testing import assert_allclose

class EOS(unittest.TestCase):
    def testConstants(self):
        from singularity_eos import thermalqs
        self.assertEqual(thermalqs.all_values, thermalqs.none |
                                               thermalqs.density |
                                               thermalqs.specific_internal_energy |
                                               thermalqs.pressure      |
                                               thermalqs.temperature   |
                                               thermalqs.specific_heat |
                                               thermalqs.bulk_modulus)

    def testIdealGas(self):
        eos = singularity_eos.IdealGas(1,1)

    def testGruneisen(self):
        eos = singularity_eos.Gruneisen(1,1,1,1,1,1,1,1,1,1)

    def testJWL(self):
        eos = singularity_eos.JWL(1,1,1,1,1,1,1)

    def testDavisReactants(self):
        eos = singularity_eos.DavisReactants(1,1,1,1,1,1,1,1,1,1,1)

    def testDavisProducts(self):
        eos = singularity_eos.DavisProducts(1,1,1,1,1,1,1,1)


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
class SpinerEOSdependsOnRhoT_Steel(unittest.TestCase):
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
        self.assertEqual(self.steelEOS_host.matid(), self.steelID)
        self.assertEqual(self.steelEOS_host.filename(), self.eosName)

    def test_reference(self):
        "We can get a reference density and temperature"
        state = self.steelEOS_host.ValuesAtReferenceState()
        state_pac = self.eospac.ValuesAtReferenceState()
        self.assertTrue(math.isclose(state.density, state_pac.density, rtol=1e-12))
        self.assertTrue(math.isclose(state.temperature, state_pac.temperature, rtol=1e-12))

    # TODO: this needs to be a much more rigorous test
    def test_quantities_from_rho_temp(self):
        "Quantities can be read from density and temperature"
        sie_pac = self.eospac.InternalEnergyFromDensityTemperature(1e0, 1e6)
        ie = self.steelEOS.InternalEnergyFromDensityTemperature(1e0, 1e6)
        self.assertTrue(math.isclose(ie, sie_pac, rtol=1e-12))

    def test_rho_of_P_T(self):
        "rho(P,T) correct for P=1atm, T=freezing"
        T = 273  # Kelvin
        P = 1e6  # barye
        lmbda = np.zeros(self.steelEOS_host.nlambda(), dtype=np.double)
        rho, sie = self.steelEOS_host.DensityEnergyFromPressureTemperature(P, T, lmbda)
        rho_pac, sie_pac = self.eospac.DensityEnergyFromPressureTemperature(P, T)
        self.assertTrue(math.isclose(rho, rho_pac, rtol=1e-12))

    def tearDown(self):
        # Failing to call finalize leads to a memory leak,
        # but otherwise behaviour is as expected.
        # It's possible to this automatically clean up with
        # some form of reference counting. If this is a priority,
        # we can re-examine.
        self.steelEOS_host.Finalize()
        self.eospac.Finalize()


@unittest.skipIf('SpinerEOSDependsRhoT' not in dir(singularity_eos) or 'EOSPAC' not in dir(singularity_eos), "No Spiner or EOSPAC support")
class SpinerEOSdependsOnRhoT_Air(unittest.TestCase):
    "SpinerEOS and EOSPAC for air can be initialized with matid"

    def setUp(self):
        self.eosName = "../materials.sp5"
        self.airEOS_host = singularity_eos.SpinerEOSDependsRhoT(self.eosName, self.airID)
        self.eospac = singularity_eos.EOSPAC(self.airID)

    def test_reference(self):
        "We can get a reference state"
        state = self.airEOS_host.ValuesAtReferenceState()
        state_pac = eospac.ValuesAtReferenceState()
        self.assertTrue(math.isclose(state.density, state_pac.density, rtol=1e-12))
        self.assertTrue(math.isclose(state.temperature, state_pac.temperature, rtol=1e-12))
   
    def test_P_from_rho_sie(self): 
        "P(rho, sie) correct for extrapolation regime"
        rho = 1
        sie = 2.43e16
        P_pac = self.eospac.PressureFromDensityInternalEnergy(rho, sie)
        P_spi = self.airEOS_host.PressureFromDensityInternalEnergy(rho, sie)
        self.assertTrue(math.isclose(P_pac, P_spi, rtol=1e-12))
    
    def tearDown(self):
        self.airEOS_host.Finalize()


@unittest.skipIf('SpinerEOSDependsRhoT' not in dir(singularity_eos) or 'EOSPAC' not in dir(singularity_eos), "No Spiner or EOSPAC support")
class EOS_init_with_matid2(unittest.TestCase):
    "EOS initialized with matid"

    def setUp(self):
        self.eosName = "../materials.sp5"
        self.DTID = 5267
        self.eos_spiner = singularity_eos.SpinerEOSDependsRhoT(self.eosName, self.DTID)
        self.eos_eospac = singularity_eos.EOSPAC(self.DTID)

    def test_inversion_for_T_rho_P(self):
        "Inversion for T(rho,P) works on host"
        P = 1e8
        rho = 1.28e-3
        output = (thermalqs.temperature | thermalqs.specific_internal_energy | thermalqs.specific_heat | thermalqs.bulk_modulus)
        state = self.eos_spiner.FillEos(rho=rho, press=P, output=output)
        state_pac = self.eos_eospac.FillEos(rho=rho, press=P, output=output)
        self.assertTrue(math.isclose(state.temperature, state_pac.temperature, rtol=1e-12))
        self.assertTrue(math.isclose(state.specific_internal_energy, state_pac.specific_internal_energy, rtol=1e-12))
        self.assertTrue(math.isclose(state.specific_heat, state_pac.specific_heat, rtol=1e-12))
    
    def tearDown(self):
        eos_spiner.Finalize()
        eos_eospac.Finalize()


@unittest.skipIf('SpinerEOSDependsRhoT' not in dir(singularity_eos) or 'EOSPAC' not in dir(singularity_eos), "No Spiner or EOSPAC support")
class EOS_init_with_matid(unittest.TestCase):
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
        lmbda = np.zeros(eos_spiner.nlambda(), dtype=np.double)
        rho, sie = self.eos_spiner.DensityEnergyFromPressureTemperature(P, T, lmbda.data())
        rho_pac, sie_pac = self.eos_eospac.DensityEnergyFromPressureTemperature(P, T, lmbda.data())
        self.assertTrue(math.isclose(rho, rho_pac))

    def tearDown(self):
        self.eos_spiner.Finalize()

if __name__ == "__main__":
    unittest.main()
