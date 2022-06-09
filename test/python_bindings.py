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
        Cv = 5.0;
        gm1 = 0.4;
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


if __name__ == "__main__":
    unittest.main()
