import unittest
import singularity_eos

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

if __name__ == "__main__":
    unittest.main()
