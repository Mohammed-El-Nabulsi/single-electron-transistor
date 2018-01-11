
import HaPPPy
import unittest
import numpy as np
from HaPPPy.TwoBody.OneParticleLoader import SpectrumData

class TwoBodyTestSuite(unittest.TestCase):
    """A test class to the TwoBody module.

    """

    def test_TwoBody_exists(self):
        """ Checks wether the One Body module exists
        """
        self.assertTrue(hasattr(HaPPPy, 'TwoBody'))

    def test_TwoBody_doCalculation(self):
        """ Checks the dummy calculation
        """  
        TBSolver = HaPPPy.TwoBody.TwoBodySolver()
        self.assertEqual(TBSolver.doCalculation(), 2.0)

    def test_oneParticleLoader(self):
        """ Checks the OneParticleLoader for self consistency
        """
        # create test data
        path = "/tmp/oneParticleTest"
        energies = np.array([1.0, 2.0])
        wave0 = np.array([0.0, 0.5, 0.9, 1.0, 1e40, -1e-40, -1.0, -0.9, -0.5, -0.0])
        wave1 = np.array([-0.0, -0.5, -0.9, -1.0, -2.0, 2.0, 1.0, 0.9, 0.5, 0.0])
        # save to file
        loader1 = SpectrumData()
        loader1.init(path, 2, 10, L=5)
        self.assertAlmostEqual(loader1.n * loader1.dx, 5, msg='wrong grid size!')
        self.assertAlmostEqual((loader1.n - 1) * loader1.dx, loader1.x1 - loader1.x0, msg='wrong grid bounds!')
        loader1.energies[:] = energies
        loader1.waves[0,:] = wave0
        loader1.waves[1,:] = wave1
        loader1.close()
        # load from file
        loader2 = SpectrumData()
        loader2.open(path)
        # compare data
        self.assertAlmostEqual(loader1.x0, loader2.x0, msg='corrupted grid bounds!')
        self.assertAlmostEqual(loader1.x1, loader2.x1, msg='corrupted grid bounds!')
        self.assertAlmostEqual(loader1.dx, loader2.dx, msg='corrupted grid size!')
        self.assertEqual(loader1.m, loader2.m, msg='wrong dataset dimensions!')
        self.assertEqual(loader1.n, loader2.n, msg='wrong dataset dimensions!')
        self.assertTrue(np.allclose(loader2.energies[:], energies), msg='corrupted eigenenergy values')
        self.assertTrue(np.allclose(loader2.waves[1,:], wave0), msg='corrupted wavefunction data')
        self.assertTrue(np.allclose(loader2.waves[0,:], wave1), msg='corrupted wavefunction data probably too low precision')
        # test normalisation
        waves = loader2.getNormalizedWaves()
        for i in range(loader2.m):
            wave = waves[i,:]
            self.assertEqual(len(wave), loader2.n, msg='result has wrong dimensions!')
            self.assertAlmostEqual(np.inner(wave, wave) * loader2.dx, 1.0, msg='incorrectly normalized!')
        loader2.close()

if __name__ == '__main__':
    two_body_suite = unittest.TestLoader().loadTestsFromTestCase(TwoBodyTestSuite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(two_body_suite)
