
import HaPPPy
import unittest
import h5py
import numpy as np
import matplotlib.pyplot as plt
from HaPPPy.TwoBody.OneParticleLoader import SpectrumData
from HaPPPy.TwoBody.TwoParticle import createTwoParticleData
import HaPPPy.TwoBody.MatrixElement as MatrixElement

# general settings
_debug = False # whether to log debug info
_useExistingData = False # whether to use data that was generated in a previous test instead of recalculating everything
_filePath = 'tests/data/twoBody.h5py' # file used to store or load test data
_obDataFile='data_group1' # file which provides the used one body data

# opens the hdf5 file for storing test data
if _useExistingData:
    repository = h5py.File(_filePath, "r")
else:
    repository = h5py.File(_filePath, "w")

def check_convergence(Val):
    """helper function to test whether a given array of test data converges to a certain value.
    
    raises an AssertionError if data is not converging
    return lowest and highest convergence factor from one element to the next.
    Values between -1 and +1 indicate convergence: close to 0 = fast, close to +/-1 = slow
    Values below -1 or above +1 indicate divergence: close to +/-1 = slow, close to +/-infinite = fast
    Positive indicates monotone whereas negative indicates oscillating convergence/divergence
    """
    dE = np.ediff1d(Val)
    L_max = -1.0
    L_min = 1.0
    for i in range(len(dE) - 1):
        e = dE[i + 1]
        if (abs(e / Val[i]) < 1e-6): break
        L = e / dE[i]
        L_max = max(L_max, L)
        L_min = min(L_min, L)
    return L_min, L_max


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
        # TBSolver = HaPPPy.TwoBody.TwoBodySolver()
        # self.assertIsNotNone(TBSolver.doCalculation('data_group1'))

    def test_opWaves_count_convergence(self):
        """ Checks whether the two-body base energy converges when increasing the amount of feed in one-body states.
        """
        # settings
        _step = 2.0
        _max = 20.0
        
        # generate / load test data
        if (_useExistingData):
            N = repository['varObCount/count'][:]
            E = repository['varObCount/energy'][:]
        else:
            obData = SpectrumData()
            obData.open(_obDataFile)
            N = np.arange(_step, min(obData.m, _max), _step)
            it = np.nditer([N, None])
            for n, e in it:
                if _debug: print("running with %d ob-states:" % n)
                obData.m = int(n)
                ev = createTwoParticleData(obData)[0]
                if _debug: print(ev)
                e[...] = ev[0]
            E = it.operands[1]
            repository.create_dataset('varObCount/count', data=N, dtype='i')
            repository.create_dataset('varObCount/energy', data=E, dtype='d')
            repository.flush()
        
        # evaluate data
        if _debug: print("summary of two-body base energies:")
        if _debug: print(E)
        L_min, L_max = check_convergence(E)
        print("Converges with %.2f < factor < %.2f per %d added wavefunctions" % (L_min, L_max, _step))
        self.assertLess(L_max, 1.0, msg='diverging!')
        self.assertGreater(L_min, -1.0, msg='diverging!')

    def test_coulomb_epsilon_convergence(self):
        """ Checks whether the two-body energy states converge when decreasing epsilon (the coulomb interaction cutoff)
        """
        #TODO write this test

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
        self.assertAlmostEqual(loader1.l, 5, msg='wrong grid size!')
        self.assertAlmostEqual(loader1.n * loader1.dx, loader1.l, msg='inconsistent grid size!')
        loader1.energies[:] = energies
        loader1.waves[:,0] = wave0
        loader1.waves[:,1] = wave1
        loader1.close()
        # load from file
        loader2 = SpectrumData()
        loader2.open(path)
        # compare data
        self.assertEqual(loader1.m, loader2.m, msg='wrong amount of energies!')
        self.assertEqual(loader1.n, loader2.n, msg='wrong amount of grid points!')
        self.assertAlmostEqual(loader1.l, loader2.l, msg='corrupted grid width!')
        self.assertAlmostEqual(loader1.dx, loader2.dx, msg='corrupted grid size!')
        self.assertTrue(np.allclose(loader2.energies[:], energies), msg='corrupted eigenenergy values')
        self.assertTrue(np.allclose(loader2.waves[:,0], wave0), msg='corrupted wavefunction data')
        self.assertTrue(np.allclose(loader2.waves[:,1], wave1), msg='corrupted wavefunction data probably too low precision')
        # test normalisation
        waves = loader2.getNormalizedWaves()
        for i in range(loader2.m):
            wave = waves[:,i]
            self.assertEqual(len(wave), loader2.n, msg='result has wrong dimensions!')
            self.assertAlmostEqual(np.inner(wave, wave) * loader2.dx, 1.0, msg='incorrectly normalized!')
        loader2.close()

if __name__ == '__main__':
    two_body_suite = unittest.TestLoader().loadTestsFromTestCase(TwoBodyTestSuite)
    unittest.TextTestRunner(verbosity=2, buffer=False).run(two_body_suite)
