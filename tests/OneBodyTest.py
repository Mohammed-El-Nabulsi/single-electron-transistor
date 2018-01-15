
import HaPPPy
import unittest

class OneBodyTestSuite(unittest.TestCase):
    """A test class to the OneBody module.

    """

    def test_OneBody_exists(self):
        """ Checks wether the One Body module exists
        """
        self.assertTrue(hasattr(HaPPPy, 'OneBody'))

    def test_OneBody_doCalculation(self):
        """ Checks the dummy calculation
        """  
        print("hi 2")
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 100)
        self.assertEqual(2.0, 2.0)

    def test_OneBody_doCalculation(self):
        """ Checks the dummy calculation
        """
        print("hi")
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 100)
        self.assertEqual(42.0, 42.0)

    def test_OneBody_calcualteHarmonocPotential(self):
        """ Checks the dummy calculation
        """  
        print("hi 2")
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 100)
        la, vabcd, i = OBSolver.calcualteHarmonocPotential(0)
        self.assertTrue((la[0] - la[1]) - ( la[1] - la[2]) < 1e-5)

if __name__ == '__main__':
    # OBSolver = HaPPPy.OneBody.OneBodySolver(100,100)
    # la, vabcd, i = OBSolver.calcualteHarmonocPotential(0)
    # OBSolver.exportData(la, vabcd, i)
    # print(vabcd)

   one_body_suite = unittest.TestLoader().loadTestsFromTestCase(OneBodyTestSuite)
   unittest.TextTestRunner(verbosity=2, buffer=True).run(one_body_suite)
