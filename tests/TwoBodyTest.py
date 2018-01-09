
import HaPPPy
import unittest

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

if __name__ == '__main__':
    two_body_suite = unittest.TestLoader().loadTestsFromTestCase(TwoBodyTestSuite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(two_body_suite)
