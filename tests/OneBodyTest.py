
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
        OBSolver = HaPPPy.OneBody.OneBodySolver()
        self.assertEqual(OBSolver.doCalculation(), 2.0)


if __name__ == '__main__':
    one_body_suite = unittest.TestLoader().loadTestsFromTestCase(OneBodyTestSuite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(one_body_suite)
