
import HaPPPy
import unittest

class MasterEquationTestSuite(unittest.TestCase):
    """A test class to the MasterEquation module.

    """

    def test_MasterEquation_exists(self):
        """ Checks wether the MasterEquation module exists
        """
        self.assertTrue(hasattr(HaPPPy, 'MasterEquation'))

    def test_MasterEquation_doCalculation(self):
        """ Checks the dummy calculation
        """  
        Calc = HaPPPy.MasterEquation.MasterEquationSolver()
        self.assertEqual(Calc.doCalculation(), 2.0)

if __name__ == '__main__':
    master_equation_suite = unittest.TestLoader().loadTestsFromTestCase(MasterEquationTestSuite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(master_equation_suite)
