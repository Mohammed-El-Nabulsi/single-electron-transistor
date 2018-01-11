
import HaPPPy
import unittest

class RatesTestSuite(unittest.TestCase):
    """A test class to the Rates module.

    """

    def test_Rates_exists(self):
        """ Checks if the Rates module exists
        """
        self.assertTrue(hasattr(HaPPPy, 'Rates'))

    def test_Rates_doCalculation(self):
        """ Checks the dummy calculation
        """  
        Calc = HaPPPy.Rates.RateCalculator()
        self.assertEqual(Calc.doCalculation(), 2.0)

if __name__ == '__main__':
    rates_suite = unittest.TestLoader().loadTestsFromTestCase(RatesTestSuite)        
    unittest.TextTestRunner(verbosity=2, buffer=True).run(rates_suite)
