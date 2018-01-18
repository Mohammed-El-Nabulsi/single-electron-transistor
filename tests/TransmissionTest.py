
import HaPPPy
import unittest

class TransmissionTestSuite(unittest.TestCase):
    """A test class to the Transmission module.

    """

    def test_Transmission_exists(self):
        """ Checks wether the Transmission module exists
        """
        self.assertTrue(hasattr(HaPPPy, 'Transmission'))

    def test_Transmission_doCalculation(self):
        """ Checks the dummy calculation
        """  
        Calc = HaPPPy.Transmission.TransmissionCalculator()
        self.assertEqual(Calc.doCalculation(), 2.0)

if __name__ == '__main__':
    transmission_suite = unittest.TestLoader().loadTestsFromTestCase(TransmissionTestSuite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(transmission_suite)
