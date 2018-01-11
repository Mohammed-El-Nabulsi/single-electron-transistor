
import HaPPPy
import unittest
import numpy as np

Energie_Einteilchen =  np.array([0.1, 3, 4]) #Energieeigenwerte von Gruppe 1
Vinput = [0, 1, 3, 6] #Vorgegebene Tunnelbarriere
Energie_Zweiteilchen = np.array([3, 4, 3.6]) #Energieegenwerte von Gruppe 2
Koeeffizienten_Zweiteilchen= np.array([[[0, 1, 0],[1, 2, 0],[0, 0, 1]],[[1, 0 ,0],[0.2, 0.5, 0],[0.1, 0.3, 1]],[[0, 0, 0.5],[1, 0, 0.1],[0, 1, 0.1]]]) #Koeeffizientenmatrizen von Gruppe 2
muL = 78
muR= 95
T= 270

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
        self.assertEqual(Calc.doCalculation(Energie_Einteilchen, Energie_Zweiteilchen, muL, muR, T, Vinput, Koeeffizienten_Zweiteilchen), 2.0)

if __name__ == '__main__':
    rates_suite = unittest.TestLoader().loadTestsFromTestCase(RatesTestSuite)        
    unittest.TextTestRunner(verbosity=2, buffer=True).run(rates_suite)
