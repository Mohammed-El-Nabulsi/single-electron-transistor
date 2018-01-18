
import HaPPPy
import unittest
import numpy as np

Energie_Einteilchen =  np.array([0.1, 3, 4]) #Energieeigenwerte von Gruppe 1
Vinput = [0, 1, 3, 6] #Vorgegebene Tunnelbarriere
Energie_Zweiteilchen = np.array([3, 4, 4, 5]) #Energieegenwerte von Gruppe 2
Koeeffizienten_Zweiteilchen= np.array([[[0, 1, 0, 0],[0, 1, 2, 0],[0, 0, 0, 1]],[[0, 1, 0 ,0],[0.5, 0.2, 0.5, 0],[0, 0.1, 0.3, 1]],[[0.5, 0, 0, 0.5],[1, 0, 0.1, 0.5],[0, 1, 0.1, 0.3]], [[0, 1, 0,3],[0.2,0.3,0,0.1],[0.1, 0.3, 0,1 ]]]) #Koeeffizientenmatrizen von Gruppe 2
muL = 10
muR= 2
T= 10

class DummyTransmission:

    def calculate_transmission(self, a, b):
        return 0.5

class DummyDensityofStates:

    def calculate_DensityofStates(self, A):
        return 0.5
    
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
        DummyTr = DummyTransmission()
        DummyDe = DummyDensityofStates()
        Calc.doCalculation(Energie_Einteilchen, Energie_Zweiteilchen, muL, muR, T, Vinput, Koeeffizienten_Zweiteilchen, DummyTr, DummyDe)
        self.assertEqual(Calc.doCalculation(Energie_Einteilchen, Energie_Zweiteilchen, muL, muR, T, Vinput, Koeeffizienten_Zweiteilchen, DummyTr, DummyDe), 2.0)	
	

if __name__ == '__main__':
    rates_suite = unittest.TestLoader().loadTestsFromTestCase(RatesTestSuite)        
    unittest.TextTestRunner(verbosity=2, buffer=True).run(rates_suite)
