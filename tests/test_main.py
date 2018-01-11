#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import HaPPPy
import unittest
import numpy as np
import matplotlib.pyplot as plt

class HappyBasicTestSuite(unittest.TestCase):
    """A test class to test HaPPPy in general.

        This class is supposed to test the general behavior of HaPPPy. This includes test cases involving more than one 
        package, as well as input and output handling. 

    """

    def test_HaPPPy_main_exists(self):
        """ Checks wether the main function of HaPPPy exists
        """
        self.assertTrue(hasattr(HaPPPy, 'main'))

    def test_HaPPPy_main_runs(self):
        """ Checks wether HaPPPy runs without exceptions if no user input was made
        """
        try:
            HaPPPy.main()
        except Exception:
            self.fail("HaPPPy.main() raised Exception unexpectedly!")

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

class RatesTestSuite(unittest.TestCase):
    """A test class to the Rates module.

    """

    def test_Rates_exists(self):
        """ Checks wether the Rates module exists
        """
        self.assertTrue(hasattr(HaPPPy, 'Rates'))

    def test_Rates_doCalculation(self):
        """ Checks the dummy calculation
        """  
        Calc = HaPPPy.Rates.RateCalculator()
        self.assertEqual(Calc.doCalculation(), 2.0)

class MasterEquationTestSuite(unittest.TestCase):
    """A test class to the MasterEquation module.

    """

    def test_MasterEquation_exists(self):
        """ Checks wether the MasterEquation module exists
        """
        self.assertTrue(hasattr(HaPPPy, 'MasterEquation'))

    def test_MasterEquation_doCalculation(self, plot_figures=False):
        """ Checks the dummy calculation
        """
        ## test program for simulate_time_development_of_propabilities
        # set-up a reasonable Γ-matrix
        Γ_l = np.array([[0, 0.5, 0], [0, 0, 0], [0, 0, 0]])
        Γ_r = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
        # choose a legitimate start value for P_0 (P_0 = P(t=0))
        P_0 = np.array([0.9, 0.1, 0])
        # simulate

        mes = HaPPPy.MasterEquation.MasterEquationSolver()
        sim_tdp, sim_cur = mes.doCalculation(1, 100, P_0, Γ_r, Γ_l)
        Ps = sim_tdp.getValues()
        self.assertTrue(Ps[99][0] < 1e-10)

        ## plot result
        if plot_figures:
            sim_tdp.quickPlot()
            sim_cur.quickPlot()

if __name__ == '__main__':
    happy_suite = unittest.TestLoader().loadTestsFromTestCase(HappyBasicTestSuite)
    one_body_suite = unittest.TestLoader().loadTestsFromTestCase(OneBodyTestSuite)
    two_body_suite = unittest.TestLoader().loadTestsFromTestCase(TwoBodyTestSuite)
    transmission_suite = unittest.TestLoader().loadTestsFromTestCase(TransmissionTestSuite)
    rates_suite = unittest.TestLoader().loadTestsFromTestCase(RatesTestSuite)        
    master_equation_suite = unittest.TestLoader().loadTestsFromTestCase(MasterEquationTestSuite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(happy_suite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(one_body_suite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(two_body_suite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(transmission_suite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(rates_suite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(master_equation_suite)
