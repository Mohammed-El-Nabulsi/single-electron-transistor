#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import HaPPPy
import unittest
from tests.HaPPPyTest import HappyBasicTestSuite
from tests.OneBodyTest import OneBodyTestSuite
from tests.TwoBodyTest import TwoBodyTestSuite
from tests.TransmissionTest import TransmissionTestSuite
from tests.RatesTest import RatesTestSuite
from tests.MasterEquationTest import MasterEquationTestSuite

if __name__ == '__main__':
    # happy_suite = unittest.TestLoader().loadTestsFromTestCase(HappyBasicTestSuite)
    one_body_suite = unittest.TestLoader().loadTestsFromTestCase(OneBodyTestSuite)
    two_body_suite = unittest.TestLoader().loadTestsFromTestCase(TwoBodyTestSuite)
    transmission_suite = unittest.TestLoader().loadTestsFromTestCase(TransmissionTestSuite)
    rates_suite = unittest.TestLoader().loadTestsFromTestCase(RatesTestSuite)        
    master_equation_suite = unittest.TestLoader().loadTestsFromTestCase(MasterEquationTestSuite)
    # unittest.TextTestRunner(verbosity=2, buffer=True).run(happy_suite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(one_body_suite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(two_body_suite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(transmission_suite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(rates_suite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(master_equation_suite)
