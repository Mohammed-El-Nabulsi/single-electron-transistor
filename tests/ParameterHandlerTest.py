#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import HaPPPy
import unittest

class ParameterHandlerTestSuite(unittest.TestCase):
    """A test class to test HaPPPy in general.

        This class is supposed to test the general behavior of HaPPPy.
        This includes test cases involving more than one 
        package, as well as input and output handling. 

    """

    def test_HaPPPy_ParameterHandler_exists(self):
        """ Checks wether the ParameterHandler of HaPPPy exists
        """
        self.assertTrue(hasattr(HaPPPy, 'ParameterHandler'))

    def test_HaPPPy_ParameterHandler_SetGetTemperature(self):
        """ Checks the set/get method of the ParameterHandler
        """
        parameter = HaPPPy.ParameterHandler()
        parameter.setTemperature(3.0)
        self.assertEqual(parameter.getTemperature(), 3.0)

        self.assertRaises(TypeError, parameter.setTemperature, ["3"])

if __name__ == '__main__':
    parameter_suite = unittest.TestLoader().loadTestsFromTestCase(ParameterHandlerTestSuite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(parameter_suite)
