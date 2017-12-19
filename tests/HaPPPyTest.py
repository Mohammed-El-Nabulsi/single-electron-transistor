#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import HaPPPy
import unittest

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

if __name__ == '__main__':
    happy_suite = unittest.TestLoader().loadTestsFromTestCase(HappyBasicTestSuite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(happy_suite)
