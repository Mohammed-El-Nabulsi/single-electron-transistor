# The following imports where disabled to prevent runtime errors caused by
# double imports with Python 3.6.x or higher. Since all listed classes are
# part of test_main.py - which is the file containing the program entry point -
# there is no need for these explicit imports.
# Details and recreation instructions:
# The error orrured when running the tests module via the command
# python3 -m tests.test_main
# with Python 3.6.3 . The error message was:
# "... RuntimeWarning: 'tests.test_main' found in sys.modules after import of
# package 'tests', but prior to execution of 'tests.test_main'; this may result
# in unpredictable behaviour."
"""
from .test_main import HappyBasicTestSuite
from .test_main import OneBodyTestSuite
from .test_main import TwoBodyTestSuite
from .test_main import TransmissionTestSuite
from .test_main import RatesTestSuite
from .test_main import MasterEquationTestSuite
"""
