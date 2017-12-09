#!/usr/bin/env python3
# coding: utf-8

"""
.. module:: HaPPPy
   :platform: Unix, Windows
   :synopsis: The HaPPPy module to simulate quantum dots.

.. moduleauthor:: Lars-Hendrik Frahm <lfrahm@physnet.uni-hamburg.de>


The HaPPPy program
==================

The **Ha**\ mburg **P**\ rogrammier **P**\ rojekt **Py**\ thon is the program written by the 
students of the course *66-527 Proseminar: Programmierprojekt* at the University of Hamburg. 
The program serves to simulate various properties of quantum dots.

"""

import sys
import datetime
import socket
import getpass

from HaPPPy.OneBody import OneBodySolver
from HaPPPy.TwoBody import TwoBodySolver
from HaPPPy.Transmission import TransmissionCalculator
from HaPPPy.Rates import RateCalculator
from HaPPPy.MasterEquation import MasterEquationSolver

__version__ = "0.0.1"

"""The main function of HaPPPy, which executed the code. 

Right now, HaPPPy does not require any command line variables
"""

def main(argv=None):
    """The main function of HaPPPy, which executed the code. 

    Args:
       argv, the command line arguments.

    Returns:
       int.  The return code::

          0 -- Success!
          1 -- No good.

    Raises:
       RuntimeError

    Use me as

    >>> import HaPPPy
    >>> HaPPPy.main()
    0

    .. warning::

      This function requires Python version 3.0.0 or higher. 
      It raises a RuntimeError when used with Python 2 interpreter.

    """

    print("################################################################################")
    print("###")
    print("### Hello and welcome to HaPPPy! " + __version__)
    print("###")
    print("### Happy is the program developed in the \"Programmier Projekt 2017/2018\"")
    print("### of the University of Hamburg ")
    print("###")
    print("################################################################################")
    print("###")
    print("### Your are running Happy at " + str(datetime.datetime.now()) + " as " + getpass.getuser() + "@" + socket.gethostname())
    print("### Your Python version is " + str(sys.version_info[0]) + "." + str(sys.version_info[1]) + "." + str(sys.version_info[2]))

    if(sys.version_info[0] < 3):
        print("### Your Python version is to low! Use Pathon 3!")
        raise RuntimeError('version', 'You need Python version 3 to run HaPPPy')

    print("###")
    print("################################################################################")

    OBSolver = OneBodySolver()
    OBSolver.doCalculation()

    TBSolver = TwoBodySolver()
    TBSolver.doCalculation()

    TMCal = TransmissionCalculator()
    TMCal.doCalculation()

    RateCal = RateCalculator()
    RateCal.doCalculation()

    MSolver = MasterEquationSolver()
    MSolver.doCalculation()

    return 0