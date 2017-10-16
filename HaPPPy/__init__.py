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

__version__ = "0.0.1"

def main(argv=None):
    """
    The main function of HaPPPy, which executed the code. 

    Right now, HaPPPy does not have any comand line variables
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
        raise RuntimeError('version', 'You need Python version 3 to run HaPPPy')
    print("###")
    print("################################################################################")