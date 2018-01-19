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
import json
import h5py

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

    # parse passed commandline arguments
    verbose = "-v" in argv
    if verbose: print("Running in verbose mode."
                     +"This is recommended for debuging only!")

    # Loading config from json file
    configdata = json.load(open('config.json'))

    L = configdata["L"]
    N = configdata["N"]
    E0 = configdata["E0"]

    # Do the one body calculation
    # TODO Group 1: Fix spelling mistakes
    OneBody = OneBodySolver(L,N)
    OneBodyEigenvalues, OneBodyEigenVectors, Info = OneBody.calcualteHarmonicPotential(E0)

    # Save file to hdf5 file
    # TODO Group 1: Allow to choose the path for hdf5
    # TODO Group 1: Also allow to importData from hdf5
    OneBody.exportData(OneBodyEigenvalues, OneBodyEigenVectors, Info)

    # TODO Group 2: Put all this into a class, say "TwoBodySolver"
    # TODO Group 2: Fix the interface with Group 1 in an elegant way
    # TODO Group 2: Make sure the module works as expected with the interface provided
    # TODO Group 2: Try to enhance performance by using BLAS routines
    from HaPPPy.TwoBody.OneParticleLoader import SpectrumData
    onebodydatapath = 'data_group1'
    obData = SpectrumData()
    obData.open(onebodydatapath)
    TwoBodyEigenvalues, TwoBodyEigenvectors = createTwoParticleData(obData)

    # Save the two body result
    # TODO Group 2: Write a function that allows to store the results in a hdf5 file as part of "TwoBodySolver"
    twobodydatapath = 'data_group2.hdf5'
    dataFile = h5py.File(twobodydatapath, "w")
    dataSet_calcInfo = dataFile.create_dataset("TwoBodyEigenvalues", data=TwoBodyEigenvalues)
    dataSet_eigvalues = dataFile.create_dataset("TwoBodyEigenvectors", data=TwoBodyEigenvectors)  
    dataFile.close()

    # TODO Group 2 or Group 3: Read two body data from hdf5
    file = h5py.File('data_group2.hdf5', "a")
    TwoBodyEigenvalues = file["TwoBodyEigenvalues"][:]
    TwoBodyEigenvectors = file["TwoBodyEigenvectors"][:]
    file.close()

    # TODO Group 3 or Group 4: Make an interface between both modules. 
    #                          How do the transmission coefficients get into the RateCalculator
    #                          Maybe think of function pointers
    Transmission = TransmissionCalculator()

    muL = configdata["muL"]
    muR = configdata["muR"]
    T = 1.0
    V = 1.0 # What is V actually?

    # TODO Group 4: There seems to be a problem with the Fermin function
    # TODO Group 4 and Group 5: Fix the interface between both modules
    #                           There is no common way to transfer the rates 
    RateCal = RateCalculator()
    Gamma_L, Gamma_R = Rates.doCalculation(OneBodyEigenvalues, TwoBodyEigenvalues,\
                                                                  muL, muR, T, V, TwoBodyEigenvectors, Transmission, DensityofStates)

    # Getting some parameters
    TMax = configdata["Tmax"]
    DT = configdata["DT"]

    # TODO Group 5: Allow to calculate stationary conductance
    # TODO Group 5: OPTIONAL Think how do we get the initial states from the config file in here
    #               Maybe divide into specific QD states. Like "empty", "double occupied", ...
    n = Gamma_L.shape[0]
    # choose a legitimate start value for P_0 (currently a dummy!)
    P_0 = np.arange(n)
    P_0 = 1 / sum(P_0) * P_0
    mes = HaPPPy.MasterEquation.MasterEquationSolver()
    sim_tdp, sim_cur = mes.simulateDynamicSloution(DT, TMax, P_0, Gamma_L, Gamma_R)
    stat_ps, stat_curs = mes.calculateStationarySloutions(Gamma_L, Gamma_R)
    ## plot/print results
    # 1st plot: time development of propabilities
    sim_tdp.quickPlot(xlabel="t", ylabel="P")
    print("static solutions (P_ij) =\n", stat_ps)
    # 2nd plot: time development of netto current
    sim_cur.quickPlot(xlabel="t", ylabel="I",legend=["$I^L$","$I^R$"])
    print("static solutions (I_i) = \n", stat_curs)

    return 0
