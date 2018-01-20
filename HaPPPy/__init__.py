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


    ## Group 5: MASTER MODULE
    if verbose: print("\nSTEP 5 : MASTER EQUATION\n")
    # TODO Group 5: Allow to calculate stationary conductance
    # Get all relevant parameters from config file and previous calculations.
    n_one = len(OneBodyEigenvalues)
    n_two = len(TwoBodyEigenvalues)
    n_tot = 1 + n_one + n_two
    ns = [1, n_one, n_two]
    BCond = configdata["BCond"]
    Epsilon = configdata["Epsilon"]
    if Epsilon == "default":
        Epsilon = None
    if verbose: print("Epsilon = " + str(Epsilon))
    if verbose: print("n_one = " + str(n_one))
    if verbose: print("n_two = " + str(n_two))
    if verbose: print("n_tot = " + str(n_tot))
    if verbose: print("BCond = " + str(BCond))
    if verbose: print("Gamma_L =\n" + str(Gamma_L))
    if verbose: print("Gamma_R =\n" + str(Gamma_R))

    # Store input values to preserve the context.
    datapath_master = 'data/group5.hdf5'
    with h5py.File(datapath_master, "w") as f:
        f.create_dataset("Epsilon", data=Epsilon)
        f.create_dataset("n_one", data=n_one)
        f.create_dataset("n_two", data=n_two)
        f.create_dataset("n_tot", data=n_tot)
        f.create_dataset("BCond", data=BCond)
        f.create_dataset("Gamma_L", data=Gamma_L)
        f.create_dataset("Gamma_R", data=Gamma_R)

    # Calculate static or dynamic solution(s).
    mes = HaPPPy.MasterEquation.MasterEquationSolver()
    if BCond == "static":
        # Static solutions are requested.
        if verbose: print("mode = STATIC")
        stat_ps, stat_curs = mes.calculateStationarySloutions(Gamma_L, Gamma_R,
                                                              ns,
                                                              verbose=verbose,
                                                             )
        print("(P_stat_ij) =\n", stat_ps)
        print("(I^k_stat_i) = \n", stat_curs)
        # Store output values.
        with h5py.File(datapath_master, "w") as f:
            f.create_dataset("stat_ps", data=stat_ps)
            f.create_dataset("stat_curs", data=stat_curs)
    else:
        # A simulation is requested.
        if verbose: print("mode = DYNAMIC")
        # Get additional parameters for simulation.
        DT = configdata["DT"]
        TMax = configdata["TMax"]
        if verbose: print("DT = " + str(DT))
        if verbose: print("TMax = " + str(TMax))
        # Calculate a legitimate value for P_0 based on the input:
        P_0 = np.zeros(n_tot)
        if BCond == "empty":
            P_0[0] = 1
        elif BCond == "one":
            P_0[1] = 1
        elif BCond == "two":
            P_0[1 + n_one] = 1
        else:
            P_0 = np.array(BCond)
        P_0 = P_0 / sum(P_0)
        if verbose: print("P_0 =\n" + str(P_0))
        # Store P_0 to preserve context.
        with h5py.File(datapath_master, "w") as f:
            f.create_dataset("P_0", data=P_0)
        # Calculate the simulation.
        sim_p, sim_cur = mes.simulateDynamicSloution(DT, TMax,
                                                     P_0,
                                                     Gamma_L, Gamma_R,
                                                     ns,
                                                     verbose=verbose,
                                                    )
        # Plot the simulations.
        # 1st plot: time development of propabilities
        sim_p.quickPlot(xlabel="t", ylabel="P",
                        xunit_symbol="ps",
                       )
        # 2nd plot: time development of netto current
        sim_cur.quickPlot(xlabel="t", ylabel="I"
                          xunit_symbol="ps", yunit_symbol="$e$",
                          legend=["$I^L$","$I^R$"],
                         )
        # Store output values.
        with h5py.File(datapath_master, "w") as f:
            f.create_dataset("sim_p", data=sim_p)
            f.create_dataset("sim_cur", data=sim_cur)

    return 0
