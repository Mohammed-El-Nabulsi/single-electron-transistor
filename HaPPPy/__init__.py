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
import argparse
import datetime
import socket
import getpass
import json
import h5py
import numpy as np

from HaPPPy.OneBody import OneBodySolver
from HaPPPy.TwoBody import TwoBodySolver
from HaPPPy.Transmission import TransmissionCalculator
from HaPPPy.Rates import RateCalculator
from HaPPPy.MasterEquation import MasterEquationSolver

from HaPPPy.TwoBody.TwoParticleLoader import TwoBodySpectrumData

__version__ = "0.0.1"

"""The main function of HaPPPy, which executed the code.

Right now, HaPPPy does not require any command line variables
"""

class DOS:

    def calculate_DensityofStates(self, energy):
        return 1.0

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

    # Constants.
    default_path_config = 'config.json'

    # Parse passed commandline arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose',
                        action="store_true",
                        help="print detailed information about each step of the calulations",
                       )
    parser.add_argument('-c', '--config',
                        type=str,
                        default=default_path_config,
                        metavar="<PATH TO CONFIG FILE>",
                        dest="path_config",
                        help="defines path to the config file (see documentation)",
                       )
    args = parser.parse_args()

    if args.verbose:
        print("Running in VERBOSE mode ...")
        print("Path to cofig file = \"" + str(args.path_config)+ "\"")

    # Loading config from json file
    configdata = json.load(open(args.path_config))

    L = configdata["L"]
    N = configdata["N"]
    E0 = configdata["E0"]

    # Do the one body calculation
    # TODO Group 1: Fix spelling mistakes
    OneBody = OneBodySolver(L,N)
    OneBodyEigenvalues, OneBodyEigenVectors, Info, kin, pot = OneBody.calcualteHarmonicPotential(E0)

    # Save file to hdf5 file
    OneBody.exportData(OneBodyEigenvalues, OneBodyEigenVectors, Info)

    TwoBody = TwoBodySolver()
    TwoBody.doCalculation(obDataFile='data_group1', tbDataFile='data_group2')

    # load TwoBody data
    TbData = TwoBodySpectrumData()
    TbData.open('data_group2')
    TwoBodyEigenvalues = TbData.energies[:]
    TwoBodyEigenvectors = TbData.coeficients[:,:,:]
    TbData.close()

    Transmission = TransmissionCalculator()

    current = []
    for i in range(len(configdata["muL"])):
        
        muL = configdata["muL"][i]
        muR = configdata["muR"][i]

        T = configdata["T"]
        V = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        dos = DOS()

        # TODO Group 4: There seems to be a problem with the Fermin function
        # TODO Group 4 and Group 5: Fix the interface between both modules
        #                           There is no common way to transfer the rates
        RateCal = RateCalculator()
        Gamma_L, Gamma_R = RateCal.doCalculation(OneBodyEigenvalues, TwoBodyEigenvalues,\
                                                muL, muR, T, V, \
                                                TwoBodyEigenvectors, Transmission, dos)


        ## Group 5: MASTER MODULE
        if args.verbose: print("\nSTEP 5 : MASTER EQUATION\n")
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
        # verbose only: Print all parmeters.
        if args.verbose: print("Epsilon = " + str(Epsilon))
        if args.verbose: print("n_one = " + str(n_one))
        if args.verbose: print("n_two = " + str(n_two))
        if args.verbose: print("n_tot = " + str(n_tot))
        if args.verbose: print("BCond = " + str(BCond))
        if args.verbose: print("muL = " + str(muL))
        if args.verbose: print("muR = " + str(muR))
        if args.verbose: print("Gamma_L =\n" + str(Gamma_L))
        if args.verbose: print("Gamma_R =\n" + str(Gamma_R))

        # Store input values to preserve the context.
        datapath_master = 'data_group5.hdf5'
        with h5py.File(datapath_master, "w") as f:
            f.create_dataset("Epsilon", data=Epsilon)
            f.create_dataset("n_one", data=n_one)
            f.create_dataset("n_two", data=n_two)
            f.create_dataset("n_tot", data=n_tot)
            f.create_dataset("BCond", data=BCond)
            f.create_dataset("Gamma_L", data=Gamma_L)
            f.create_dataset("Gamma_R", data=Gamma_R)

        # Calculate static or dynamic solution(s).
        mes = MasterEquationSolver()
        if BCond == "static":
            # Static solutions are requested.
            if args.verbose: print("mode = STATIC")
            stat_ps, stat_curs = mes.calculateStationarySloutions(Gamma_L, Gamma_R,
                                                                ns,
                                                                verbose=args.verbose,
                                                                )
            print("(P_stat_ij) =\n", stat_ps)
            print("(I^k_stat_i) = \n", stat_curs)
            print("muL",muL,"muR",muR,sum(stat_curs[0]))
            current.append(sum(stat_curs[0]))
            # Store output values.
            with h5py.File(datapath_master, "w") as f:
                f.create_dataset("stat_ps", data=stat_ps)
                f.create_dataset("stat_curs", data=stat_curs)
        else:
            # A simulation is requested.
            if args.verbose: print("mode = DYNAMIC")
            # Get additional parameters for simulation.
            DT = configdata["DT"]
            TMax = configdata["TMax"]
            if args.verbose: print("DT = " + str(DT))
            if args.verbose: print("TMax = " + str(TMax))
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
            if args.verbose: print("P_0 =\n" + str(P_0))
            # Store P_0 to preserve context.
            with h5py.File(datapath_master, "w") as f:
                f.create_dataset("P_0", data=P_0)
            # Calculate the simulation.
            sim_p, sim_cur = mes.simulateDynamicSloution(DT, TMax,
                                                        P_0,
                                                        Gamma_L, Gamma_R,
                                                        ns,
                                                        verbose=args.verbose,
                                                        )
            # Plot the simulations.
            # 1st plot: time development of propabilities
            sim_p.quickPlot(xlabel="t", ylabel="P",
                            xunit_symbol="ps",
                        )
            # 2nd plot: time development of netto current
            sim_cur.quickPlot(xlabel="t", ylabel="I",
                            xunit_symbol="ps",
                            )
            # Store output values.
            with h5py.File(datapath_master, "w") as f:
                f.create_dataset("sim_p", data=sim_p)
                f.create_dataset("sim_cur", data=sim_cur)

    return 0
