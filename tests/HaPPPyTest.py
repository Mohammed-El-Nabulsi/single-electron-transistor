#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import HaPPPy
import unittest

class HappyBasicTestSuite(unittest.TestCase):
    """A test class to test HaPPPy in general.

        This class is supposed to test the general behavior of HaPPPy.
        This includes test cases involving more than one 
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


    # n = 30
    # l = 50

    # energyOffset = 0.0

    # # Do the one body caclulation
    # from HaPPPy.OneBody import OneBodySolver
    # OBSolver = HaPPPy.OneBody.OneBodySolver(l, n)
    # la, v , Info = OBSolver.calcualteHarmonicPotential(energyOffset)
    # OBSolver.exportData(la, v, Info)
    
    # print(Info)

    # # # Plot the ground state for validation
    # # import matplotlib.pyplot as plt
    # # import numpy
    # # plt.plot(numpy.linspace(-l/2.0,l/2.0,n), v[:,0])
    # # plt.show() 

    # # Do the two body calculation
    # from HaPPPy.TwoBody import TwoBodySolver
    # from HaPPPy.TwoBody.TwoParticle import createTwoParticleData
    # from HaPPPy.TwoBody.OneParticleLoader import SpectrumData
    # onebodydatapath = 'data_group1'

    # import h5py
    # # obData = SpectrumData()
    # # obData.open(onebodydatapath)
    # # E, Q = createTwoParticleData(obData)

    # # twobodydatapath = 'data_group2.hdf5'
    # # dataFile = h5py.File(twobodydatapath, "w")
    # # dataSet_calcInfo = dataFile.create_dataset("E", data=E)
    # # dataSet_eigvalues = dataFile.create_dataset("Q", data=Q)  
    # # dataFile.close()

    # file = h5py.File('data_group2.hdf5', "a")
    # E = file["E"][:]
    # Q = file["Q"][:]
    # file.close()

    # muL = 0.0
    # muR = 0.0
    # T = 1 # What?
    # V = 1 # TODO what is that?
    # Rates = HaPPPy.Rates.RateCalculator()
    # Gamma_01L,Gamma_01R,Gamma_10L,Gamma_10R,Gamma_12L,Gamma_12R,Gamma_21L,Gamma_21R = Rates.doCalculation(la, E, muL, muR, T, V, Q)

    # # instantiate a MasterEquationSolver
    # mes = HaPPPy.MasterEquation.MasterEquationSolver()

    # TMax = 

    # # parameters of this tets:
    # ε = mes.getε()
    # a = 1 # > 0
    # p0 = 0.7 # > 0
    # p1 = 0.1 # > 0
    # p2 = 0.2 # > 0, p0 + p1 + p2 = 1
    # Δt = 1
    # t_max = 100
    # # set-up a reasonable Γ-matrix
    # Γ_L = np.array([[0, a, 0], [0, 0, 0], [0, 0, 0]])
    # Γ_R = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    # # choose a legitimate start value for P_0 (P_0 = P(t=0))
    # P_0 = np.array([p0, p1, p2])
    # # simulate
    # sim_time_dev_prop, sim_current = = mes.doCalculation(Δt, t_max, P_0, Γ_L, Γ_R)

    # ## plot results
    # # 1st plot: time development of propabilities
    # sim_time_dev_prop.quickPlot(xlabel="t", ylabel="P")
    # # 2nd plot: time development of netto current
    # sim_current.quickPlot(xlabel="t", ylabel="I")

