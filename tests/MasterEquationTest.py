
import HaPPPy
import unittest
import numpy as np

class MasterEquationTestSuite(unittest.TestCase):
    """A test class to the MasterEquation module.

    """

    def test_MasterEquation_exists(self):
        """ Checks wether the MasterEquation module exists
        """
        self.assertTrue(hasattr(HaPPPy, 'MasterEquation'))

    def test_MasterEquation_simulateDynamicSloution(self,
                                                    plot_figures=True,
                                                    verbose=False,
                                                   ):
        """
        This method tests the bahaviour of the HaPPPy.MasterEquation module
        by checking some simple edge cases.
        """

        do_test_1 = True
        do_test_2 = True
        do_test_3 = True
        do_test_4 = True
        do_test_5 = True

        if do_test_1:
            ## test 1:

            # instanciate a MasterEquationSolver
            mes = HaPPPy.MasterEquation.MasterEquationSolver()
            # parameters of this tets:
            ε = mes.getε()
            δ = 1E-2
            a = 1 # > 0
            p0 = 0.1 # > 0
            p1 = 0.7 # > 0
            p2 = 0.2 # > 0, p0 + p1 + p2 = 1
            Δt = 1
            t_max = 100
            # set-up a reasonable Γ-matrix
            Γ_L = np.array([[0, a, 0], [0, 0, 0], [0, 0, 0]])
            Γ_R = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
            # choose a legitimate start value for P_0 (P_0 = P(t=0))
            P_0 = np.array([p0, p1, p2])
            # simulate
            sim_tdp, sim_cur = mes.simulateDynamicSloution(
                                   Δt, t_max, P_0, Γ_L, Γ_R, verbose=verbose
                                  )
            stat_ps, stat_curs = mes.calculateStationarySloutions(
                                    Γ_L, Γ_R, verbose=verbose
                                   )
            # plot result
            if plot_figures:
                sim_tdp.quickPlot(xlabel="t", ylabel="P")
                sim_cur.quickPlot(xlabel="t", ylabel="I")
            # check validity
            self.assertTrue(sim_tdp.valid())
            self.assertTrue(sim_cur.valid())
            # check development of porpabilities
            Ps = sim_tdp.getValues()
            #     a. check conservation of total probability
            for P in Ps:
                self.assertTrue(1 - ε <= sum(P) <= 1 + ε)
            #     b. check values after a 'long' time period
            #        expected: propability of state 1 'shifts' to state 2
            #                  state 3 is constant
            self.assertTrue(abs(Ps[t_max][0] - (p0 + p1)) <= δ)
            self.assertTrue(abs(Ps[t_max][1])  <= δ)
            self.assertTrue(abs(Ps[t_max][2] - p2) <= ε)

        if do_test_2:
            ## test 2:
            # symmetric Γ = Γ_L - Γ_R --> uniform distibution of propability

            # instanciate a MasterEquationSolver
            mes = HaPPPy.MasterEquation.MasterEquationSolver()
            # parameters of this tets:
            ε = mes.getε()
            Δt = 1
            t_max = 100
            n = 50
            # set-up a reasonable Γ-matrix
            Γ_L = np.array([[(i + j)/n**2 for j in range(n)] for i in range(n)])
            Γ_R = 2 * Γ_L
            # choose a legitimate start value for P_0 (P_0 = P(t=0))
            P_0 = np.array([i for i in range(n)])
            P_0 = P_0 / sum(P_0)
            # simulate
            sim_tdp, sim_cur = mes.simulateDynamicSloution(
                                   Δt, t_max, P_0, Γ_L, Γ_R, verbose=verbose
                                  )
            stat_ps, stat_curs = mes.calculateStationarySloutions(
                                    Γ_L, Γ_R, verbose=verbose
                                   )
            # plot result
            if plot_figures:
                sim_tdp.quickPlot(xlabel="t", ylabel="P")
                sim_cur.quickPlot(xlabel="t", ylabel="I")
            # check validity
            self.assertTrue(sim_tdp.valid())
            self.assertTrue(sim_cur.valid())
            # check development of porpabilities
            Ps = sim_tdp.getValues()
            #     a. check conservation of total probability
            for P in Ps:
                self.assertTrue(1 - ε <= sum(P) <= 1 + ε)
            #     b. check values after a 'long' time period
            #        expected: uniform distribution
            for P_i in Ps[t_max]:
                self.assertTrue(1/n - ε <= P_i <= 1/n + ε)

        if do_test_3:
            ## test 3:

            # instanciate a MasterEquationSolver
            mes = HaPPPy.MasterEquation.MasterEquationSolver()
            # parameters of this tets:
            ε = mes.getε()
            Δt = 1
            t_max = 100
            n = 50
            # set-up a reasonable Γ-matrix
            Γ_L = np.zeros((n, n))
            Γ_R = np.zeros((n, n))
            # choose a legitimate start value for P_0 (P_0 = P(t=0))
            P_0 = np.array([i for i in range(n)])
            P_0 = P_0 / sum(P_0)
            # simulate
            sim_tdp, sim_cur = mes.simulateDynamicSloution(
                                   Δt, t_max, P_0, Γ_L, Γ_R, verbose=verbose
                                  )
            stat_ps, stat_curs = mes.calculateStationarySloutions(
                                    Γ_L, Γ_R, verbose=verbose
                                   )
            # plot result
            if plot_figures:
                sim_tdp.quickPlot(xlabel="t", ylabel="P")
                sim_cur.quickPlot(xlabel="t", ylabel="I")
            # check validity
            self.assertTrue(sim_tdp.valid())
            self.assertTrue(sim_cur.valid())
            # check development of porpabilities
            Ps = sim_tdp.getValues()
            #     a. check conservation of individual probabilities
            for P in Ps:
                self.assertTrue(-ε <= sum(abs(P - P_0)) <= ε)
            # check current
            Is = sim_cur.getValues()
            #     a. check if current is constntly 0
            for I in Is:
                self.assertTrue(-ε <= abs(I) <= +ε)

        if do_test_4:
            ## test 4:
            # test if algorith can handle large inputs
            # (like test 2 but with large n to test accuracy)

            # instanciate a MasterEquationSolver
            mes = HaPPPy.MasterEquation.MasterEquationSolver()
            # parameters of this tets:
            ε = mes.getε()
            Δt = 1
            t_max = 100
            n = 100
            # set-up a reasonable Γ-matrix
            Γ_L = np.array([[(i + j)/n**2 for j in range(n)] for i in range(n)])
            Γ_R = 2 * Γ_L
            # choose a legitimate start value for P_0 (P_0 = P(t=0))
            P_0 = np.array([i for i in range(n)])
            P_0 = P_0 / sum(P_0)
            # simulate
            sim_tdp, sim_cur = mes.simulateDynamicSloution(
                                   Δt, t_max, P_0, Γ_L, Γ_R, verbose=verbose
                                  )
            stat_ps, stat_curs = mes.calculateStationarySloutions(
                                    Γ_L, Γ_R, verbose=verbose
                                   )
            # check validity
            self.assertTrue(sim_tdp.valid())
            self.assertTrue(sim_cur.valid())
            # check development of porpabilities
            Ps = sim_tdp.getValues()
            #     a. check conservation of total probability
            for P in Ps:
                self.assertTrue(1 - ε <= sum(P) <= 1 + ε)
            #     b. check values after a 'long' time period
            #        expected: uniform distribution
            for P_i in Ps[t_max]:
                self.assertTrue(1/n - ε <= P_i <= 1/n + ε)

        if do_test_5:
            ## test 5:
            # test tolerance behaviour

            # parameters of this tets:
            ε = 1E-100 # ridiculously precise
            Δt = 1
            t_max = 100
            n = 10
            # set-up a reasonable Γ-matrix
            Γ_L = np.array([[(i + j)/n**2 for j in range(n)] for i in range(n)])
            Γ_R = 2 * Γ_L
            # choose a legitimate start value for P_0 (P_0 = P(t=0))
            P_0 = np.array([i for i in range(n)])
            P_0 = P_0 / sum(P_0)
            # instanciate a MasterEquationSolver
            mes = HaPPPy.MasterEquation.MasterEquationSolver(ε=ε)
            # simulate
            sim_tdp, sim_cur = mes.simulateDynamicSloution(
                                   Δt, t_max, P_0, Γ_L, Γ_R, verbose=verbose
                                  )
            stat_ps, stat_curs = mes.calculateStationarySloutions(
                                    Γ_L, Γ_R, verbose=verbose
                                   )
            # plot result
            if plot_figures:
                sim_tdp.quickPlot(xlabel="t", ylabel="P")
                sim_cur.quickPlot(xlabel="t", ylabel="I")
            # check validity
            self.assertTrue(not sim_tdp.valid())
            self.assertTrue(not sim_cur.valid())

if __name__ == '__main__':
    master_equation_suite = unittest.TestLoader().loadTestsFromTestCase(MasterEquationTestSuite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(master_equation_suite)
