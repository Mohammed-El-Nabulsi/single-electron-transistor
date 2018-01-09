
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
import math


__docformat__ = 'reStructuredText'

class MasterEquationSolver:
    """ Solves the master equation

    TODO: Add more documentation

    """


    def __init__(self):
        """ The constructor.
        """

        print("Hello from the MasterEquationSolver, what can I do for you?")

    def doCalculation(self, Γ_L, Γ_R, P_0, t_max, t_delta, ε=1E-10, check_tolerance=True):
        Γ = Γ_L + Γ_R
        sim_successful, sim = self.simulate_time_development_of_propabilities(P_0, Γ, t_max, t_delta, ε, check_tolerance)
        if sim_successful:
            I_t = self.calculate_current(Γ_L, Γ_R, sim, ε)

        return I_t, sim

    def simulate_time_development_of_propabilities(self,
                                                   P_0,
                                                   Γ,
                                                   t_max,
                                                   t_delta,
                                                   ε=1E-10,
                                                   check_tolerance=True,
                                                   verbose=False
                                                  ):
        """Simulates the time development of propabilities P(t) for an ``n`` state system with transition rates Γ.

        .. todo::

            Create extra return type for simulation.

        This method returns a numerical approximation of the solution to the following problem.

        Let

        .. math::

            n \in \mathbb{N}, t_{max} \in \mathbb{R}^+_0, t_{delta} \in \mathbb{R}^+ \\\\ \Γ \in Mat(n, n, \mathbb{R}^+_0), \\vec{P} : [0, t_{max}] \\to { \left ( \mathbb{R}^+_0 \\right ) }^n, t \mapsto \\vec{P}(t)

        The differential equation of first order with constant coefficients to be solved is stated as

        .. math::

            \\frac{d P_{\\alpha}}{d t} = \sum_{\\beta} \Γ_{\\beta \\rightarrow \\alpha} P_{\\beta} - \sum_{\\beta} \Γ_{\\alpha \\rightarrow \\beta} P_{\\alpha}

        where

        .. math::

            \Γ_{\\alpha \\rightarrow \\beta} \equiv \Γ_{\\alpha \\beta}

        denotes the rate from state :math:`\\alpha` to :math:`\\beta`.

        The returned solution for :math:`\\vec{P}` is discrete and is therefore consists of pairs :math:`(t_i, \\vec{P}(t_i))` where :math:`t_i \in \{n \cdot t_{delta} | n \in \mathbb{Z} \land n \cdot t_{delta} \leq t_{max} \}`

        :param P_0: Start value of propabilities. Must be either a list or a ``nx1`` matrix.
        :type P_0: numpy.array
        :param Γ: Matrix containing the transition rates where :code:`Γ[i][j]` ≣ :math:`\Γ_{i \\rightarrow j}`. Must be a ``nxn`` matrix.
        :type Γ: numpy.array
        :param t_max: Last point in time to be simulated. Must be >= 0.
        :type t_max: float
        :param t_delta: Length of the time tintervall between two simulated events. Must be > 0.
        :type t_delta: float
        :param ε: Tolerance value. Sum of all probabilities must be equals to 1 within this tolerance.
        :type ε: float
        :param check_tolerance: Enables initial and successive check as described for ε.
        :type check_tolerance: bool
        :param verbose: Enables the verbose mode: Calculation will be printed in detail to support debugging.
        :type verbose: bool

        :return: Returns a pair. The first value is called ``sim_successful`` and is ``True`` iff no issue occurred during the simulation. The second value is a list of pairs representing the result. The first values track te time and the second contains the probabilities as vector given as a list of numbers.
        :rtype: (bool, list)

        :example: .. code-block:: python
                :emphasize-lines: 1,3,6-16

                import numpy as np
                import matplotlib.pyplot as plt
                from HaPPPy.MasterEquation.MasterEquationSolver import MasterEquationSolver as MES

                ## test program for simulate_time_development_of_propabilities
                # set-up a reasonable Γ-matrix
                Γ = np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
                # choose a legitimate start value for P_0 (P_0 = P(t=0))
                P_0 = np.array([0.9, 0.1, 0])
                # simulate
                sim_successful, sim = MES.simulate_time_development_of_propabilities(
                                   P_0,
                                   Γ,
                                   t_max=20,
                                   t_delta=0.5
                                  )

                ## plot result
                # seperate x value (time) from y values (P) to be able to plot it
                ts = [sim[i][0] for i in range(len(sim))]
                Ps = [sim[i][1] for i in range(len(sim))]
                n = len(sim[0][1]) # dimension of P
                legend_names = ["$P_" + str(i) + "$" for i in range(n)]
                plt.plot(ts, Ps)
                plt.xlabel("t")
                plt.ylabel("P")
                plt.legend(legend_names)
                plt.grid()
                if not sim_successful:
                     plt.text(0, 0, "Simulation failed!")
                plt.show()

            The relevant lines of code for the simulation to work are highlighted. To give a real live example and to demonstarte the usage of the result some code to plot the result is added.

        """
        # Don't remove the last dot since it is a workaround to shpinx's
        # code-block interpreter.

        # if necessary: reformat P_0 as matrix (otherwise P_0 could not be multiplicated with matricies)
        if P_0.ndim == 1:
            P_0 = np.array([P_0]).transpose()

        if verbose:
            print("P(t=0) = \n", P_0)
            print("Γ = \n", Γ)
            if check_tolerance:
                print("ε = ", ε)

        ## input checks
        # warn if tolerance value is unreasonable
        if check_tolerance and ε <= 0:
            raise RuntimeError("ε must be a positive number > 0. \nε = " + str(ε))
        # Γ must be a nxn matrix
        if Γ.ndim != 2 or Γ.shape[0] != Γ.shape[1] or not (Γ >= 0).all():
            raise RuntimeError("Γ must be a square matrix with coefficients >= 0. \nΓ = \n" + str(Γ))
        # P_0 must be a nx1 matrix (with n matching Λ)
        if P_0.ndim != 2 or P_0.shape[0] != Γ.shape[0] or P_0.shape[1] != 1:
            raise RuntimeError("P_0 must be a "
                               + str(Γ.shape[0])
                               + "x1 matrix (aka. a 'dotable vector')! \nP_0 = \n"
                               + str(P_0)
                              )
        # P_0 must have coefficients >= 0
        if not (P_0 >= 0).all():
            raise RuntimeError("P_0 must have coefficients >= 0. \nP_0 = \n" + str(P_0))
        # coefficients of P_0 must add up to 1
        P_sum = sum(P_0)
        if check_tolerance and (P_sum < 1 - ε or P_sum > 1 + ε):
            raise RuntimeError("Coefficients of P_0 must add up to 1 (within tolerance ε = "
                               + str(ε)
                               + "). \nP_0 = \n"
                               + str(P_0)
                               + "\n Σ = "
                               + str(P_sum)
                              )
        # simulated time intervals must be positive or negative
        if t_max < 0 or t_delta <= 0:
            raise RuntimeError("Simulated time intervals must be finite and positive ("
                               + "t_max >= 0 and t_delta > 0).\n"
                               + "t_max = " + str(t_max)
                               + ", t_delta = " + str(t_delta)
                              )


        ## simulation
        # track if simulation had any issues
        sim_successful = True

        # set-up Λ-matrix
        Λ_in = Γ
        Λ_out = np.diag([sum([Γ[i][j] for j in range(Γ.shape[0])]) for i in range(Γ.shape[0])])
        Λ = Λ_in - Λ_out
        if verbose: print("Λ_in = \n", Λ_in)
        if verbose: print("Λ_out = \n", Λ_out)
        if verbose: print("Λ = \n", Λ)

        # find eigenvector basis of Λ (Λ_evecs_T transforms to eigenvevtor basis of Λ)
        (Λ_evals, Λ_evecs) = lin.eig(Λ)
        Λ_evecs_T = Λ_evecs.transpose()
        if verbose: print("Λ.eigenvalues = \n", Λ_evals)
        if verbose: print("Λ.eigenvectors = \n", Λ_evecs_T)

        # get P_0 in eigenvector basis of Λ
        P_evb = Λ_evecs_T.dot(P_0)
        if verbose: print("P(t=0) in Λ.eigenvectorbase = \n", P_evb)

        # get Λ eigenvector basis of Λ
        Λ_evb = np.diag(Λ_evals)
        if verbose: print("Λ in Λ.eigenvectorbase = \n", Λ_evb)

        # get inverse of Λ_evecs_T (transforms back from eigenvevtor basis of Λ)
        Λ_evecs_T_inv = np.linalg.inv(Λ_evecs_T)
        if verbose: print("Λ.eingenvectors^-1 = \n", Λ_evecs_T_inv)

        # check if creation inversion was successful
        if (np.dot(Λ_evecs_T_inv, Λ_evecs_T) != np.identity(Λ.shape[0])).all():
            raise RuntimeError("Can not invert Λ = \n" + str(Λ))

        verbose = False
        # time development
        sim = []
        P_0_evb = np.dot(Λ_evecs_T, P_0)
        for t in np.arange(0, t_max + t_delta, t_delta):
            if verbose: print("\nt = ", t)
            #exp_tΛ_evb = np.diag(np.power(Λ_evals, t)) # Does not behave well due to numerical issues!
            exp_tΛ_evb = np.diag(np.exp(t * Λ_evals))
            if verbose: print("exp(" + str(t) + " * Λ) in Λ.eigenvectorbase = \n", exp_tΛ_evb)
            P_evb_t = np.dot(exp_tΛ_evb, P_0_evb)
            if verbose: print("P(t=" + str(t) + ") in Λ.eigenvectorbase = \n", P_evb_t)
            P_t = np.dot(Λ_evecs_T_inv, P_evb_t)
            if verbose: print("P(t=" + str(t) + ") = \n", P_t)
            sim.append((t, [P_t[i][0] for i in range(P_t.shape[0])]))
            # check if sum of coefficients of P_t is still 1
            P_sum = sum(P_t)
            if check_tolerance and (P_sum < 1 - ε or P_sum > 1 + ε):
                print("Warning! Calculation aborted Coefficients of P_(t="
                      + str(t)
                      + ") must add up to 1 (within tolerance ε = "
                      + str(ε)
                      + "). \nP_(t="
                      + str(t)
                      + ") = \n"
                      + str(P_t)
                      + "\n Σ = "
                      + str(P_sum[0])
                     )
                sim_successful = False
                break

        # Use either:
        #   return sim_successful ? sim : None
        # or
        #   return sim_successful, sim
        # to return the simulated data. The last option allowes to analyse the
        # simulation up to the point of failure but must be analysed more carefully.
        return sim_successful, sim

    def calculate_current(self, Γ_L, Γ_R, sim, ε=1E-10, verbose=False):
        #Γ_L has to be an nxn Matrix
        
        #Γ_R has to be an nxn Matrix

        I_t = []
        Γ = Γ_L - Γ_R
        if verbose: print("Γ: \n", Γ)

        for i in range(len(sim)):
            I_t.append((sim[i][0], sum(np.dot(Γ, sim[i][1]))))
            if verbose: print("I_t[" + str(i) + "= \n", I_t[i])

        return I_t

#### test program for simulate_time_development_of_propabilities
### set-up a reasonable Γ-matrix
##Γ_l = np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
##Γ_r = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
### choose a legitimate start value for P_0 (P_0 = P(t=0))
##P_0 = np.array([0.9, 0.1, 0])
### simulate
##MES = MasterEquationSolver()
##I_t, sim = MES.doCalculation(
##                   Γ_l,
##                   Γ_r,
##                   P_0,
##                   check_tolerance = False,
##                   t_max=20000000000,
##                   t_delta=1000000
##                  )
##
#### plot result
### seperate x value (time) from y values (P) to be able to plot it
##ts = [sim[i][0] for i in range(len(sim))]
##Ps = [sim[i][1] for i in range(len(sim))]
##n = len(sim[0][1]) # dimension of P
##legend_names = ["$P_" + str(i) + "$" for i in range(n)]
##plt.plot(ts, Ps)
##plt.xlabel("t")
##plt.ylabel("P")
##plt.legend(legend_names)
##plt.grid()
##plt.show()
##
##ts = [I_t[i][0] for i in range(len(sim))]
##Ps = [I_t[i][1] for i in range(len(sim))]
##plt.plot(ts, Ps)
##plt.xlabel("t")
##plt.ylabel("P")
##plt.legend(legend_names)
##plt.grid()
##plt.show()

