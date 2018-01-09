
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
import math


__docformat__ = 'reStructuredText'

class MasterEquationSolver:

    def doCalculation(self,
                      P_0,
                      Γ_L,
                      Γ_R,
                      t_max,
                      Δt,
                      ε=1E-10,
                      check_tolerance=True
                     ):
        """Simulates the time development of propabilities P(t) for an ``n`` state system with transition rates Γ including the netto current.

        .. todo::

            - Missing e in current simulation?
            - Extra return type for each simulation recommended.

        This method returns numerical approximations of the solution to the following problems.

        **1st Problem: Time Development of Propabilities**

        Let

        .. math::

            n \in \mathbb{N}, t_{max} \in \mathbb{R}^+_0, \Delta t \in \mathbb{R}^+ \\\\ \Gamma = \Gamma^L + \Gamma^R \in Mat(n, n, \mathbb{R}^+_0), \\vec{P} : [0, t_{max}] \\to { \left ( \mathbb{R}^+_0 \\right ) }^n, t \mapsto \\vec{P}(t)

        The differential equation of first order with constant coefficients to be solved is stated as

        .. math::

            \\frac{d P_{\\alpha}}{d t} = \sum_{\\beta} \Gamma_{\\beta \\rightarrow \\alpha} P_{\\beta} - \sum_{\\beta} \Gamma_{\\alpha \\rightarrow \\beta} P_{\\alpha}

        where

        .. math::

            \Gamma_{\\alpha \\rightarrow \\beta} \equiv \Gamma_{\\alpha \\beta}

        and

        .. math::

            \\vec{P_0} = \\vec{P}(t=0)

        is the boundary condition.

        denotes the rate from state :math:`\\alpha` to :math:`\\beta`.

        The solution of this problem is returned as a discrete version of :math:`\\vec{P}`. Therefore it consists of pairs :math:`(t_i, \\vec{P}(t_i))` where :math:`t_i \in \{n \cdot \Delta t | n \in \mathbb{N}_0 \land n \cdot \Delta t \leq t_{max} \}`.

        **2nd Problem: Netto Current**

        From the time develpment of probabilities the netto current is derived.

        In addition to the 1st problem let

        .. math::

            I : [0, t_{max}] \\to \mathbb{R}, t \mapsto I(t) \\\\
            I^k : [0, t_{max}] \\to \mathbb{R}, t \mapsto I^k(t), k \in \\{ L, R \\} \\\\
            I^k = e \sum_{\\alpha, \\beta} \Gamma^k_{\\alpha \\rightarrow \\beta} P_{\\alpha} \\\\
            I = I^L - I^R = e \sum_{\\alpha, \\beta} \\left ( \Gamma^L_{\\alpha \\rightarrow \\beta} - \Gamma^R_{\\alpha \\rightarrow \\beta} \\right ) P_{\\alpha}

        With :math:`e \\approx 1.602 \cdot 10^{-19} Coulomb`.

        :math:`I^L,I^R` descibe the current the the *left* resp. *right* barrier, hence :math:`I` describes the netto flow.

        The solution of this problem is returned as a discrete version as well. It consists of pairs :math:`(t_i, I(t_i))` where :math:`t_i \in \{n \cdot \Delta t | n \in \mathbb{N}_0 \land n \cdot \Delta t \leq t_{max} \}`. (Just like in the first problem.)

        :param P_0: Start value of propabilities where :code:`P_0` ≣ :math:`\\vec{P_0}`. Must be either a list or a ``nx1`` matrix.
        :type P_0: numpy.array
        :param Γ_L: Matrix containing the transition rates regarding electrons going to the *left* where :code:`Γ_L[i][j]` ≣ :math:`\Gamma^L_{i \\rightarrow j}`. Must be a ``nxn`` matrix.
        :type Γ_L: numpy.array
        :param Γ_R: Matrix containing the transition rates regarding electrons going to the *right* where :code:`Γ_R[i][j]` ≣ :math:`\Gamma^R_{i \\rightarrow j}`. Must be a ``nxn`` matrix.
        :type Γ_R: numpy.array
        :param t_max: Last point in time to be simulated. Must be >= 0.
        :type t_max: float
        :param Δt: Length of the time tintervall between two simulated events. Must be > 0.
        :type Δt: float
        :param ε: Tolerance value. Sum of all probabilities must be equals to 1 within this tolerance.
        :type ε: float
        :param check_tolerance: Enables initial and successive check as described for ε.
        :type check_tolerance: bool

        :return: Returns :code:`sim_time_dev_prop, sim_current` where :code:`sim_time_dev_prop[i][0]` ≣ :math:`i \cdot \Delta t`, :code:`sim_time_dev_prop[i][1][j]` ≣ :math:`P_j(i \cdot \Delta t)`, :code:`sim_current[i][0]` ≣ :math:`i \cdot \Delta t` and :code:`sim_current[i][1]` ≣ :math:`I(i \cdot \Delta t)`.
        :rtype: (list, list)

        :example: .. code-block:: python
                :emphasize-lines: 1,3,7-8,10-18

                import numpy as np
                import matplotlib.pyplot as plt
                from HaPPPy.MasterEquation.MasterEquationSolver import MasterEquationSolver as MES

                ## simple simulation with the MasterEquationSolver class
                # set-up a reasonable Γ-matrix
                Γ_L = np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
                Γ_R = np.array([[0, 0.0, 0.5], [0.0, 0, 0.5], [0.0, 0.5, 0]])
                # choose a legitimate start value for P_0 (P_0 = P(t=0))
                P_0 = np.array([0.9, 0.0, 0.1])
                # simulate
                mes = MES()
                sim_time_dev_prop, sim_current = mes.doCalculation(P_0,
                                                                   Γ_L,
                                                                   Γ_R,
                                                                   t_max=10,
                                                                   Δt=0.1
                                                                  )

                ## plot results
                # 1st plot: time development of propabilities
                # seperate x value (time) from y values (P) to be able to plot it
                ts = [sim_time_dev_prop[i][0] for i in range(len(sim_time_dev_prop))]
                Ps = [sim_time_dev_prop[i][1] for i in range(len(sim_time_dev_prop))]
                # plot
                plt.plot(ts, Ps)
                plt.xlabel("t")
                plt.ylabel("P")# create legend
                n = len(sim_time_dev_prop[0][1]) # dimension of P
                legend_names = ["$P_" + str(i) + "$" for i in range(n)]
                plt.legend(legend_names)
                plt.grid()
                plt.show()

                # 2nd plot: time development of netto current
                # seperate x value (time) from y values (P) to be able to plot it
                ts = [sim_current[i][0] for i in range(len(sim_current))]
                Is = [sim_current[i][1] for i in range(len(sim_current))]
                # plot
                plt.plot(ts,Is)
                plt.xlabel("t")
                plt.ylabel("I")
                # create legend
                legend_names = ["$I$"]
                plt.legend(legend_names)
                plt.grid()
                plt.show()
                plt.show()


            The relevant lines of code for the simulation to work are highlighted. To give a real live example and to demonstarte the usage of the result some code to plot the result is added.

        """
        # Don't remove the last dot/paragrph since it is a workaround to shpinx's
        # code-block interpreter.


        ## type conversions
        # if necessary: reformat P_0 as matrix (otherwise P_0 could not be multiplicated with matricies)
        if P_0.ndim == 1:
            P_0 = np.array([P_0]).transpose()

        ## input checks
        # warn if tolerance value is unreasonable
        if check_tolerance and ε <= 0:
            raise RuntimeError("ε must be a positive number > 0. \nε = " + str(ε))
        # P_0 must be a nx1 matrix (with n matching Λ)
        if P_0.ndim != 2 or P_0.shape[1] != 1:
            raise RuntimeError("P_0 must be a "
                               + str(Γ.shape[0])
                               + "x1 matrix (aka. a 'dotable vector')! \nP_0 = \n"
                               + str(P_0)
                              )
        # Γ_L must be a nxn matrix
        if Γ_L.ndim != 2 or Γ_L.shape[0] != P_0.shape[0] or Γ_L.shape[1] != P_0.shape[0] or not (Γ_L >= 0).all():
            raise RuntimeError("Γ_L must be a nxn matrix with coefficients >= 0. \nΓ_L = \n" + str(Γ_L))
        # Γ_R must be a nxn matrix
        if Γ_R.ndim != 2 or Γ_R.shape[0] != P_0.shape[0] or Γ_R.shape[1] != P_0.shape[0] or not (Γ_R >= 0).all():
            raise RuntimeError("Γ_R must be a nxn matrix with coefficients >= 0. \nΓ_R = \n" + str(Γ_R))
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
        if t_max < 0 or Δt <= 0:
            raise RuntimeError("Simulated time intervals must be finite and positive ("
                               + "t_max >= 0 and Δt > 0).\n"
                               + "t_max = " + str(t_max)
                               + ", Δt = " + str(Δt)
                              )

        ## calculation
        Γ = Γ_L + Γ_R
        sim_successful, sim = self.__simulate_time_development_of_propabilities(P_0, Γ, t_max, Δt, ε, check_tolerance)
        if sim_successful:
            I_t = self.__simulate_current(Γ_L, Γ_R, sim, ε)

        return sim, I_t

    def __simulate_time_development_of_propabilities(self,
                                                   P_0,
                                                   Γ,
                                                   t_max,
                                                   Δt,
                                                   ε=1E-10,
                                                   check_tolerance=True,
                                                   verbose=False
                                                  ):
        """Simulates the time development of propabilities P(t) for an :math:`n`-state system with transition rates Γ.

        .. todo::

            Create extra return type for simulation.

        This is a private function and should not be called from outside the MasterEquationSolver class! (Hence it does not perform further tests than MasterEquationSolver.)

        :return: Returns :code:`sim_sucessfull, sim_time_dev_prop`. Iff :code:`sim_sucessfull==True` than :code:`sim_time_dev_prop` has values like describet in MasterEquationSolver.doCalculation (see return) otherwise sim_time_dev_prop is undefined.
        :rtype: (bool, list)

        See documentation of MasterEquationSolver.doCalculation (1st problem) for more details.

        """

        # print input values (if verbose)
        if verbose:
            print("P(t=0) = \n", P_0)
            print("Γ = \n", Γ)
            if check_tolerance:
                print("ε = ", ε)

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
        # simulate time development discretly
        sim = []
        P_0_evb = np.dot(Λ_evecs_T, P_0)
        for t in np.arange(0, t_max + Δt, Δt):
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

    def __simulate_current(self, Γ_L, Γ_R, sim_time_dev_prop, ε=1E-10, verbose=False):
        """Simulates the time current I(t) for an :math:`n`-state system with transition rates Γ.

        .. todo::

            Create extra return type for simulation.

        This is a private function and should not be called from outside the MasterEquationSolver class! (Hence it does not perform further tests than MasterEquationSolver.)

        :return: Returns :code:`sim_current`. :code:`sim_time_dev_prop` needs to be valid for :code:`sim_current` to be valid!
        :rtype: (bool, list)

        See documentation of MasterEquationSolver.doCalculation (2nd problem) for more details.

        """

        # calculate ΔΓ (direction sensivtive form)
        ΔΓ = Γ_L - Γ_R
        if verbose: print("Γ: \n", Γ)

        # simulate the current discretly
        sim = []
        for i in range(len(sim_time_dev_prop)):
            sim.append((sim_time_dev_prop[i][0], sum(np.dot(ΔΓ, sim_time_dev_prop[i][1]))))
            if verbose: print("I_t[" + str(i) + "= \n", sim[i])

        return sim
