
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
import math


__docformat__ = 'reStructuredText'

class MasterEquationSolver:
    """
    Simulates the time development of propabilities :math:`\\vec{P}(t)` for an :math:`n`-state system with transition rates :math:`\\Gamma` including the netto current :math:`I(t)`.

    This method returns numerical approximations of the solutions to the following problems.

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

    denotes the rate from state :math:`\\alpha` to :math:`\\beta` and

    .. math::

        \\vec{P_0} = \\vec{P}(t=0)

    is the boundary condition..

    The solution of this problem is returned as a discrete simulation of :math:`\\vec{P}`. (See HaPPPy.MasterEquation.Simulation for more details.)

    **2nd Problem: Netto Current**

    From the time develpment of probabilities the netto current is derived.

    In addition to the 1st problem let

    .. math::

        I : [0, t_{max}] \\to \mathbb{R}, t \mapsto I(t) \\\\
        I^k : [0, t_{max}] \\to \mathbb{R}, t \mapsto I^k(t), k \in \\{ L, R \\} \\\\
        I^k = q \sum_{\\alpha, \\beta} \Gamma^k_{\\alpha \\rightarrow \\beta} P_{\\alpha} \\\\
        I = I^L - I^R = q \sum_{\\alpha, \\beta} \\left ( \Gamma^L_{\\alpha \\rightarrow \\beta} - \Gamma^R_{\\alpha \\rightarrow \\beta} \\right ) P_{\\alpha}

    It is assumed that the charge per particle :math:`q = 1`.

    :math:`I^L,I^R` descibe the current through the *left* resp. *right* barrier, hence :math:`I` describes the netto flow from left to right.

    The solution of this problem is returned as a discrete simulation as well. (See HaPPPy.MasterEquation.Simulation for more details.)

    :example: .. code-block:: Python
            :emphasize-lines: 1-2,7-8,10, 12-13

            import HaPPPy
            import numpy as np
            import matplotlib.pyplot as plt

            ## simple simulation with the MasterEquationSolver class
            # set-up a reasonable Γ-matrix
            Γ_L = np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
            Γ_R = np.array([[0, 0.0, 0.5], [0.0, 0, 0.5], [0.0, 0.5, 0]])
            # choose a legitimate start value for P_0 (P_0 = P(t=0))
            P_0 = np.array([0.9, 0.0, 0.1])
            # simulate
            mes = HaPPPy.MasterEquation.MasterEquationSolver()
            sim_time_dev_prop, sim_current = mes.doCalculation(0.1, 10, P_0, Γ_L, Γ_R)

            ## plot results
            # 1st plot: time development of propabilities
            sim_time_dev_prop.quickPlot(xlabel="t", ylabel="P")
            # 2nd plot: time development of netto current
            sim_current.quickPlot(xlabel="t", ylabel="I")


        The relevant lines of code for the simulation to work are highlighted. To give a real live example and to demonstarte the usage of the result some code to plot the result is added.

    """
    # Don't remove the last dot/paragrph since it is a workaround to shpinx's
    # code-block interpreter.

    def __init__(self, ε=None):
        """
        :param ε: Tolerance value. Sum of all probabilities must be equals to 1 within this tolerance. If :code:`ε == None` the fedault value is used (see getDefaultε()).
        :type ε: float
        """
        if ε == None:
            self.setε(self.getDefaultε())
        else:
            self.ε = ε

    def getε(self):
        """
        :return: Returns ε - a tolerance value: The sum of all probabilities must be equals to 1 within this tolerance.
        :rtype: float
        """
        return self.ε

    def setε(self, ε):
        """
        :param ε: Tolerance value: The sum of all probabilities must be equals to 1 within this tolerance.
        :type ε: float
        """
        self.ε = ε

    @staticmethod
    def getDefaultε():
        """
        :return: Returns the default value of ε - a tolerance value: The sum of all probabilities must be equals to 1 within this tolerance.
        :rtype: float
        """
        return 1E-10

    def doCalculation(self, Δt, t_max, P_0, Γ_L, Γ_R, check_tolerance=True):
        """
        See class description of Happy.MasterEquationSolver for details.

        :param P_0: Start value of propabilities where :code:`P_0` ≣ :math:`\\vec{P_0}`. Must be either a list or a ``nx1`` matrix.
        :type P_0: numpy.array or numpy.ndarray
        :param Γ_L: Matrix containing the transition rates regarding electrons tunneling trough the *left* barrier where :code:`Γ_L[i][j]` ≣ :math:`\Gamma^L_{i \\rightarrow j}`. Must be a ``nxn`` matrix.
        :type Γ_L: numpy.ndarray
        :param Γ_R: Matrix containing the transition rates regarding electrons tunneling trough the *right* barrier where :code:`Γ_R[i][j]` ≣ :math:`\Gamma^R_{i \\rightarrow j}`. Must be a ``nxn`` matrix.
        :type Γ_R: numpy.ndarray
        :param t_max: Last point in time to be simulated. Must fulfil :code:`t_max >= 0`.
        :type t_max: float
        :param Δt: Length of the time tintervall between two simulated events. Must fulfil :code:`Δt > 0`.
        :type Δt: float

        :return: Returns :code:`sim_time_dev_prop, sim_current` where both values are a HaPPPy.MasterEquation.Simulation of :math:`\\vec{P}` rep. :math:`I`.
        :rtype: (HaPPPy.MasterEquation.Simulation, HaPPPy.MasterEquation.Simulation)
        """

        ## type conversions
        # if necessary: reformat P_0 as nx1 matrix (otherwise P_0 could not be multiplicated with matricies)
        if P_0.ndim == 1:
            P_0 = np.array([P_0]).transpose()

        ## input checks
        # get dimension
        n = P_0.shape[0]
        # warn if tolerance value is unreasonable
        if check_tolerance and self.ε <= 0:
            raise RuntimeError("ε must be a positive number > 0. \nε = " + str(self.ε))
        # P_0 must be a nx1 matrix (with n matching Λ)
        if P_0.ndim != 2 or P_0.shape[1] != 1:
            raise RuntimeError("P_0 must be a "
                               + str(n)
                               + "x1 matrix (aka. a 'dotable vector')! \nP_0 = \n"
                               + str(P_0)
                              )
        # Γ_L must be a nxn matrix
        if Γ_L.ndim != 2 or Γ_L.shape[0] != P_0.shape[0] or Γ_L.shape[1] != P_0.shape[0] or not (Γ_L >= 0).all():
            raise RuntimeError("Γ_L must be a "
                               + str(n) + "x" + str(n)
                               + " matrix with coefficients >= 0."
                               + "\nΓ_L = \n" + str(Γ_L)
                              )
        # Γ_R must be a nxn matrix
        if Γ_R.ndim != 2 or Γ_R.shape[0] != P_0.shape[0] or Γ_R.shape[1] != P_0.shape[0] or not (Γ_R >= 0).all():
            raise RuntimeError("Γ_R must be a "
                               + str(n) + "x" + str(n) +
                               " matrix with coefficients >= 0."
                               + "\nΓ_R = \n" + str(Γ_R)
                              )
        # P_0 must have coefficients >= 0
        if not (P_0 >= 0).all():
            raise RuntimeError("P_0 must have coefficients >= 0. \nP_0 = \n" + str(P_0))
        # coefficients of P_0 must add up to 1
        P_sum = sum(P_0)
        if check_tolerance and (P_sum < 1 - self.ε or P_sum > 1 + self.ε):
            raise RuntimeError("Coefficients of P_0 must add up to 1 (within tolerance ε = "
                               + str(self.ε)
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
        sim_tdp_successful, sim_tdp = self.__simulateTimeDevelopmentOfPropabilities(Δt, t_max, P_0, Γ, check_tolerance)
        if sim_tdp_successful:
            sim_cur = self.__simulateCurrent(Γ_L, Γ_R, sim_tdp)

        return sim_tdp, sim_cur

    def __simulateTimeDevelopmentOfPropabilities(self,
                                                 Δt,
                                                 t_max,
                                                 P_0,
                                                 Γ,
                                                 check_tolerance=True,
                                                 verbose=False
                                                ):
        """
        Simulates the time development of propabilities P(t) for an :math:`n`-state system with transition rates Γ.

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
                print("ε = ", self.ε)

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
        valid = True
        sim = Simulation(Δt, t_max)
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
            sim.append(np.array([P_t[i][0] for i in range(P_t.shape[0])]))
            # check if sum of coefficients of P_t is still 1
            P_sum = sum(P_t)
            if check_tolerance and (P_sum < 1 - self.ε or P_sum > 1 + self.ε):
                print("Warning! Calculation aborted Coefficients of P_(t="
                      + str(t)
                      + ") must add up to 1 (within tolerance ε = "
                      + str(self.ε)
                      + "). \nP_(t="
                      + str(t)
                      + ") = \n"
                      + str(P_t)
                      + "\n Σ = "
                      + str(P_sum[0])
                     )
                valid = False
                break

        if valid:
            sim.validate()

        # Use either:
        #   return sim_successful ? sim : None
        # or
        #   return sim_successful, sim
        # to return the simulated data. The last option allowes to analyse the
        # simulation up to the point of failure but must be analysed more carefully.
        return sim_successful, sim

    def __simulateCurrent(self,
                           Γ_L,
                           Γ_R,
                           sim_time_dev_prop,
                           verbose=False
                          ):
        """
        Simulates the time current I(t) for an :math:`n`-state system with transition rates Γ.

        .. todo::

            Create extra return type for simulation.

        This is a private function and should not be called from outside the MasterEquationSolver class! (Hence it does not perform further tests than MasterEquationSolver.)

        :return: Returns :code:`sim_current`. :code:`sim_time_dev_prop` needs to be valid for :code:`sim_current` to be valid!
        :rtype: (bool, list)

        See documentation of MasterEquationSolver.doCalculation (2nd problem) for more details.

        """

        valid = sim_time_dev_prop.valid()
        sim = Simulation(sim_time_dev_prop.getΔt(), sim_time_dev_prop.getT_max())

        # calculate ΔΓ (direction sensivtive form)
        ΔΓ = Γ_L - Γ_R
        if verbose: print("ΔΓ: \n", ΔΓ)

        # simulate the current discretly
        for i in range(len(sim_time_dev_prop)):
            I_i = sum(sum(np.dot(ΔΓ, sim_time_dev_prop[i][1])))
            sim.append(I_i)
            if verbose: print("I_t[" + str(i) + "= \n", I_i)

        if valid:
            sim.validate()

        return sim

class Simulation():
    """
    This class represents a simulation crated by HaPPPy.MasterEquation.MasterEquationSolver.doCalculation().
    A function :math:`f : [0, t_{max}] \\to X, t \mapsto f(t)` is approximated by calculating its discrete values :math:`f(n \\cdot \\Delta t)` where :math:`n \in \mathbb{N}_0 \land n \cdot \Delta t \leq t_{max}`.

    If :code:`sim` is a valid HaPPPy.MasterEquation.Simulation then:

    :code:`sim.getTimeBins()` ≣ :math:`\\{n \\cdot \\Delta t | n \in \mathbb{N}_0 \land n \cdot \Delta t \leq t_{max} \\}`

    and

    :code:`sim.getValues()` ≣ :math:`\\{f(n \\cdot \\Delta t) | n \in \mathbb{N}_0 \land n \cdot \Delta t \leq t_{max} \\}`.

    Both in increasing order of :math:`n`.
    """

    def __init__(self, Δt, t_max):
        """
        :param Δt: :code:`Δt` ≣ :math:`\\Delta t > 0`
        :type: float
        :param t_max: :code:`t_max` ≣ :math:`t_{max} \\geq 0`
        :type: float
        """
        self.Δt = Δt
        self.t_max = t_max
        self.__values= []
        self.__valid = False

    def getΔt(self):
        """
        :return: Returns :code:`Δt` ≣ :math:`\\Delta t`.
        :rtype: float
        """
        return self.Δt

    def getT_max(self):
        """
        :return: Returns :code:`t_max` ≣ :math:`t_{max}`.
        :rtype: float
        """
        return self.t_max

    def getTimeBins(self):
        """
        :return: All values :math:`n \\cdot \\Delta t` where :math:`n \in \mathbb{N}_0 \land n \cdot \Delta t \leq t_{max}` with increasing :math:`n`.
        :rtype: numpy.ndarray
        """
        return np.arange(0, self.t_max + self.Δt, self.Δt)

    def getValues(self):
        """
        :return: All values :math:`f(n \\cdot \\Delta t)` where :math:`f` is the simulated function and :math:`n \in \mathbb{N}_0 \land n \cdot \Delta t \leq t_{max}` with increasing :math:`n`.
        :rtype: numpy.ndarray
        """
        return np.array(self.__values)

    def append(self, value):
        # for internal use only
        self.__values.append(value)

    def valid(self):
        """
        :return: True iff simulation was not abroted.
        :rtype: bool
        """
        return self.__valid

    def validate(self):
        # for internal use only
        self.__valid = True

    def __getitem__(self, n):
        """
        :param n: :code:`n` ≣ :math:`n`.
        :return: Returns :math:`f(n \\cdot \\Delta t)`.
        :rtype: float or numpy.ndarray
        """
        return self.__values[n]

    def __len__(self):
        """
        :return: Number of (current) calculated values in simulation.
        :rtype: integer
        """
        return len(self.__values)

    def __repr__(self):
        ts = self.getTimeBins()
        vs = self.getValues()
        return "Sim<" + str([(ts[i], vs[i]) for i in range(len(self))]) + ">"

    def __str__(self):
        return repr(self)

    def quickPlot(self, title=None, xlabel=None, ylabel=None, xunit=None, yunit=None):
        """
        Simple plotting method to quickly get an overview on the simulation.

        :param title: The title of the plot. (optional)
        :type title: string
        :param xlabel: The symbol to retresent the parameter. (optional)
        :type xlabel: string
        :param ylabel: The symbol to retresent the function values. (optional)
        :type ylabel: string
        :param xunit: The unit the paramter is meassured in. (optional)
        :type xunit: string
        :param yunit: The unit the function valuesare meassured in. (optional)
        :type yunit: string

        :example: See the example given at the documentation of HaPPPy.MasterEquation.MasterEquationSolver.
        """
        # auquire values
        ts = self.getTimeBins()
        vs = self.getValues()
        if (not self.valid())and len(vs) < len(ts):
            ts = ts[:len(vs)]

        # plot
        plt.plot(ts, vs)
        # validity
        if not self.valid():
            if title == None:
                title = " (not valid)"
            else:
                title = str(title) + " (not valid)"
        # title
        if title != None:
            plt.title(str(title))
        # labels
        if xlabel != None:
            if xunit != None:
                plt.xlabel(str(xlabel) + "/" + str(xunit))
            else:
                plt.xlabel(str(xlabel))
        if ylabel != None:
            if yunit != None:
                plt.ylabel(str(ylabel) + "/" + str(yunit))
            else:
                plt.ylabel(str(ylabel))
        # legend
        if ylabel != None and (type(vs[0]) == list or type(vs[0]) == np.ndarray) and len(vs[0]) > 1:
            n = len(vs[0]) # dimension of v
            legend = ["$" + str(ylabel) + "_" + str(i) + "$" for i in range(n)]
            plt.legend(legend)
        # grid
        plt.grid()
        # show
        plt.show()
