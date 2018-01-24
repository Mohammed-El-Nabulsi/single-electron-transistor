
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
import math
import copy


__docformat__ = 'reStructuredText'

class MasterEquationSolver:
    """

    Simulates the time development of propabilities :math:`\\vec{P}(t)` or finds
    stationary solutions :math:`\\vec{P}_{stat}` for an :math:`n`-state quantum
    dot transistor with two reservoirs with transition rates :math:`\\Gamma_{L}`
    and :math:`\\Gamma_{R}` including the netto current :math:`I(t)` or the
    sationary netto current :math:`I_{stat}`.


    calculateStationarySloutions() returns all possible **sationary** solutions
    and simulateDynamicSloution() returns a numerical approximated and discrete
    **dynamic** simulation of both quantities.

    The problem of calculating the sationary solutions or a dynamic simulation
    is divided into two steps.

    **1st Problem: Propabilities**

    Let :math:`n_i \in \mathbb{N}` be the number of states of :math:`i`
    particles in the quantum dot where :math:`i` runns from :math:`0` to
    :math:`m` and :math:`n = \\sum_{i=0}^m n_i`.
    The propability of the system to be in state :math:`\\alpha \\in
    \\{1, \\dots n\\}` is called :math:`P_{\\alpha}(t)` which can variy with the
    time :math:`t`. They are summarized as a vector :math:`\\vec{P}(t)`.
    The rates for any such state to transit to any other state is represented as
    :math:`\Gamma = \Gamma^L + \Gamma^R + \Gamma^0 \in
    Mat(n, n, \mathbb{R}^+_0)` where :math:`\Gamma_{\\alpha \\beta}
    \equiv \Gamma_{\\alpha \\rightarrow \\beta}` denotes the rate from state
    :math:`\\alpha` to :math:`\\beta`. There are two different matricies
    :math:`\\Gamma^L` and :math:`\\Gamma^R` to describe the transitions whitch
    include exchanges of particles with the two reservoirs
    (called *left* and *right*). The third matrix :math:`\Gamma^0` is
    *optional* and describes the rates for relaxation and stimualtion which do
    not require any interaction with the reservoirs. Since for this problem only
    the total rates matter they are added up. (The idividual rates are only
    relevant to the second problem.)

    The time development of propabilities are derived from the master equation
    (a differential equation of first order with constant coefficients).

    .. math::

        \\frac{d P_{\\alpha}}{d t}
        = \sum_{\\beta} \Gamma_{\\beta \\rightarrow \\alpha} P_{\\beta}
        - \sum_{\\beta} \Gamma_{\\alpha \\rightarrow \\beta} P_{\\alpha}

    It states that the propability of the quantum dot to be in state
    :math:`\\alpha` changes within an infinitesimal time interval depending on
    the current state and the trasition rates to another state. All influences
    can be sorted into two categories. The left sum describes the total *gain*
    in propability whereas the right sum describes the total *loss*. If a state
    :math:`\\beta` has some propability and a transition rate
    :math:`\Gamma_{\\beta \\rightarrow \\alpha}` other than :math:`0`, state
    :math:`\\alpha` gains some propability proportional to this rate. All
    these gains add up. Likewise sate :math:`\\alpha` loses propability to
    any other state :math:`\\beta` proporitional to
    :math:`\Gamma_{\\alpha \\rightarrow \\beta}` (hence the minus sign).

    This problem can be restated as a simple linear equation:

    .. math::
        \\begin{pmatrix}
            \\dot{P}_1 \\\\ \\dot{P}_2 \\\\ \\vdots \\\\ \\dot{P}_n
        \\end{pmatrix}
        =
        \\underbrace{ \\left(
            \\begin{pmatrix}
                \\Gamma_{11} & \\Gamma_{21} & \\dots  & \\Gamma_{n1} \\\\
                \\Gamma_{12} & \\Gamma_{22} & \\dots  & \\Gamma_{n2} \\\\
                \\vdots      & \\vdots      & \\ddots & \\vdots      \\\\
                \\Gamma_{1n} & \\Gamma_{2n} & \\dots  & \\Gamma_{nn}
            \\end{pmatrix}
            -
            \\begin{pmatrix}
                \\sum^{n}_{\\beta = 1}{\\Gamma_{1 \\beta}} & 0 & \\dots & 0 \\\\
                0 & \\sum^{n}_{\\beta = 1}{\\Gamma_{2 \\beta}} & \\dots & 0 \\\\
                \\vdots & \\vdots & \\ddots & \\vdots \\\\
                0 & 0 & \\dots & \\sum^{n}_{\\beta = 1}{\\Gamma_{n \\beta}}
            \\end{pmatrix}
        \\right) }_{:= \\Lambda}
        \\begin{pmatrix}
            P_1 \\\\ P_2 \\\\ \\vdots \\\\ P_n
        \\end{pmatrix}

    When time development is considered the solution to this differntial
    equation is:

    .. math::
        \\vec{P}(t) = \\exp(\\Lambda t) \\cdot \\vec{P}(t=0)

    This equation becomes almost trivial when :math:`\\Lambda` is diagonalized:
    Let :math:`\\vec{v}_i(t=0)` be an eigenvector of :math:`\\Lambda` with
    eigenvalue :math:`\\lambda_i`. It's solution is:

    .. math::
        \\vec{v}_i(t) = \\exp(\\lambda_i t) \\cdot \\vec{v}_i(t=0)

    Note that sationary solutions :math:`\\vec{P}_{stat,i}` have eigenvalues
    :math:`\\lambda_i = 0`. All other eigenvectors have eigenvalues
    :math:`\\lambda_i < 0` and therefore decrese exponentially over time.

    **Stationary** solutions require :math:`\\vec{P}(t)` to not vary in time
    which is expressed as :math:`\\dot{\\vec{P}} = 0`. Therfore all linear
    superpositions of :math:`\\vec{P}_{stat,i}` with
    :math:`\\dot{\\vec{P}} = \\Lambda \\vec{P}_{stat,i} = 0` are stationary.
    (Notice that all linear superpositions need to fulfil
    :math:`0 \\leq P_i \\leq 1` and the normalization condition
    :math:`\\sum^n_{i=1}{P_i} = 1` to be valid.)

    A basis of normalized sationary solutions is returned as a matrix.

    **Dynamic** simulations are calculated from a (valid) distribution of
    propability :math:`\\vec{P}(t=0)`. It approximates :math:`\\vec{P}(t)`
    in the time interval :math:`0 \\leq t \\leq t_{max} \in \mathbb{R}^+_0` in
    discrete steps of time :math:`t = i \\cdot \\Delta t` where
    :math:`\\Delta t \\in \mathbb{R}^+` and :math:`i \\in \\mathbb{N}`.

    A dynamic simulation is returned as a special Simulation type.
    (See HaPPPy.MasterEquation.Simulation for more details.)

    **2nd Problem: Netto Current**

    For each point in time the netto current :math:`I^k(t)` flowing trough each
    barrier of the quantum dot can be calculated. This does not depend on the
    type of solution (stationary or dynamic).

    For each barrier between the quantum dot and a reservoir a current
    :math:`I^k` results from the tunneling processes of charged perticles
    through them. Each particle carries as charge :math:`q`. The flow from the
    left to the right is considered to be positive (if :math:`q` > 0).

    The netto current for both barriers are:

    .. math::
        I^L = q \\sum_{\\alpha, \\beta}
                   \\sigma_{\\alpha \\rightarrow \\beta}
                   \\Gamma^L_{\\alpha \\rightarrow \\beta}
                   P_{\\alpha}

    and

    .. math::
        I^R = -q \\sum_{\\alpha, \\beta}
                    \\sigma_{\\alpha \\rightarrow \\beta}
                    \\Gamma^R_{\\alpha \\rightarrow \\beta}
                    P_{\\alpha}


    where

    .. math::
        \\sigma_{\\alpha \\rightarrow \\beta}
        =   \\begin{cases}
                +1 & | \\text{ if the transition }
                       \\alpha \\rightarrow \\beta
                       \\text{ adds a particle (*)}\\\\
                 0 & | \\text{ if the transition }
                       \\alpha \\rightarrow \\beta
                       \\text{ adds or removes no particle (*)}\\\\
                -1 & | \\text{ if the transition }
                       \\alpha \\rightarrow \\beta
                       \\text{ removes a particle (*)}
            \\end{cases}

    (*) with respect to the left barrier.

    During the calculation it is assumed that the charge per particle is
    :math:`q = 1` which is equal to divide the given formulae by :math:`q`.

    **Stationary** solutions contains both total currents for each basis vector
    of the stationary solution from the first problem (in the same order).

    These solutions are returned as a list.

    **Dynamic** simulations are calculated in the same manner as in the
    simulation of problem one. Both currents are returned.

    The solution of the simulation is returned as a Simulation as well.
    (See HaPPPy.MasterEquation.Simulation for more details.)

    :example: .. code-block:: Python
            :emphasize-lines: 1-2, 6-8, 10, 12-16

            import HaPPPy
            import numpy as np

            ## simple simulation with the MasterEquationSolver class
            # set-up a reasonable Γ-matrix
            Γ_L = np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
            Γ_R = np.array([[0, 0.0, 0.5], [0.0, 0, 0.5], [0.0, 0.5, 0]])
            ns = [1, 2]
            # choose a legitimate start value for P_0 (P_0 = P(t=0))
            P_0 = np.array([0.9, 0.0, 0.1])
            # simulate
            Δt = 1
            t_max = 100
            mes = HaPPPy.MasterEquation.MasterEquationSolver()
            sim_tdp, sim_cur = mes.simulateDynamicSloution(Δt, t_max, P_0, Γ_L, Γ_R, ns)
            stat_ps, stat_curs = mes.calculateStationarySloutions(Γ_L, Γ_R, ns)
            ## plot/print results
            # 1st plot: time development of propabilities
            sim_tdp.quickPlot(xlabel="t", ylabel="P")
            print("static solutions (P_ij) =\\n", stat_ps)
            # 2nd plot: time development of netto current
            sim_cur.quickPlot(xlabel="t", ylabel="I",legend=["$I^L$","$I^R$"])
            print("static solutions (I_i) = \\n", stat_curs)


        The relevant lines of code for the simulation to work are highlighted.
        To give a real live example and to demonstarte the usage of the result
        some code to plot the result is added.

    """
    # Don't remove the last dot/paragrph since it is a workaround to shpinx's
    # code-block interpreter.

    def __init__(self, ε=None):
        """
        :param ε: A tolerance value. Sum of all probabilities must be equals to
                  1 within this tolerance. If :code:`ε == None` the default
                  value is used (see getDefaultε()).
        :type ε: float
        """
        if ε == None:
            self.setε(self.getDefaultε())
        else:
            self.ε = ε

    def getε(self):
        """
        :return: Returns ε - a tolerance value: The sum of all probabilities
                 must be equals to 1 within this tolerance.
        :rtype: float
        """
        return self.ε

    def setε(self, ε):
        """
        :param ε: A tolerance value: The sum of all probabilities must be equals
                  to 1 within this tolerance.
        :type ε: float
        """
        self.ε = ε

    @staticmethod
    def getDefaultε():
        """
        :return: Returns the default value of ε - a tolerance value: The sum of
                 all probabilities must be equals to 1 within this tolerance.
        :rtype: float
        """
        return 1E-10

    def calculateStationarySloutions(self,
                                     Γ_L, Γ_R,
                                     ns,
                                     Γ_0 = None,
                                     check_tolerance=True,
                                     verbose=False,
                                    ):
        """
        Calculates all possible stationary solutions :math:`\\vec{P}_{stat,i}`
        (with :math:`\\Lambda \\cdot \\vec{P}_{stat,i} = 0`). The result is
        returned as a :math:`n \\times k`-matrix :math:`\\Pi` formated
        like:

        .. math::
            \\Pi = \\begin{pmatrix}
            \\vec{P}_{stat,1} & \\vec{P}_{stat,2} & \\dots & \\vec{P}_{stat,k}
            \\end{pmatrix}

        Note: The total amount of stationary solutions
        :math:`k \\in \\mathbb{N}` depends on the input.

        For more details see class description of Happy.MasterEquationSolver.

        :param Γ_L: Matrix containing the transition rates regarding electrons
                    tunneling trough the *left* barrier where
                    :code:`Γ_L[i,j]` ≣ :math:`\Gamma^L_{i \\rightarrow j}`.
                    Must be a ``nxn`` matrix.
        :type Γ_L: numpy.ndarray
        :param Γ_R: Matrix containing the transition rates regarding electrons
                    tunneling trough the *right* barrier where
                    :code:`Γ_R[i,j]` ≣ :math:`\Gamma^R_{i \\rightarrow j}`.
                    Must be a ``nxn`` matrix.
        :type Γ_R: numpy.ndarray
        :param ns: List of numbers of states :math:`n_i` with the same number
                   :math:`i` of electrons in inreasing order of :math:`i`:

                   :code:`ns = [n_0, n_1, …, n_k]`
                   ≣ :math:`(n_0, n_1, \\dots, n_k)`.

                   E.g.
                   :code:`ns = [<numer of empty states>, <number of one body states>, ...]`
        :type ns: list(int)
        :param Γ_0: Matrix containing the transition rates for relaxation or
                    stimulation where
                    :code:`Γ_0[i,j]` ≣ :math:`\Gamma^0_{i \\rightarrow j}`.
                    Must be a ``nxn`` matrix. (optinal)
        :type Γ_0: numpy.ndarray

        :return: Returns :code:`stat_ps, stat_curs` where
                 :code:`stat_ps[i,j]` ≣ :math:`P_{stat,ij}` and
                 :code:`stat_curs[i,k]` ≣ :math:`I^k_{stat,i}`.
        :rtype: numpy.ndarray(float), numpy.ndarray(float)
        """

        ## STEP 1: PARAMETER CHECKS
        # In this step some general checks are applied to the input parameters
        # to warn the user if the passed parameters are unreasonable and raise
        # an RuntimeError exception if so.

        # The dimension n is extracted from an arbitrary matrix Γ_i.
        n = Γ_L.shape[0]
        # Check and warn if Γ_L is not a nxn matrix with positive coeffiecents.
        if (Γ_L.ndim != 2
            or Γ_L.shape[0] != n
            or Γ_L.shape[1] != n
            or not (Γ_L >= 0).all()
           ):
            raise RuntimeError("Γ_L must be a "
                               + str(n) + "x" + str(n)
                               + " matrix with coefficients >= 0."
                               + "\nΓ_L = \n" + str(Γ_L)
                              )
        # Check and warn if Γ_R is not a nxn matrix with positive coeffiecents.
        if (Γ_R.ndim != 2
            or Γ_R.shape[0] != n
            or Γ_R.shape[1] != n
            or not (Γ_R >= 0).all()
           ):
            raise RuntimeError("Γ_R must be a "
                               + str(n) + "x" + str(n) +
                               " matrix with coefficients >= 0."
                               + "\nΓ_R = \n" + str(Γ_R)
                              )

        ## STEP 2: CALCULATION
        # In this step the actual calculation happens.

        # Calculate the eigenvectors and eigenvalues of Λ.
        Λ_evals, Λ_evecs = self.__calculateLambda(Γ_L, Γ_R, Γ_0,
                                                  check_tolerance,
                                                  verbose,
                                                 )

        # Calls the simulation method of the propabilities.
        stat_ps = self.__calculateStationaryPossibilitySolutions(Λ_evals,
                                                                Λ_evecs,
                                                                check_tolerance,
                                                                verbose,
                                                                )
        if verbose: print("stat_ps =\n" + str(stat_ps))

        # Calls the simulation method of the netto current.
        # Γ_L and Γ_R are passed seperatly to preserve the information of
        # direction.
        stat_curs = self.__calculateStationaryCurrents(Γ_L, Γ_R,
                                                       ns,
                                                       stat_ps,
                                                       verbose,
                                                      )
        if verbose: print("stat_curs =\n" + str(stat_curs))

        # Return stationary solutions.
        return stat_ps, stat_curs

    def simulateDynamicSloution(self,
                                Δt, t_max,
                                P_0,
                                Γ_L, Γ_R,
                                ns,
                                Γ_0 = None,
                                check_tolerance=True,
                                verbose=False,
                               ):
        """
        Calculates two simulations: the propability :math:`\\vec{P}(t)` and
        the current :math:`I(t)`.

        Raises an exception of the type RuntimeError when the calculation must
        be abored due to bad or unreasonable parameters. If the simulation is
        aborded during the calculation process later on an incomplete answer is
        still returned and this can be check via the method :code:`valid()` (see
        HaPPPy.MasterEquation.Simulation.valid()). The invalid simulation
        contains all calculated datapoints up to the point in time where the
        calculation failed. Leaving the tolverance of :math:`\\epsilon`
        is considered as a failure. E.g. it is assured that
        :math:`\\left| \\sum^n_{i=1} P_i - 1 \\right| < \\epsilon`.

        For more details see class description of Happy.MasterEquationSolver.

        :param P_0: Start value of propabilities where
                    :code:`P_0` ≣ :math:`\\vec{P_0}`.
                    Must be either a list or a ``nx1`` matrix.
        :type P_0: numpy.array or numpy.ndarray
        :param Γ_L: Matrix containing the transition rates regarding electrons
                    tunneling trough the *left* barrier where
                    :code:`Γ_L[i,j]` ≣ :math:`\Gamma^L_{i \\rightarrow j}`.
                    Must be a ``nxn`` matrix.
        :type Γ_L: numpy.ndarray
        :param Γ_R: Matrix containing the transition rates regarding electrons
                    tunneling trough the *right* barrier where
                    :code:`Γ_R[i,j]` ≣ :math:`\Gamma^R_{i \\rightarrow j}`.
                    Must be a ``nxn`` matrix.
        :type Γ_R: numpy.ndarray
        :param ns: List of numbers of states :math:`n_i` with the same number
                   :math:`i` of electrons in inreasing order of :math:`i`:

                   :code:`ns = [n_0, n_1, …, n_k]`
                   ≣ :math:`(n_0, n_1, \\dots, n_k)`.

                   E.g.
                   :code:`ns = [<numer of empty states>, <number of one body states>, ...]`
        :type ns: list(int)
        :param Γ_0: Matrix containing the transition rates for relaxation or
                    stimulation where
                    :code:`Γ_0[i,j]` ≣ :math:`\Gamma^0_{i \\rightarrow j}`.
                    Must be a ``nxn`` matrix. (optinal)
        :type Γ_0: numpy.ndarray
        :param t_max: Last point in time to be simulated. Must fulfil
                      :code:`t_max >= 0`.
        :type t_max: float
        :param Δt: Length of the time tintervall between two simulated events.
                   Must fulfil :code:`Δt > 0`.
        :type Δt: float

        :return: Returns :code:`sim_time_dev_prop, sim_current` where both
                 values are a HaPPPy.MasterEquation.Simulation of
                 :math:`\\vec{P}(t)` rep. :math:`I(t)`.
        :rtype: (HaPPPy.MasterEquation.Simulation,
                 HaPPPy.MasterEquation.Simulation)
        """

        ## STEP 1: TYPE CONVERSION OF PARAMETERS
        # In this step parameters are manipulated to match the desired type.

        # P_0 can be passed as either a vector with n components or as a nx1
        # matrix - both are a numpy.ndarray but of different shapes. Since any
        # vector - one-dimensional numpy.ndarray - can not be transposed, the
        # matrix representation is used from now on. Hence a vector is
        # transformed to a nx1 matrix.
        if P_0.ndim == 1:
            P_0 = np.array([P_0]).transpose()

        ## STEP 2: PARAMETER CHECKS
        # In this step some general checks are applied to the input parameters
        # to warn the user if the passed parameters are unreasonable and raise
        # an RuntimeError exception if so.

        # The dimension n is extracted from the P_0 vector.
        n = P_0.shape[0]
        # Check and warn if tolerance ε value is unreasonable.
        if check_tolerance and self.ε <= 0:
            raise RuntimeError("ε must be a positive number > 0."
                               + "\nε = " + str(self.ε)
                              )
        # Check and warn if ns does not add up to n.
        if sum(ns) != n:
            raise RuntimeError("The total number of states must be "
                               + str(n)
                               + "."
                               + "\n ns = "
                               + str(ns)
                               + "\n"
                               + " Σ = "
                               + str(sum(ns))
                               + " != "
                               + str(n)
                              )
        # Check and warn if P_0 is not a nx1 matrix.
        if P_0.ndim != 2 or P_0.shape[1] != 1:
            raise RuntimeError("P_0 must be a "
                               + str(n)
                               + "x1 matrix (aka. a 'dotable vector')!\nP_0 =\n"
                               + str(P_0)
                              )
        # Check and warn if Γ_L is not a nxn matrix with positive coeffiecents.
        if (Γ_L.ndim != 2
            or Γ_L.shape[0] != n
            or Γ_L.shape[1] != n
            or not (Γ_L >= 0).all()
           ):
            raise RuntimeError("Γ_L must be a "
                               + str(n) + "x" + str(n)
                               + " matrix with coefficients >= 0."
                               + "\nΓ_L = \n" + str(Γ_L)
                              )
        # Check and warn if Γ_R is not a nxn matrix with positive coeffiecents.
        if (Γ_R.ndim != 2
            or Γ_R.shape[0] != n
            or Γ_R.shape[1] != n
            or not (Γ_R >= 0).all()
           ):
            raise RuntimeError("Γ_R must be a "
                               + str(n) + "x" + str(n) +
                               " matrix with coefficients >= 0."
                               + "\nΓ_R = \n" + str(Γ_R)
                              )
        # Check and warn if Γ_0 is not a nxn matrix with positive coeffiecents.
        # Γ_0 is optional!
        if (Γ_0 != None and (Γ_0.ndim != 2
                             or Γ_0.shape[0] != n
                             or Γ_0.shape[1] != n
                             or not (Γ_0 >= 0).all()
                            )
           ):
            raise RuntimeError("Γ_0 must be a "
                               + str(n) + "x" + str(n) +
                               " matrix with coefficients >= 0."
                               + "\nΓ_0 = \n" + str(Γ_0)
                              )
        # Check and warn if coefficients of P_0 are not positive.
        if not (P_0 >= 0).all():
            raise RuntimeError("P_0 must have coefficients >= 0.\nP_0 ="
                               + "\n" + str(P_0)
                              )
        # Check and warn if coefficients of P_0 do not add up to 1.
        P_sum = sum(P_0)
        if check_tolerance and (P_sum < 1 - self.ε or P_sum > 1 + self.ε):
            raise RuntimeError("Coefficients of P_0 must add up to 1 "
                               + "(within tolerance ε = " + str(self.ε) + ")."
                               + "\nP_0 =\n" + str(P_0)
                               + "\n Σ = " + str(P_sum)
                              )
        # Check and warn if t_max is negative.
        if t_max < 0:
            raise RuntimeError("Simulated time must be positive or zero "
                               + "(t_max >= 0) but passed value is "
                               + "t_max = " + str(t_max)
                              )
        # Check and warn if Δt is negative or zero.
        if Δt <= 0:
          raise RuntimeError("Simulated time intervals must be positive "
                             + "(Δt > 0) but passed value is "
                             + "Δt = " + str(Δt)
                            )

        ## STEP 3: SIMULATION
        # In this step the actual simulation happens.

        # Calculate the eigenvectors and eigenvalues of Λ.
        Λ_evals, Λ_evecs = self.__calculateLambda(Γ_L, Γ_R, Γ_0,
                                                  check_tolerance,
                                                  verbose,
                                                 )
        if verbose: print("Λ_evals =\n" + str(Λ_evals))
        if verbose: print("Λ_evecs =\n" + str(Λ_evecs))

        # Calls the simulation method of the propabilities.
        sim_tdp = self.__simulateDynamicPossiblitySolution(Δt,
                                                 t_max,
                                                 P_0,
                                                 Λ_evals, Λ_evecs,
                                                 check_tolerance,
                                                 verbose,
                                                )
        if verbose: print("sim_tdp =\n" + str(sim_tdp))

        # Calls the simulation method of the netto current.
        # Γ_L and Γ_R are passed seperatly to preserve the information of
        # direction.
        sim_cur = self.__simulateDynamicCurrent(Γ_L, Γ_R, ns, sim_tdp, verbose)
        if verbose : print("sim_cur =\n" + str(sim_cur))

        # Both simulations are returnd - mind the order!
        return sim_tdp, sim_cur

    def __calculateLambda(self, Γ_L, Γ_R, Γ_0, check_tolerance, verbose):
        """
        Takes both matricies :math:`\\Gamma^k`and calculates :math:`\\Lambda` in
        diagonalized form: a list of eigenvalues :math:`\\lambda_i` and a
        :math:`n \\times n`-matrix containing all eigenvectors
        :math:`\\vec{P}_i` formated like:

        .. math::
            \\Eta =
            \\begin{pmatrix}
                \\vec{P}_1 & \\vec{P}_2 & \\dots & \\vec{P}_n
            \\end{pmatrix}

        where

        .. math::
            \\Lambda \\vec{P}_i = \\lambda_i \\vec{P}_i

        This is a private function and should not be called from outside the
        MasterEquationSolver class! (Hence it does not perform further tests.)

        :return: Returns :code:`Λ_evals, Λ_evecs` ≣
                         :math:`(\\lambda_1, \\lambda_2, \\dots, \\lambda_n),
                         \\Eta`.
        :rtype: numpy.ndarray(float)

        See documentation of MasterEquationSolver.doCalculation (1st problem)
        for more details.

        """

        ## Disclaimer: The usage of `verbose` is for debugging only. If set to
        #              `True` most parts of the calculation are printed to the
        #              console. All expressions of code realting to this process
        #              begin  with `if verbose: ...`. Removing those expressions
        #              does not alter the behaviour of this method if verbose is
        #              set to `False`.

        # Γ stores the total rates - in a non-directional manner.
        Γ = Γ_L + Γ_R
        if Γ_0 != None:
            Γ += Γ_0

        # verbose only: Print Γ (and ε if tolerence checks are enabled).
        if verbose:
            print("Γ = \n", Γ)
            if check_tolerance:
                print("ε = ", self.ε)

        ## STEP 1: CALCULATION OF THE COEFFICIENTS OF THE DIFFERENTIAL EQUAITON
        #          (MASTER EQUATION) IN EIGENVEVTORBASIS

        ## STEP 1.a: CALCULATION OF THE COEFFICIENTS OF THE MASTER EQUAITON
        # Clalculate the coefficients of the master equation in the given basis.
        # Λ_in holds all positive coefficients in master equation - the inflow.
        # (Λ_in is effectivly the left sum in master equation in the
        # documentation)
        Λ_in = Γ
        # And Λ_in holds all negative coefficients in master equation - the
        # outflow. (Λ_out is effectivly the right sum in master equation in the
        # documentation)
        Λ_out = np.diag([sum([Γ[j][i]
                        for j in range(Γ.shape[0])])
                        for i in range(Γ.shape[0])]
                       )
        # Combine both in to one matrix - the netto flow.
        Λ = Λ_in - Λ_out
        # verbose only: Prints all newly calculated quantities.
        if verbose: print("Λ_in = \n", Λ_in)
        if verbose: print("Λ_out = \n", Λ_out)
        if verbose: print("Λ = \n", Λ)

        ## STEP 1.b: RESTATE THE MASTEREQUATION IN EIGENVECTORBASIS
        # Find an eigenvector basis of Λ and their eigenvalues. (Λ_evecs_T
        # transforms from the given basis to an eigenvevtor basis of Λ.)
        # Convention: Let v be an eigenvector of Λ with eigen value λ: Λv = λv
        #             Then: Λ_evals = (λ_1, λ_2, ..., λ_n)
        #             and   Λ_evecs = (v_1, v_2, ..., v_n)
        # Convention: The calculated eigenvetorbasis of Λ by this step is called
        #             `the` eigenvector basis from now on.
        (Λ_evals, Λ_evecs) = lin.eig(Λ)
        # verbose only: Print the eigenvector basis and the eigenvalues.
        if verbose: print("Λ.eigenvalues = \n", Λ_evals)
        if verbose: print("Λ.eigenvectors = \n", Λ_evecs)

        ## STEP 1.c: PLAUSIBILITY CHECK (Λ SHOULD BE INERVIBLE)
        # Claulate the inverse of Λ in the eigenvector base.
        Λ_evecs_inv = np.linalg.inv(Λ_evecs)
        # verbose only: Print the inverse of Λ in the eigenvector base.
        if verbose: print("Λ.eingenvectors^-1 = \n", Λ_evecs_inv)
        # Check an warn if this matric does not truly invert Λ.
        if (np.dot(Λ_evecs_inv, Λ_evecs) != np.identity(Λ.shape[0])).all():
            raise RuntimeError("Could not invert Λ = \n" + str(Λ))

        return Λ_evals, Λ_evecs

    def __calculateSigma(self, ns, verbose, fast=False):
        """
        Calculates the σ_αβ as described in the master equation (see
        documentation).

        :param ns: List of numbers of states :math:`n_i` with the same number
                   :math:`i` of electrons in inreasing order of :math:`i`:

                   :code:`ns = [n_0, n_1, …, n_k]`
                   ≣ :math:`(n_0, n_1, \\dots, n_k)`.

                   E.g.
                   :code:`ns = [<numer of empty states>, <number of one body states>, ...]`
        :type ns: list(int)
        :param fast: Set this to :code:`True` if it is not necessary to set the
                     diagonal blocks in :math:`σ`
                     (where :math:`\\sigma_{i \\rightarrow j} = 0`) to zero.

        :return: Returns :code:`σ` where :code:`σ[i,j]`
                 ≣ :math:`\\sigma_{i \\rightarrow j}`.
        """
        # Construct a Matrix σ containing the signs representing the direction
        # of current for each rate. If the total numer of particles in the
        # channel increases - e.g. with rate Γ_αβ - its contribution has an
        # oppsite sign that the reverse process - e.g. with rate Γ_βα. For
        # relaxation of stimulation precesses - whrere the total number of
        # perticles does not change - no current flows through a barrier. All
        # of these 'signs' are zero.
        #
        # σ =
        #      n_0 |     n_1   |  ⋯ |  n_k
        # n_0 ⎧ 0  |  1  1  1  | ⋯  |  1  1⎫
        #     ------------------------------
        #     ⎪-1  |  0  0  0  | ⋯  |  1  1⎪
        # n_1 ⎪-1  |  0  0  0  | ⋯  |  1  1⎪
        #     ⎪-1  |  0  0  0  | ⋯  |  1  1⎪
        #     ------------------------------
        #     ⎪ ⋮  |  ⋮  ⋮  ⋮  | ⋱  |  ⋮  ⋮⎪
        #     ------------------------------
        # n_k ⎪-1  | -1 -1     | ⋯  |  0  0⎪
        #     ⎩-1  | -1 -1     | ⋯  |  0  0⎭
        # Set-up the ones.
        n_tot = np.sum(ns)
        σ = np.tri(n_tot).transpose() - np.tri(n_tot)
        # Make blocks on diagonal zero.
        # This step can be omitted if all diagonal blocks in Γ_L and Γ_R
        # are zero dy default. (It is a time-consuming step!)
        if not fast:
            i = 0
            for n_i in ns:
                σ[i:i+n_i, i:i+n_i] = np.zeros((n_i, n_i))
                i += n_i
        # verbose only: Pront σ.
        if verbose:
            print("σ<fast=" + str(fast) + "> =\n" + str(σ))

        return σ

    def __calculateStationaryPossibilitySolutions(self,
                                                  Λ_evals,
                                                  Λ_evecs,
                                                  check_tolerance,
                                                  verbose,
                                                 ):
        """
        Calculates all possible stationary solutions :math:`\\vec{P}_{stat,i}`
        (with :math:`\\Lambda \\cdot \\vec{P}_{stat,i} = 0`). The result is
        returned as a :math:`n \\times k`-matrix :math:`\\vec{\\Pi}` formated
        like:

        .. math::
            \\vec{\\Pi} = \\begin{pmatrix}
            \\vec{P}_{stat,1} & \\vec{P}_{stat,2} & \\dots & \\vec{P}_{stat,k}
            \\end{pmatrix}

        Note: The total amount of stationary solutions
        :math:`k \\in \\mathbb{N}` depends on the input.

        This is a private function and should not be called from outside the
        MasterEquationSolver class! (Hence it does not perform further tests.)

        For more details see class description of Happy.MasterEquationSolver and
        Happy.MasterEquationSolver.__calculateLambda().

        :param Λ_evals: A list of all eigenvalues :math:`\\lambda_i` of
                        :math`\\Lambda`.
        :type Λ_evals: numpy.ndarray
        :param Λ_evecs: The matrix :math:`\\Eta` containing all eigenvectors.
        :type Λ_evecs: numpy.ndarray

        :return: Returns :code:`Π` ≣ :math:`\\vec{\\Pi}`.
        :rtype: numpy.ndarray(float)
        """

        ## FIND STATIONARY SOLUTIONS
        # Find the stationary solution, by extracting the eigenvector
        # corresponding to the eigenvalue 0.
        stat_ps = np.compress((abs(Λ_evals) < self.ε), Λ_evecs, axis=1)
        # Renormalize all sationary Ps such that their components add up to 1.
        N = np.diag(np.power(np.sum(stat_ps, axis=0), -1))
        stat_ps = np.dot(stat_ps, N)

        return stat_ps

    def __simulateDynamicPossiblitySolution(self,
                                            Δt,
                                            t_max,
                                            P_0,
                                            Λ_evals, Λ_evecs,
                                            check_tolerance,
                                            verbose
                                           ):
        """
        Simulates the time development of the propabilities :math:`P(t)`.

        If the simulation is aborded during the calculation process later on an
        incomplete answer is still returned and can be check via the method
        :code:`valid() (see HaPPPy.MasterEquation.Simulation.valid()).
        The invalid simulation contains all calculated datapoints up to the
        point in time where the calculation failed. Leaving the tolverance value
        :math:`\\epsilon` is considered as a failure. E.g. it is assured that
        :math:`\\left| \\sum^n_{i=1} P_i - 1 \\right| < \\epsilon`.

        This is a private function and should not be called from outside the
        MasterEquationSolver class! (Hence it does not perform further tests.)

        For more details see class description of Happy.MasterEquationSolver.

        :param P_0: Start value of propabilities where
                    :code:`P_0` ≣ :math:`\\vec{P_0}`.
                    Must be either a list or a ``nx1`` matrix.
        :type P_0: numpy.array or numpy.ndarray
        :param Λ_evals: A list of all eigenvalues :math:`\\lambda_i` of
                        :math`\\Lambda`.
        :type Λ_evals: numpy.ndarray
        :param Λ_evecs: The matrix :math:`\\Eta` containing all eigenvectors.
        :type Λ_evecs: numpy.ndarray
        :param t_max: Last point in time to be simulated. Must fulfil
                      :code:`t_max >= 0`.
        :type t_max: float
        :param Δt: Length of the time intervall between two simulated events.
                   Must fulfil :code:`Δt > 0`.
        :type Δt: float

        :return: Returns a Simulation of :math:`\\vec{P}(t).`
        :rtype: HaPPPy.MasterEquation.Simulation
        """

        # verbose only: Print P_0.
        if verbose:
            print("P(t=0) =\n" + str(P_0))

        # Stores whether the simulation was completed witout anny issue.
        valid = True
        # Stores the resulting simulation.
        sim = Simulation(Δt, t_max)
        # Calculate P(t=0) in the eigenvector basis.
        P_evb_0 = np.dot(np.linalg.inv(Λ_evecs), P_0)
        # For each desired point in time the solution is calculated discretly.
        for t in sim.getBins():
            # verbose only: Print the current point in time.
            if verbose: print("\nt = ", t)

            # In this eigenvector basis Λ is diagonal and the master equation
            # becomes: d/dt v_i = Λ_ii v_i = λ_i v_i  (where i = 1,2,...,n)
            # These differential equations are decoulped and therefore can be
            # solved independently.
            # The solutions are: v_i(t) = v_i(t=0) exp(λ_i t)
            # Which can be restated as: v(t) = exp(Λt) v(t=0)
            # where v is P in the eigenvectorbasis
            # and exp(Λt) = diag(exp(λ_1 t), exp(λ_2 t), ..., exp(λ_n t))
            # Calcualte exp(Λt) in the eigenvectorbasis.
            exp_tΛ_evb = np.diag(np.exp(t * Λ_evals))
            # verbose only: Print P(t) in the eigenvectorbasis.
            if verbose:
                print("exp(" + str(t) + " * Λ) "
                      + "in Λ.eigenvectorbase = \n", exp_tΛ_evb
                     )
            # Calculate P(t) in the eigenvectorbasis - which is v(t).
            P_evb_t = np.dot(exp_tΛ_evb, P_evb_0)
            # verbose only: Print P(t) in the eigenvectorbasis.
            if verbose:
                print("P(t=" + str(t) + ") in "
                      + "Λ.eigenvectorbase = \n", P_evb_t
                     )

            # Calculate P(t) in the origianl basis.
            P_t = np.dot(Λ_evecs, P_evb_t)
            # verbose only: Print P(t) in the origianl basis.
            if verbose: print("P(t=" + str(t) + ") = \n", P_t)

            # Extend the simualtion by appending the new P(t) as a vector.
            sim.append(np.array([P_t[i][0] for i in range(P_t.shape[0])]))

            # Ceck and warn if the sum of all P_i do not add up to 1 within a
            # tolerance of ε (if requested).
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
                # If the sum of all P_i is outside the tolerance range abort
                # the simulation. (This also prevents costly and unrasonable
                # blow-ups of any P_i.) The simulation is from now on considered
                # as inavlid.
                valid = False
                break

        # Simlations are by default invalid and hece the this simulation is
        # validated iff the simulation had no issues.
        if valid:
            sim.validate()

        # Return the (valid or invalid) simulation.
        return sim

    def __calculateStationaryCurrents(self,
                                      Γ_L, Γ_R,
                                      ns,
                                      stat_ps,
                                      verbose,
                                     ):
        """
        Calculates the current :math:`I_{stat,i}` to all stationary solutions
        :math:`\\vec{P}{stat,i}`.

        This is a private function and should not be called from outside the
        MasterEquationSolver class! (Hence it does not perform further tests
        than MasterEquationSolver.)

        :param Γ_L: Matrix containing the transition rates regarding electrons
                    tunneling trough the *left* barrier where
                    :code:`Γ_L[i,j]` ≣ :math:`\Gamma^L_{i \\rightarrow j}`.
                    Must be a ``nxn`` matrix.
        :type Γ_L: numpy.ndarray
        :param Γ_R: Matrix containing the transition rates regarding electrons
                    tunneling trough the *right* barrier where
                    :code:`Γ_R[i,j]` ≣ :math:`\Gamma^R_{i \\rightarrow j}`.
                    Must be a ``nxn`` matrix.
        :type Γ_R: numpy.ndarray
        :param ns: List of numbers of states :math:`n_i` with the same number
                   :math:`i` of electrons in inreasing order of :math:`i`:

                   :code:`ns = [n_0, n_1, …, n_k]`
                   ≣ :math:`(n_0, n_1, \\dots, n_k)`.

                   E.g.
                   :code:`ns = [<numer of empty states>, <number of one body states>, ...]`
        :type ns: list(int)

        :return: :code:`stat_curs` where :code:`stat_curs[i][k]` ≣
                 :math:`I^k_{stat,i}`
        :rtype: numpy.ndarray(float)

        See documentation of MasterEquationSolver.doCalculation (2nd problem)
        for more details.
        """

        # Construct a Matrix σ containing the signs representing the direction
        # of current for each rate.
        σ = self.__calculateSigma(ns, verbose)
        # Calculate σΓ_k - a direction sensivtive form of the rates matricies.
        σΓ_L =  np.multiply(σ,Γ_L)
        σΓ_R = -np.multiply(σ,Γ_R)
        # verbose only: Print all σΓ.
        if verbose: print("σΓ_L: \n", σΓ_L)
        if verbose: print("σΓ_R: \n", σΓ_R)

        # Calculate the stationary current(s).
        # If no stationary solutions could be calculated there are no currents
        # to calculate.
        if stat_ps.shape[1] > 0:
            stat_curs = (np.sum(np.dot(σΓ_L, stat_ps), axis=0),
                         np.sum(np.dot(σΓ_R, stat_ps), axis=0))
        else:
            stat_curs = []

        return np.array(stat_curs)


    def __simulateDynamicCurrent(self,
                                Γ_L, Γ_R,
                                ns,
                                sim_time_dev_prop,
                                verbose
                               ):
        """
        Simulates the time current :math:`I(t)` for an :math:`n`-state
        system with transition rates :math:`\\Gamma^i`.

        This is a private function and should not be called from outside the
        MasterEquationSolver class! (Hence it does not perform further tests
        than MasterEquationSolver.)

        :param Γ_L: Matrix containing the transition rates regarding electrons
                    tunneling trough the *left* barrier where
                    :code:`Γ_L[i,j]` ≣ :math:`\Gamma^L_{i \\rightarrow j}`.
                    Must be a ``nxn`` matrix.
        :type Γ_L: numpy.ndarray
        :param Γ_R: Matrix containing the transition rates regarding electrons
                    tunneling trough the *right* barrier where
                    :code:`Γ_R[i,j]` ≣ :math:`\Gamma^R_{i \\rightarrow j}`.
                    Must be a ``nxn`` matrix.
        :type Γ_R: numpy.ndarray
        :param ns: List of numbers of states :math:`n_i` with the same number
                   :math:`i` of electrons in inreasing order of :math:`i`:

                   :code:`ns = [n_0, n_1, …, n_k]`
                   ≣ :math:`(n_0, n_1, \\dots, n_k)`.

                   E.g.
                   :code:`ns = [<numer of empty states>, <number of one body states>, ...]`
        :type ns: list(int)

        :return: Returns a Simulation of :math:`I(t)`
        :rtype: HaPPPy.MasterEquation.Simulation

        See documentation of MasterEquationSolver.doCalculation (2nd problem)
        for more details.
        """

        # The validity is inherited from the propabilities simulation.
        valid = sim_time_dev_prop.valid()
        # Copy all parameters from the propabilities simulation.
        sim = Simulation(sim_time_dev_prop.getΔt(),sim_time_dev_prop.getT_max())

        # Construct a Matrix σ containing the signs representing the direction
        # of current for each rate.
        σ = self.__calculateSigma(ns, verbose)
        # Calculate σΓ_k - a direction sensivtive form of the rates matricies.
        σΓ_L =  np.multiply(σ,Γ_L)
        σΓ_R = -np.multiply(σ,Γ_R)
        # verbose only: Print all σΓ.
        if verbose: print("σΓ_L: \n", σΓ_L)
        if verbose: print("σΓ_R: \n", σΓ_R)

        # For each existing discrete solution of P(t) calculate I_k(t).
        # (See the documentation of this module for more deails about this
        # calculation.)
        for i in range(len(sim_time_dev_prop)):
            # Get the current propability distribution as dotable vector.
            P_i = np.matrix(sim_time_dev_prop[i]).transpose()
            # verbose only: Print P(t).
            if verbose: print("P[" + str(i) + "] =\n = ", P_i)
            # Calculate all I(t) like in the formula in the documentation.
            I_L_i = (sum(np.dot(σΓ_L, P_i))).item(0,0)
            I_R_i = (sum(np.dot(σΓ_R, P_i))).item(0,0)
            # verbose only: Print all I_k(t).
            if verbose: print("I_L_t[" + str(i) + "]= \n", I_L_i)
            if verbose: print("I_R_t[" + str(i) + "]= \n", I_R_i)
            # Extend the simulation of I_k(t) by this new discrete values.
            sim.append(np.array([I_L_i, I_R_i]))

        # Make the simulation valid - iff no issues occurred.
        if valid:
            sim.validate()

        # Return the (valid or invaid) simulation.
        return sim

class Simulation():
    """
    This class represents a simulation crated by
    HaPPPy.MasterEquation.MasterEquationSolver.

    Let :math:`X` be an arbitrary set. A function
    :math:`f : [0, t_{max}] \\to X, t \mapsto f(t)` is approximated by
    calculating its discrete values :math:`f(n \\cdot \\Delta t)`
    in the time intervall :math:`0 \\leq t = n \\cdot \\Delta t \\leq t_{max}`
    where :math:`n \in \mathbb{N}_0, \Delta t \in \mathbb{R}^+` and
    :math:`t_{max} \in \mathbb{R}^+_0`.

    If :code:`sim` is a valid HaPPPy.MasterEquation.Simulation then:

    :code:`sim.getParameters()` ≣ :math:`(0, \\Delta t, 2 \\cdot \\Delta t,
    \\dots, m \\cdot \\Delta t)`

    and

    :code:`sim.getValues()` ≣ :math:`(f(0), f(\\Delta t), f(2 \\cdot \\Delta t),
    \\dots, f(m \\cdot \\Delta t))`.

    where :math:`m = \\lfloor \\frac{t_{max}}{\\Delta t} \\rfloor`
    ≣ :code:`len(sim) - 1`.

    Note: :math:`m` might be smaller for invalid simulations (see valid()).
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

    def getBins(self):
        """
        :return: A list of all values :math:`n \\cdot \\Delta t` where
                 :math:`n \in \mathbb{N}_0 \land n \cdot \Delta t \leq t_{max}`
                 with increasing :math:`n`.
        :rtype: numpy.ndarray
        """
        return np.arange(0, self.t_max + self.Δt, self.Δt)

    def getParameters(self, start=0, stop=None, step=1):
        """
        :return: A list of all **valid** values :math:`n \\cdot \\Delta t` where
                 :math:`n \in \mathbb{N}_0 \land n \cdot \Delta t \leq t_{max}`
                 with increasing :math:`n`.
        :rtype: numpy.ndarray
        """
        if stop == None or stop > len(self):
            stop = len(self)

        return self.Δt * np.arange(start, stop, step)

    def getValues(self, start=0, stop=None, step=1):
        """
        :return: A list of all **valid** values :math:`f(n \\cdot \\Delta t)`
                 where :math:`f` is the simulated function and :math:`n \in
                 \mathbb{N}_0 \land n \cdot \Delta t \leq t_{max}` with
                 increasing :math:`n`.
        :rtype: numpy.ndarray
        """
        if stop == None or stop > len(self):
            stop = len(self)

        return np.array([self.__values[i] for i in range(start, stop, step)])

    def append(self, value):
        # for internal use only
        self.__values.append(value)

    def valid(self):
        """
        :return: True iff simulation was not aborted.
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
        ts = self.getParameters()
        vs = self.getValues()
        return ("Sim<valid="
                + str(self.valid()) + ","
                + str([(ts[i], vs[i]) for i in range(len(self))])
                + ">"
               )

    def __str__(self):
        return repr(self)

    def quickPlot(self,
                  title=None,
                  xlabel=None, ylabel=None,
                  xunit=1, yunit=1,
                  xunit_symbol=None, yunit_symbol=None,
                  legend=None,
                  start=0,
                  stop=None,
                  step=1,
                 ):
        """
        Simple plotting method to quickly get an overview on the simulation.

        :param title: The title of the plot. (optional)
        :type title: string
        :param xlabel: The symbol to retresent the parameter. (optional)
        :type xlabel: string
        :param ylabel: The symbol to retresent the function values. (optional)
                       If the function values are vetors ylabel is treated as a
                       LATEX expression representing a symbol and automatic
                       indicies are added.
        :type ylabel: string
        :param xunit: The unit the paramter is meassured in. All x-values are
                      divided by this value before plotting. (optional)
        :type xunit: float
        :param yunit: The unit the function valuesare meassured in. All y-values
                      are divided by this value before plotting. (optional)
        :type yunit: float
        :param xunit_symbol: The unit the paramter is meassured in. (optional)
        :type xunit_symbol: string
        :param yunit_symbol: The unit the function valuesare meassured in. (optional)
        :type yunit_symbol: string
        :param legend: A list of strings representing the symbols of each graph.
                       Vectorial graphs are ploted component-wise and labeled in
                       increasing order of their index. (optional)
        :type legend: list(string)
        :param start: The :math:`i`th value to start the plot with. (optional)
        :type start: int
        :param stop: The :math:`i`th value to stop the plot with.
                     :math:`i` is not included. (optional)
        :type stop: int
        :param start: The distance between two ploted bins. E.g. :code:`step=1`
                      plots all values and :code:`i=2` plots every odd indexed
                      value. (optional)
        :type step: int

        :example: See the example given at the documentation of
                  HaPPPy.MasterEquation.MasterEquationSolver.
        """
        # Aquire all values related to the discrete representation of the
        # function.
        ts = self.getParameters(start, stop, step)
        vs = self.getValues(start, stop, step)

        # Plot the function.
        plt.plot(ts / xunit, vs / yunit)
        # Add a validity mark to the title if necessary.
        if not self.valid():
            if title == None:
                title = " (not valid)"
            else:
                title = str(title) + " (not valid)"
        # Add a title if requested.
        if title != None:
            plt.title(str(title))
        # Add labels to one or both axes if requested.
        # It is possible to add an optional unit to each axis.
        if xlabel != None:
            if xunit_symbol != None:
                plt.xlabel("$" + str(xlabel) + " / " + str(xunit_symbol) + "$")
            else:
                plt.xlabel("$" + str(xlabel) + "$")
        if ylabel != None:
            if yunit_symbol != None:
                plt.ylabel("$" + str(ylabel) + " / " + str(yunit_symbol) + "$")
            else:
                plt.ylabel("$" + str(ylabel) + "$")
        # Add a legend if requested.
        # If the function values are vetors ylabel is treated as a LATEX
        # expression representing a symbol and automatic indicies are added.
        if (legend == None
            and ylabel != None
            and (type(vs[0]) == list or type(vs[0]) == np.ndarray)
            and len(vs[0]) > 1
           ):
            n = len(vs[0]) # dimension of v
            legend = ["${" + str(ylabel) + "}_{" + str(i) + "}$"
                      for i in range(n)
                     ]
        if legend != None:
            plt.legend(legend)
        # Add a grid.
        plt.grid()
        # Draw the figure.
        plt.show()
