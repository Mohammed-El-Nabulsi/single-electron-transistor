import numpy as np
import numpy.linalg as lin


__docformat__ = 'reStructuredText'

class MasterEquationSolver:
    """ Solves the master equation

    TODO: Add more documentation

    """


    def __init__(self):
        """ The constructor.
        """

        print("Hello from the MasterEquationSolver, what can I do for you?")

    def doCalculation(self):
        """ Some dummy calculation
        
        Returns:
            double.  The result

        """

        return 1.0+1.0

    def get_span_of_stable_solution(Γ, ε=1E-10, verbose=False):
        """
        Finds span of stable solutions of the master equation.

        This method returns a numerical approximation of the solution to the following problem.

        Let

        .. math::

            n \in \mathbb{N}, \Gamma \in Mat(n, n, \mathbb{R}^+_0), \\vec{P} \in \mathbb{R} \\to { \left ( \mathbb{R}^+_0 \\right ) }^n, t \mapsto \\vec{P}(t)

        The master eqation (a differential equation of first order with constant coefficients) is stated as

        .. math::

            \\frac{d P_{\\alpha}}{d t} = \sum_{\\beta} \Gamma_{\\beta \\rightarrow \\alpha} P_{\\beta} - \sum_{\\beta} \Gamma_{\\alpha \\rightarrow \\beta} P_{\\alpha} \\\\ \\alpha, \\beta \in \{ 1,2, \dots , n \}

        where

        .. math::

            \Gamma_{\\alpha \\rightarrow \\beta} \equiv \Gamma_{\\alpha \\beta}

        denotes the rate from state :math:`\\alpha` to :math:`\\beta`. Stable solutions satisfy :math:`\\frac{d P_{\\alpha}}{d t} = 0, \\alpha \in \{ 1,2, \dots , n \}` and hence can be abbreviated to a single vector :math:`\\vec{P} \in { \left ( \mathbb{R}^+_0 \\right ) }^n`.

        The returned value is a list of vectors which span the space of possible stable solutions when normalized with respect to :math:`\sum_{i=1}^{n} P_i = 1`.

        :param Γ: Matrix containing the transition rates where :code:`Γ[i][j]` ≣ :math:`\Gamma_{i \\rightarrow j}`. Must be a ``nxn`` matrix.
        :type Γ: numpy.array
        :param ε: Tolerance value. Due to numerical approximations during calculations a tolerance value is needed.
        :type ε: float
        :param verbose: Enables the verbose mode: Calculation will be printed in detail to support debugging.
        :type verbose: bool

        :return: Returns a list of vectors which span the space of possible stable solution to the master eqation when normalized with respect to :math:`\sum_{i=1}^{n} P_i = 1`. Each vector of the list satisfies the mentioned normalization.

        :example: .. code-block:: python

                import numpy as np
                from HaPPPy.MasterEquation.MasterEquationSolver import MasterEquationSolver as MES

                ## test program for get_span_of_stable_solution
                # set-up reasonable Γ-matrices
                Γ_1 = np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
                Γ_2 = np.array([[0, 0.0, 0.0], [0.0, 0, 0.0], [0.0, 0.0, 0]])
                # calculate
                print(MES.get_span_of_stable_solution(Γ_1))
                print(MES.get_span_of_stable_solution(Γ_2))

                # Output:
                # [array([ 0.33333333,  0.33333333,  0.33333333])]
                # [array([ 1.,  0.,  0.]), array([ 0.,  1.,  0.]), array([ 0.,  0.,  1.])]

            .

        TODO: Extend description of ε and add more details about the calculation.
              Alternative: Narrow properties of Γ: Γ must have exactly one eigenvalue equals to 0. It follows that there would be only one stable solution!

        """

        ## input checks
        # warn if tolerance value is unreasonable
        if ε <= 0:
            raise RuntimeError("ε must be a positive number > 0. \nε = " + str(ε))
        # Γ must be a nxn matrix
        if Γ.ndim != 2 or Γ.shape[0] != Γ.shape[1] or not (Γ >= 0).all():
            raise RuntimeError("Γ must be a square matrix with coefficients >= 0. \nΓ = \n" + str(Γ))

        if verbose:
            print("Γ = \n", Γ)
            print("ε = ", ε)

        ## calculation

        # set-up Λ-matrix
        Λ_in = Γ
        Λ_out = np.diag([sum([Γ[i][j] for j in range(Γ.shape[0])]) for i in range(Γ.shape[0])])
        Λ = Λ_in - Λ_out
        if verbose: print("Λ_in = \n", Λ_in)
        if verbose: print("Λ_out = \n", Λ_out)
        if verbose: print("Λ = \n", Λ)

        # find eigenvector basis of Λ (Λ_evecs_T transforms to eigenvevtor basis of Λ)
        (Λ_evals, Λ_evecs) = lin.eigh(Λ)
        Λ_evecs_T = Λ_evecs.transpose()
        if verbose: print("Λ.eigenvalues = \n", Λ_evals)
        if verbose: print("Λ.eigenvectors = \n", Λ_evecs_T)

        # get Λ eigenvector basis of Λ
        Λ_evb = np.diag(Λ_evals)
        if verbose: print("Λ in Λ.eingenvetorbase = \n", Λ_evb)

        # get inverse of Λ_evecs_T (transforms back from eigenvevtor basis of Λ)
        Λ_evecs_T_inv = np.linalg.inv(Λ_evecs_T)
        if verbose: print("Λ.eingenvectors^-1 = \n", Λ_evecs_T_inv)

        # check if creation of inversion was successful
        if (np.dot(Λ_evecs_T_inv, Λ_evecs_T) != np.identity(Λ.shape[0])).all():
            raise RuntimeError("Can not invert Λ = \n" + str(Λ))

        # find eigenvalues with eigenvector which do not diverge or become 0 when t --> ∞
        Λ_evals_zero_indices = []
        for i in range(len(Λ_evals)):
            if abs(Λ_evals[i]) <= ε:
                Λ_evals_zero_indices.append(i)
        if verbose:
            print("indices of eigenvetors with eigenvalue 0 (within tolerance) = ", Λ_evals_zero_indices)

        # get the corresponding eigenvectors
        Λ_lasting_evecs = []
        for i in Λ_evals_zero_indices:
            Λ_lasting_evecs.append(Λ_evecs[:,i])
        if verbose:
            print("eigenvectors with eigenvalue 0 (within tolerance) = \n", Λ_lasting_evecs)

        # renoralize such that sum(P_i) = 1
        for i in range(len(Λ_lasting_evecs)):
            Λ_lasting_evecs[i] = (1 / sum(Λ_lasting_evecs[i])) * Λ_lasting_evecs[i]
        if verbose:
            print("normalized eigenvectors with eigenvalue 0 (within tolerance) = \n", Λ_lasting_evecs)

        return Λ_lasting_evecs
