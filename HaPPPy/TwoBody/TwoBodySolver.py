from .TwoParticle import createTwoParticleData

class TwoBodySolver:
    """ Solves the two body problem

    TODO: Add more documentation

    """


    def __init__(self):
        """ The constructor.
        """

        print("Hello from the TwoBodySolver, what can I do for you?")

    def doCalculation(self):
	""" Calculate and return the two-electron eigenvalues and eigenvectors from the single-electron eigenfunctions. Returns [E,Q], where E is an array of (scalar) eigen values and Q is an array of matrices of the shape [i,j,n], where Q[i,j,n] is the coefficient of the nth eigenvector belonging to the |i,j> product basis function. These arrays do not have a special ordering, but Q[:,:,n] is the eigenvector belonging to E[n]."""
        return createTwoParticleData()
