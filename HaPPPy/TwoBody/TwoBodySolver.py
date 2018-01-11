from .TwoParticle import createTwoParticleData
from .OneParticleLoader import SpectrumData

class TwoBodySolver:
    """ Solves the two body problem

    TODO: Add more documentation

    """


    def __init__(self):
        """ The constructor.
        """

        print("Hello from the TwoBodySolver, what can I do for you?")

    def doCalculation(self, obDataFile='HaPPPy/OneBody/data_group1'):
        """ Calculate and return the two-electron eigenvalues and eigenvectors from the single-electron eigenfunctions.
        
        Arguments:
        obDataFile -- path to the file containing the one body data to use as input (filename without '.hdf5' ending)
        
        Returns [E,Q], where E is an array of (scalar) eigen values and Q is an array of matrices of the shape [i,j,n], where Q[i,j,n] is the coefficient of the nth eigenvector belonging to the |i,j> product basis function. These arrays do not have a special ordering, but Q[:,:,n] is the eigenvector belonging to E[n].
        """
        obData = SpectrumData()
        obData.open(obDataFile)
        return createTwoParticleData(obData)
