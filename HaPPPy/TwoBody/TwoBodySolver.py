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

    def doCalculation(self, obDataFile='HaPPPy/OneBody/data_group1', tbDataFile='HaPPPy/TwoBody/data_group2'):
        """ Calculate and return the two-electron eigenvalues and eigenvectors from the single-electron eigenfunctions.
        
        Arguments:
        obDataFile -- path to the file containing the one body data to use as input (filename without '.hdf5' ending)
        tbDataFile -- path to the file exporting the resulting two body data (filename without '.hdf5' ending)
        
        Returns [E,Q], where E is an array of (scalar) eigen values and Q is an array of matrices of the shape [n,i,j], where Q[n,i,j] is the coefficient of the nth eigenvector belonging to the |i,j> product basis function. hese arrays are sorted by energy in ascending order, with Q[n,:,:] as the eigenvector belonging to E[n].
        """
        obData = SpectrumData()
        obData.open(obDataFile)
        result = createTwoParticleData(obData)
        ev = result[0]
        Q = result[1]
        data = TwoBodySpectrumData()
        data.init(tbDataFile, len(ev), Q.shape[1])
        data.energies[:] = ev
        data.coeficients[:,:,:] = Q
        data.close()
        return result

    def getSavedData(self):
        """ Returns the same as the last doCalculation(...,saveOutput=true) call without all of the time consuming calculation """
        #TODO
        return [[1,2,2,1],[[[1,1],[1,1]],[[2,2],[2,2]],[[1,2],[1,2]],[[2,1],[2,1]]]]
