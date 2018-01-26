from .TwoParticle import createTwoParticleData
from .OneParticleLoader import SpectrumData
from .TwoParticleLoader import TwoBodySpectrumData

class TwoBodySolver:
    """
    
    
    
    """


    def __init__(self):
        """ The constructor.
        """

        print("Hello from the TwoBodySolver, what can I do for you?")

    def doCalculation(
            self, obDataFile='HaPPPy/OneBody/data_group1', 
            tbDataFile='HaPPPy/TwoBody/data_group2'):
        """
        
        Calculates and returns the two-electron eigenvalues and 
        eigenvector coefficient matrices from the single-electron 
        eigenfunctions.
        
        :param obDataFile: Path to the file containing the one-body
         data to use as input (filename without '.hdf5' ending)
        :type obDataFile: str
        :param tbDataFile: Path for exporting the resulting two-body
         data (filename without '.hdf5' ending)
        :type tbDataFile: str
        :return: Returns :code:`[tb_energ,tb_coeff]`, where
         :code:`tb_energ` is an array containing the energy eigenvalues
         of the two-particle Hamiltonian in ascending order and
         :code:`tb_coeff[i,j,k]` ≡ :math:`Q^i_{j,k}`, with
         :code:`tb_coeff[i,:,:]` being the eigenstate belonging to
         :code:`tb_energ[i]`.
        :rtype: numpy.ndarray(numpy.float64),
         numpy.ndarray(numpy.float64)
        """
        obData = SpectrumData()
        obData.open(obDataFile)
        result = createTwoParticleData(obData)
        tb_energ = result[0]
        tb_coeff = result[1]
        data = TwoBodySpectrumData()
        data.init(tbDataFile, len(tb_energ), tb_coeff.shape[1])
        data.energies[:] = tb_energ
        data.coeficients[:,:,:] = tb_coeff
        data.close()
        return result

    def getSavedData(self, tbDataFile='HaPPPy/TwoBody/data_group2'):
        """
        
        Reads the previously generated two-electron eigenvalues and 
        eigenvector coefficient matrices from a file compatible with 
        the TwoBodySpectrumData class.
        
        :param tbDataFile: Path to the data file.
        :type tbDataFile: str
        :return: Returns :code:`[tb_energ,tb_coeff]`, where
         :code:`tb_energ` is an array containing the energy eigenvalues
         of the two-particle Hamiltonian in ascending order and
         :code:`tb_coeff[i,j,k]` ≡ :math:`Q^i_{j,k}`.
        :rtype: numpy.ndarray(numpy.float64),
         numpy.ndarray(numpy.float64)
        
        """
        tbData = TwoBodySpectrumData()
        tbData.open(tbDataFile)
        tb_energ = tbData.energies[:]
        tb_coeff = tbData.coeficients[:,:,:]
        return [tb_energ, tb_coeff]
