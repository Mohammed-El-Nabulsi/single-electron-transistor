import h5py
import numpy as np
import math

#private constants holding the dataset names
_EN_NAME = "eigenvalues"
_EV_NAME = "eigenvectors"
_OPT_NAME = "settings"
_POT_NAME = "potential"

class SpectrumData:
    """
    
    Used to store and load data of the one-particle problem in 
    hdf5-files.
    
    Fields available after initializing with :meth:`open` or 
    :meth:`init`\:
    
    - :code:`file` -- reference to the hdf5 file handle
    - :code:`energies` -- dataset containing the energy eigenvalues in 
      [meV] as floating point numbers: energies[i] stores the i-th 
      eigenvalue
    - :code:`waves` -- dataset containing the wavefunctions as 
      rasterized floating point arrays of probabilities in [nm^(-1/2)]:
      waves[i,j] stores the i-th grid entry of the j-th eigenvector
    - :code:`potential` -- dataset containing the potential that caused
      the one body states in [meV]: potential[i] stores potential
      energy at grid index i
    - :code:`m` -- the number of eigenvalues/-functions
    - :code:`n` -- the number of grid points used to represent the
      wavefunctions
    - :code:`dx` -- the spacial distance between grid points in [nm]: 
      (x1-x0) = dx * (n-1)
    - :code:`l` -- the spacial width of the grid in [nm]: 
      :code:`l = dx * n`
    
    Preconditions to data:
    
    - all the rasterized wavefunctions :code:`v` should be spatially 
      normalized: :code:`sum x in v of (|x|^2 * dx) = 1.0`
    - all the one-particle states should be sorted by energy in 
      ascending order: :code:`energies[i+1] >= energies[i]`

    """
    
    def open(self, filename):
        """
        
        Initialize the data by loading a hdf5-file
        
        :param filename: Path to the file containing the data (without
         the '.hsf5' ending).
        :type filename: str

        """
        self.file = h5py.File(filename + ".hdf5", "a")

        self.energies = self.file[_EN_NAME]
        self.waves = self.file[_EV_NAME]
        self.potential = self.file[_POT_NAME]
        par = self.file[_OPT_NAME]

        self.n = len(self.waves[:,0])
        self.m = len(self.energies)
        self.l = float(par[1][1])

        self.dx = self.l / self.n
        
    def init(self, filename, m, n, dx=None, L=None, info=None):
        """

        Initialize with empty data and create the datasets for a new 
        hdf5-file
        
        :param filename: Path to the file to store the data in (without
         the '.hsf5' ending).
        :type filename: str
        :param m: The number of eigenvalues/-functions.
        :type m: int
        :param n: The number of grid points used to represent the
         wavefunctions.
        :type n: int
        :param dx: The spacial distance between grid points in [nm]
         (optional): (n-1) * dx = x1-x0.
        :type dx: float
        :param L: The spacial width of the whole grid in [nm]
         (optional): L = n * dx
        :type L: float
        :param info: A 2D string array containing arbitrary
         configuration data in the form info[i,0] = key, info[i,1] =
         value (optional): usually info[0,1] strores n and info[1,1]
         stores L 
        :type info: numpy.ndarray
        Of the arguments **dx** and **L** only one must be given 
        (alternatively **L** could also be given via info); otherwise a
        ValueError is raised!

        """
        if info is not None:
            L = float(info[1][1])
        if L:
            dx = L / n
        elif dx:
            L = dx * n
        else:
            raise ValueError("Insuficient grid parameters given!")
        self.m = m
        self.n = n
        self.dx = dx
        self.l = L
        self.file = h5py.File(filename + ".hdf5", "w")
        self.energies = self.file.create_dataset(_EN_NAME, (m,), dtype='d')
        self.waves = self.file.create_dataset(_EV_NAME, (n, m), dtype='d')
        self.potential = self.file.create_dataset(_POT_NAME, (n,), dtype='d')
        if info is None:
            info = np.array(
                [["n-grids point" ,str(self.n)],
                ["l-lenght of potential",str(self.l)]]).astype('S9')
        self.file.create_dataset(_OPT_NAME, data=info)
        self.file.flush()

    def close(self):
        """

        Close the underlying hdf5-file and save all modified data.
        
        After this the datasets won't be accessible!
        
        """
        self.file.close()
    
    def getNormalizedWaves(self):
        """

        Helper method to get the eigenvectors normalized (in case the 
        data is potentially not nomalized yet)
        
        :return: a 2D numpy array containing all eigenvectors :code:`v`
         spatially normalized: :code:`sum x in v of (|x|^2 * self.dx) = 1.0`
        :rtype: numpy.ndarray

        """
        waves = np.empty((self.n, self.m))
        for i in range(self.m):
            vec = self.waves[:,i]
            waves[:,i] = vec / math.sqrt(np.inner(vec, vec) * self.dx)
        return waves

