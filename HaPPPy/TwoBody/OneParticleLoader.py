import h5py
import numpy as np
import math

#private constants holding the dataset names
_EN_NAME = "eigenvalues"
_EV_NAME = "eigenvectors"
_OPT_NAME = "settings"

class SpectrumData:
	"""Used to store and load data of one-particle wave energy spectrums in hdf5-files.
	
	Fields available after initializing with 'open()' or 'init()':
	file -- reference to the hdf5 file handle
	energies -- dataset containing the energy eigenvalues in [meV] as floating point numbers: energies[i] stores the i-th eigenvalue
	waves -- dataset containing the wavefunctions as rasterized floating point arrays of probabilities in [nm^(-1/2)]: waves[i,j] stores the i-th grid entry of the j-th eigenvector
	m -- the number of eigenvalues/-functions
	n -- the number of grid points used to represent the wavefunctions
	dx -- the spacial distance between grid points in [nm]: (x1-x0) = dx * (n-1)
	l -- the spacial width of the grid in [nm]: l = dx * n
	"""
	
	def open(self, filename):
		"""Initialize the data by loading a hdf5-file
		
		Arguments:
		filename -- path to the file containing the data (without the '.hsf5' ending)
		"""
		self.file = h5py.File(filename + ".hdf5", "a")

		self.energies = self.file[_EN_NAME]
		self.waves = self.file[_EV_NAME]
		par = self.file[_OPT_NAME]

		self.n = len(self.waves[:,0])
		self.m = len(self.energies)
		self.l = float(par[1][1])

		self.dx = self.l / self.n
		
	def init(self, filename, m, n, dx=None, L=None, info=None):
		"""Initialize with empty data and create the datasets for a new hdf5-file
		
		Arguments:
		filename -- path to the file to store the data in (without the '.hsf5' ending)
		m -- the number of eigenvalues/-functions
		n -- the number of grid points used to represent the wavefunctions
		dx -- the spacial distance between grid points in [nm] (optional): (n-1) * dx = x1-x0
		L -- the spacial width of the whole grid in [nm] (optional): L = n * dx
		info -- a 2D string array containing arbitrary configuration data in the form info[i,0] = key, info[i,1] = value (optional): usually info[0,1] strores n and info[1,1] stores L
		of the arguments dx and L only one must be given (alternatively L could also be given via info) otherwise a ValueError is raised!
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
		if info is None:
			info = np.array([["n-grids point" ,str(self.n)],["l-lenght of potential",str(self.l)]]).astype('S9')
		self.file.create_dataset(_OPT_NAME, data=info)
		self.file.flush()

	def close(self):
		"""Close the underlying hdf5-file and save all modified data
		
		After this the datasets won't be accessible anymore!
		"""
		self.file.close()
	
	def getNormalizedWaves(self):
		""""Helper method to get the eigenvectors normalized (in case the data is potentially not nomalized yet)
		
		return a 2D numpy array containing all eigenvectors (v) spatially normalized: sum x in v of (|x|^2 * self.dx) = 1.0
		"""
		waves = np.empty((self.n, self.m))
		for i in range(self.m):
			vec = self.waves[:,i]
			waves[:,i] = vec / math.sqrt(np.inner(vec, vec) * self.dx)
		return waves

