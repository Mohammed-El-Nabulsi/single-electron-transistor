import h5py

#private constants holding the dataset names
_EN_NAME = "eigenvalues"
_EV_NAME = "eigenvectors"
_OPT_NAME = "settings"

class SpectrumData:
	"""Used to store and load data of one-particle wave energy spectrums in hdf5-files.
	
	Fields available after initializing with 'open()' or 'init()':
	file -- reference to the hdf5 file handle
	energies -- dataset containing the energy eigenvalues in [meV] as floating point numbers: energies[i] stores the i-th eigenvalue
	waves -- dataset containing the wavefunctions as rasterized floating point arrays of probabilities in [1/nm]: waves[i,j] stores the j-th grid entry of the i-th eigenvector
	m -- the number of eigenvalues/-functions
	n -- the number of grid points used to represent the wavefunctions
	dx -- the spacial distance between grid points in [nm]: (x1-x0) = dx * (n-1)
	x0 -- the spacial location of the 0-th grid point in [nm]
	x1 -- the spacial location of the (n-1)-th grid point in [nm]
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
		self.x0 = par[0]
		self.x1 = par[1]
		self.n = len(self.waves[0,:])
		self.m = len(self.energies)
		self.dx = (self.x1 - self.x0) / (self.n - 1)
	
	def init(self, filename, m, n, x0=None, x1=None, dx=None, L=None):
		"""Initialize with empty data and create the datasets for a new hdf5-file
		
		Arguments:
		filename -- path to the file to store the data in (without the '.hsf5' ending)
		m -- the number of eigenvalues/-functions
		n -- the number of grid points used to represent the wavefunctions
		x0 -- the spacial location of the 0-th grid point in [nm] (optional)
		x1 -- the spacial location of the (n-1)-th grid point in [nm] (optional)
		dx -- the spacial distance between grid points in [nm] (optional): (n-1) * dx = x1-x0
		L -- the spacial width of the whole grid in [nm] (optional): L = n * dx
		of the arguments x0, x1, dx and L either dx, L or both x0 and x1 must be given otherwise a ValueError is raised!
		"""
		if L:
			dx = L / n
		if dx:
			if x0:
				x1 = x0 + (n-1) * dx
			elif x1:
				x0 = x1 - (n-1) * dx
			else:
				x1 = (n-1) * dx / 2.0
				x0 = -x1
		elif x0 and x1:
			dx = (x1-x0) / (n-1)
		else:
			raise ValueError("Insuficient grid parameters given!")
		self.m = m
		self.n = n
		self.x0 = x0
		self.x1 = x1
		self.dx = dx
		self.file = h5py.File(filename + ".hdf5", "w")
		self.energies = self.file.create_dataset(_EN_NAME, (m,))
		self.waves = self.file.create_dataset(_EV_NAME, (m, n))
		par = self.file.create_dataset(_OPT_NAME, (2,))
		par[0] = x0
		par[1] = x1
		self.file.flush()

	def close(self):
		"""Close the underlying hdf5-file and save all modified data
		
		After this the datasets won't be accessible anymore!
		"""
		self.file.close()

