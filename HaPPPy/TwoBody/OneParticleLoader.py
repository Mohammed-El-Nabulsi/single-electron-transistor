import h5py

_EN_NAME = "eigenvalues_group1"
_EV_NAME = "eigenvectors_group1"
_OPT_NAME = "settings"

class SpectrumData:

	def open(self, filename):
		self.file = h5py.File(filename + ".hdf5", "a")
		self.energies = file[_EN_NAME]
		self.waves = file[_EV_NAME]
		par = file[_OPT_NAME]
		self.x0 = par[0]
		self.x1 = par[1]
		self.n = len(waves[0,:])
		self.m = len(energies)
		self.dx = (x1 - x0) / (n - 1)
	
	# m: number of wave functions
	# n: number of raster indices used to represent wave functions
	# x0: x pos at first raster index
	# x1: x pos at last raster index
	def init(self, filename, m, n, x0, x1):
		self.m = m
		self.n = n
		self.x0 = x0
		self.x1 = x1
		self.dx = (x1 - x0) / (n - 1)
		self.file = h5py.File(filename + ".hdf5", "w")
		self.energies = self.file.create_dataset("energy", shape=(m))
		self.waves = self.file.create_dataset("waves", shape=(m, n))
		par = self.file.create_dataset("settings", shape=(2))
		par[0] = x0
		par[1] = x1
		self.file.flush()

	def close(self):
		file.close()

