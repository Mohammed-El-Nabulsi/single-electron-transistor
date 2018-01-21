import h5py
import numpy as np
import math

#private constants holding the dataset names
_EN_NAME = "eigenvalues"
_CM_NAME = "coefMat"
_OPT_NAME = "settings"

class TwoBodySpectrumData:
	"""Used to store and load data of two-particle wave energy spectrums in hdf5-files.
	
	Fields available after initializing with 'open()' or 'init()':
	file -- reference to the hdf5 file handle
	energies -- dataset containing the energy eigenvalues in [meV] as floating point numbers: energies[i] stores the i-th eigenvalue
	coeficients -- dataset containing the coeficient matrices: coeficients[i,k,l] contains the |k,l> product basis function entry in the coeficient matrix for the i-th eigenvalue
	m -- the number of eigenvalues & matrices
	n -- matrix size
	"""
	
	def open(self, filename):
		"""Initialize the data by loading a hdf5-file
		
		Arguments:
		filename -- path to the file containing the data (without the '.hsf5' ending)
		"""
		self.file = h5py.File(filename + ".hdf5", "a")

		self.energies = self.file[_EN_NAME]
		self.coeficients = self.file[_CM_NAME]
		self.par = self.file[_OPT_NAME]
		self.m = par[0]
		self.n = par[1]
		
	def init(self, filename, m, n):
		"""Initialize with empty data and create the datasets for a new hdf5-file
		
		Arguments:
		filename -- path to the file to store the data in (without the '.hsf5' ending)
		m -- the number of eigenvalues & matrices
		n -- matrix size
		"""
		self.m = m
		self.n = n
		self.file = h5py.File(filename + ".hdf5", "w")
		self.energies = self.file.create_dataset(_EN_NAME, (m,), dtype='d')
		self.coeficients = self.file.create_dataset(_CM_NAME, (n, m, m), dtype='d')
		self.file.create_dataset(_OPT_NAME, (2,), data=[m, n])
		self.file.flush()

	def close(self):
		"""Close the underlying hdf5-file and save all modified data
		
		After this the datasets won't be accessible anymore!
		"""
		self.file.close()

