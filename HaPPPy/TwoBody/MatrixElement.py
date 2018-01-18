import numpy as np
#import matplotlib.pyplot as plt
import math

from numpy import sin, exp

_epsilon = 1e-1
_coulomb_consts = 1.440343e3 #[nm*meV]

def getMatrixElement(dX, n, Phi1, Phi2, Phi3, Phi4):
	"""Compute the coulomb energy matrix element for two pairs of single-electron wavefunctions.
	
	Arguments:
	dX -- spacial distance between the grid points used to model the wavefunctions in [nm]
	n -- number of grid points used to model the wavefunctions
	Phi1 -- wavefunction of electron A from first entry
	Phi2 -- wavefunction of electron B from first entry
	Phi3 -- wavefunction of electron A from second entry
	Phi4 -- wavefunction of electron B from second entry
	Return energy of coulomb interaction in [meV]
	"""
	sum = 0.0
	for i1 in range(n):
		for i2 in range(n):
			if i1 == i2:
				r = _epsilon * dX
			else:
				r = abs(i1 - i2) * dX
			sum += Phi1[i1] * Phi2[i2] / r * Phi3[i1] * Phi4[i2]
	return sum * dX * dX * _coulomb_consts

#### testing ####

def testFunction(x0, x1, n, m):
	L = x1 - x0
	y = np.zeros(m)
	t0 = np.sqrt(2 / L)
	for i in range(m):
		x = L * (float(i) + 0.5) / float(m) + x0
		y[i] = t0 * np.sin(n * np.pi / L * x)
	return y

#x0 = -2.0
#x1 = 2.0
#Phi1 = testFunction(x0, x1, 1)
#Phi2 = testFunction(x0, x1, 2)
#Phi3 = testFunction(x0, x1, 2)
#Phi4 = testFunction(x0, x1, 5)
#print(getMatrixElement(Phi1, Phi2, Phi3, Phi4))
