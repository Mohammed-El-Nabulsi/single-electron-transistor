#Input: energy in meV. internal calculations in meV and 10 femtoseconds.
import numpy as np

from scipy.constants import codata

from A_Potential01 import Potential
from B_Gauss00 import GaussianWave
from C_Fourier03 import Fourier

me   = codata.value("electron mass energy equivalent in MeV") * 1e8 ;  # Convert to milli eV/c^2
hbar = codata.value("Planck constant over 2 pi in eV s")      * 1e19;  # Convert to 10femtoseconds*millielectronvolts

a = np.sqrt(me)/hbar

#new units
sqrt_two_me_over_hbar = np.sqrt(2*me)/hbar
print(sqrt_two_me_over_hbar)

#Main

#set testparameters
barriere = np.array([2,2,2,2,2])
dy = 1/100 #Abstand zwischen zwei sample points
energy = 1

#testing A_Potential00

potential = Potential(barriere, dy)


#print(potential.dx)
#print(potential.barrier)
#print(potential.gauss_index_width)
gauss_width = potential.gauss_width
N = potential.potential.size
#print(N)
#print(N)
#print()
#print(potential.potential[:10])
#print()
positions = potential.position_grid
#print(potential.pos_grid_width)
#print(positions[:10])
#print(potential.gauss_symmetry_index)
gauss_symmetry_point = potential.gauss_symmerey_point
#print(gauss_symmetry_point)


#finished, everything works fine here.

#testing interaction between gauss and potential class

gauss = GaussianWave(positions,gauss_symmetry_point,gauss_width,energy)

#print(gauss.width)
#print(gauss.symm)
print(gauss.k0)
#print(gauss.x_grid.size)
#print(gauss.x_probability)
#gauss.plot_x_package()
#gauss_package = gauss.x_package
#print(gauss_package)

#interaction test finished, everything works fine here.

#testing Fourier
fourier = Fourier(N,dy)
#print(fourier.get_wavenumbers().size)
fourier.plot_dft(gauss.x_package)

#Anm.: Fourier Trafo des Gausspaktes ist anscheiend wieder Gauss, aber k-werte und amplithuden stimmen keinesfalls!
