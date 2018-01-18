import numpy as np
from scipy.constants import codata

import matplotlib.pyplot as plt

me   = codata.value("electron mass energy equivalent in MeV") * 1e8 ;  # Convert to milli eV/c^2
hbar = codata.value("Planck constant over 2 pi in eV s")      * 1e19;  # Convert to 10*fs*millielectronvolts



class GaussianWave():
    
    def __init__(self, x_grid, symm, width, energy):

        self.width = width
        self.symm = symm
        self.x_grid = x_grid
        
        self.k0 = np.sqrt(2*me*energy)/hbar
        
        self.x_package = self.create_gauss_x_package()
        
        self.x_probability = np.multiply(self.x_package, self.x_package.conj()).real
    
    def create_gauss_at_x(self,x):
        norm_x = ( 2/(np.pi*self.width**2) )**(1/4)
        return norm_x * np.exp(1j*self.k0*(x-self.symm)) * np.exp(- ( (x-self.symm)/self.width )**2 )
    
    def create_gauss_x_package(self):
        return np.array([self.create_gauss_at_x(pos) for pos in self.x_grid])
    
    def plot_x_package(self):
        
        #psi_abs_squared_x = np.multiply(self.x_package, self.x_package.conj()).real
        
        x_package = plt.plot(self.x_grid, np.abs(self.x_package))
        plt.title ("Propabilitydensity of gaussian package at time zero in position-space")
        plt.xlabel("position grid")
        plt.ylabel("probabilitydensity")
        
        plt.show(x_package)
