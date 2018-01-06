import numpy as np
from scipy.constants import codata

import matplotlib.pyplot as plt

me   = codata.value("electron mass energy equivalent in MeV") * 1e9 ;  # Convert to milli eV
hbar = codata.value("Planck constant over 2 pi in eV s")      * 1e15;  # Convert to ps


from HaPPPy.Transmission.Modules.PotentialPreparation import Potential

class GaussianWave():
    
    def __init__(self, width, symm, energy, x_grid, k_grid):

        self.width = width
        self.symm = symm
        self.x_grid = x_grid
        
        self.k_grid = k_grid
        self.k0 = np.sqrt(2*me*energy)/hbar
        
        self.x_package = self.create_gauss_x_package()
        self.k_package = self.create_gauss_k_package()
    
    def create_gauss_at_x(self,x):
        norm_x = ( 2/(np.pi*self.width**2) )**(1/4)
        return norm_x * np.exp(1j*self.k0*(x-self.symm)) * np.exp(- ( (x-self.symm)/self.width )**2 ) #verschiebung um symm in zweiter exp vergessen
    
    def create_gauss_x_package(self):
        return np.array([self.create_gauss_at_x(pos) for pos in self.x_grid])
    
    def plot_x_package(self):
        
        psi_abs_squared_x = np.multiply(self.x_package, self.x_package.conj()).real
        
        x_package = plt.plot(self.x_grid, psi_abs_squared_x)
        plt.title ("Propabilitydensity of gaussian package at time zero in position-space")
        plt.xlabel("position grid")
        plt.ylabel("probabilitydensity")
        
        plt.show()
        
    def create_gauss_at_k(self, k):
        norm_k = np.sqrt(self.width)/(2*np.pi)**(1/4)
        return norm_k*np.exp( -(self.width*(k-self.k0)/2 )**2 )*np.exp( -1j*(k-self.k0)*self.symm)
    
    def create_gauss_k_package(self):
        return np.array([self.create_gauss_at_k(wnum) for wnum in self.k_grid])
    
    def plot_k_package(self):
        
        psi_abs_squared_k = np.multiply(self.k_package, self.k_package.conj()).real
        
#        print(psi_abs_squared_k[98])
#        print(psi_abs_squared_k[99])
#        print(psi_abs_squared_k[100])
#        print(max(psi_abs_squared_k)-psi_abs_squared_k[99])
        
        k_package = plt.plot(self.k_grid, psi_abs_squared_k)
        plt.title("Propabilitydensity of gaussian package at time zero in wavenumber-space")
        plt.xlabel("wavenumber grid")
        plt.ylabel("probabilitydensity")
        
        plt.show()
