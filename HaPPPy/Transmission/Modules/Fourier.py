import numpy as np
from scipy.fftpack import fft, fftfreq, fftshift

import matplotlib.pyplot as plt

class Fourier:
    
    def __init__(self, x_grid, dx):
        self.N = x_grid.size
        self.x_grid = x_grid
        self.dx = dx
        
        self.k = self.get_wavenumbers()
        
    def get_wavenumbers(self):
        xf = fftfreq(self.N, self.dx)
        return fftshift(xf)
    
    def dft(self, x_func):
        k_func = fft(x_func)
        return fftshift(k_func)
    
    def plot_dft_abs_squared(self, x_func):
        dft_func = self.dft(x_func)
        dft_abs_squared = np.multiply(dft_func, dft_func.conj()).real
        
        plot_dft = plt.plot(self.k, dft_abs_squared)
        plt.title ("Propabilitydensity of the transformed wave in wavenumber-space")
        plt.xlabel("wavenumber grid")
        plt.ylabel("probabilitydensity")
        
        plt.show()
