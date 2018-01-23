import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq, fftshift, ifft

class Fourier:
    def __init__(self,x):
        """
        Performs the Transformation between position- and wavenumberspace.
        
        Parameters
        ----------
        x : Array
            Grid in positionspace
        """

        self.N = x.size
        self.dx = x[1]-x[0]
        self.x = x
        self.k_stand = 2 * np.pi * np.fft.fftfreq(self.N, self.dx)
        self.k = fftshift(self.k_stand)
            
    def dft(self,f):
        """
        Parameter
        ---------
        f : Array
        Wavefunction in positionspace
        
        Returns
        -------
        F : Array (dtype=complex_), len(f)
        Wavefunction in wavenumberspace
        """
        yf = np.fft.fft(f,norm='ortho')
        F = fftshift(yf)
        return F

    def plot_dft(self,f):
        Fplot = self.dft(f)
        plt.plot(self.k,np.abs(Fplot))
        plt.grid()
        plt.show()
       
    def idft(self,F):
        """
        Parameter
        ---------
        f : Array
        Wavefunction in wavenumberspace
        
        Returns
        -------
        F : Array (dtype=complex_), len(f)
        Wavefunction in positionspace
        """
        yf = np.fft.ifft(F,norm='ortho')
        return yf
    
    def plot_idft(self,F):

        fplot = self.idft(F)
        plt.plot(self.x,np.abs(fplot))
        plt.grid()
        plt.show()
