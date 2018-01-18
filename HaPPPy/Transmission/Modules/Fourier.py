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
        yf = np.fft.fft(f,norm='ortho')#, norm='ortho')# * 1/self.N #self.dx/np.sqrt(2*np.pi)
        F = fftshift(yf) #* self.dx/np.sqrt(2*np.pi) #1/self.N
#        F[0] = 0
#        F[-1] = 0
        #print(self.dx*np.sum(np.abs(F)))
        return F

    def plot_dft(self,f):
#        yf = np.fft.fft(f)#, norm='ortho')
#        Fplot = fftshift(yf)
        Fplot = self.dft(f)
        plt.plot(self.k,np.abs(Fplot))
        plt.grid()
        plt.show()
       # print(max(abs(Fplot)))
       
    def idft(self,F):
        yf = np.fft.ifft(F,norm='ortho')
        #f = fftshift(yf)#*self.N
#        yf[0]=0
#        yf[-1]=0
        return yf
    
    def plot_idft(self,F):
#        yf = np.fft.fft(f)#, norm='ortho')
#        Fplot = fftshift(yf)
        fplot = self.idft(F)
        plt.plot(self.x,np.abs(fplot))
        plt.grid()
        plt.show()
       # print(max(abs(Fplot)))
       
#x = np.arange(0,15,0.1)   
#f = 0.5*np.sin(1*2*np.pi  *x) #+ np.exp(1j*20*x) #np.exp(2*1j*2*np.pi*x-x**2)#
#
#fourier = Fourier(x)
#a=fourier.dft(f)
#fourier.plot_idft(a)

#a = fourier.dft(f)

#fourier.plot_idft(a)
