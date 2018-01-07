import numpy as np
from scipy.fftpack import fft, fftfreq, fftshift
import matplotlib.pyplot as plt

class Fourier:
    
    def __init__(self, N, dx):
        """
        Parameters
        ----------
        N : Integer
            len(x_grid), number of position sample points
        dx : float
            distance between two neighbouring samplepoints in position space grid
        """
        
        self.N = N
        self.dx = dx
        self.k = self.get_wavenumbers()
        
    def get_wavenumbers(self):
        xf = fftfreq(self.N, self.dx)
        return fftshift(xf)

    def dft(self, x_func):
        k_func = fft(x_func)
        return fftshift(k_func)*1/self.N
        
    def plot_dft(self, x_func):
        dft = self.dft(x_func)
        plt.plot(self.k,np.multiply(dft,dft.conj()).real)
        plt.grid()
        plt.show()


#for Testing: (works!)
#fourier = Fourier(500,1/500)
#
#x = np.linspace(0, 1, 500)
#y = np.exp(50.0 * 1.j * 2.0*np.pi*x) + 0.6*np.exp(-25.0 * 1.j * 2.0*np.pi*x)
#
#fourier.plot_dft(y)
#
### number of signal points
##N = 200
##L = 1
### sample spacing
##dx = L / (N)
##
##x = np.linspace(0.0, N*dx, N)
##
##print(N*dx)
##
##
##y = np.exp(50.0 * 1.j * 2.0*np.pi*x) + 0.6*np.exp(-25.0 * 1.j * 2.0*np.pi*x)
##
##yf = fft(y)
##xf = fftfreq(N, dx)
##xf = fftshift(xf)
##yplot = fftshift(yf)
##import matplotlib.pyplot as plt
###plt.plot(x,np.abs(y))
##plt.show()
##plt.plot(xf, 1.0/(N) * np.abs(yplot))
##plt.grid()
##plt.show()
