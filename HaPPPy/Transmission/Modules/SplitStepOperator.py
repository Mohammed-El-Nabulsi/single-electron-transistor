import numpy
import cmath
from scipy.constants import codata

from HaPPPy.Transmission.Modules.Fourier import Fourier

dt = 10**-15
me   = codata.value("electron mass energy equivalent in MeV") * 1e9 ;  # Convert to milli eV
hbar = codata.value("Planck constant over 2 pi in eV s")      * 1e15;  # Convert to p    
    
class SplitStepOperator:
    def __init__(self, V, psi, x, x0, L):
       
        def create_v_element(v_i):
            return cmath.exp(-1j*dt*v_i/hbar)
        
        def create_t_element(k_i):
            return cmath.exp(-1j*hbar*dt*(k_i**2)/(4*me))
       
        self.fourier = Fourier()

        k = self.fourier.waveNumbers(psi, x0, L)

        self.v_elements = [create_v_element(v_i) for v_i in V]
        self.t_elements = [create_t_element(k_i) for k_i in k]

        self.x0 = x0
        self.L = L

    def use(self, psi):
        '''
        Applying Split Step Method
        
        Parameters
        ----------
        psi : Array
            A sequence of wavefunction values in position space 
        x0 : float
            lowest sample point argument in position space         
        L : float
            range of the sample points in position space
            
        Returns
        -------
        psi_x_new : array, shape(len(psi))
            Array containing the wavefunction values after one timestep
        '''   
        x0 = self.x0
        L = self.L

        psi_k = self.fourier.FT(psi, x0, L)
        
        psi_tk = numpy.multiply(psi_k, self.t_elements,  dtype=numpy.complex64)
        
        psi_x = self.fourier.IFT(psi_tk, x0, L)
        
        psi_vx = numpy.multiply(psi_x, self.v_elements, dtype=numpy.complex64)
        
        psi_k_new = self.fourier.FT(psi_vx, x0, L)
        
        psi_tk_new = numpy.multiply(psi_k_new, self.t_elements, dtype=numpy.complex64)
        
        psi_x_new = self.fourier.IFT(psi_tk_new, x0, L)
        
        return psi_x_new
