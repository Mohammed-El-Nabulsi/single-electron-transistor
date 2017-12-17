import numpy
import cmath
from scipy.constants import codata
from numpy.fft import fft, ifft

dt = 10**-15
me   = codata.value("electron mass energy equivalent in MeV") * 1e9 ;  # Convert to milli eV
hbar = codata.value("Planck constant over 2 pi in eV s")      * 1e15;  # Convert to p    
    
class SplitStepOperator:
    def __init__(self, V, x):
        def create_v_element(v_i):
            return cmath.exp(-1j*dt*v_i/hbar)
        
        def create_t_element(k_i):
            return cmath.exp(-1j*hbar*dt*(k_i**2)/(4*me))

        k = fft(x)

        self.v_elements = [create_v_element(v_i) for v_i in V]
        self.t_elements = [create_t_element(k_i) for k_i in k]


    def use(self, psi):
        psi_k = fft(psi, norm='ortho')
        
        psi_tk = numpy.multiply(psi_k, self.t_elements,  dtype=numpy.complex64)
        
        psi_x = ifft(psi_tk, norm='ortho')
        
        psi_vx = numpy.multiply(psi_x, self.v_elements, dtype=numpy.complex64)
        
        psi_k_new = fft(psi_vx, norm='ortho')
        
        psi_tk_new = numpy.multiply(psi_k_new, self.t_elements, dtype=numpy.complex64)
        
        psi_x_new = ifft(psi_tk_new, norm='ortho')
        
        return psi_x_new
    
