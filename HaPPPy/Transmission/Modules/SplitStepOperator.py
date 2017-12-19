import numpy
import cmath
from scipy.constants import codata

from HaPPPy.Transmission.Modules.GaussianWave import Fourier

Fourier = Fourier()

dt = 10**-15
me   = codata.value("electron mass energy equivalent in MeV") * 1e9 ;  # Convert to milli eV
hbar = codata.value("Planck constant over 2 pi in eV s")      * 1e15;  # Convert to p    
    
class SplitStepOperator:
    def __init__(self, V, psi, x, x0, L):
        def create_v_element(v_i):
            return cmath.exp(-1j*dt*v_i/hbar)
        
        def create_t_element(k_i):
            return cmath.exp(-1j*hbar*dt*(k_i**2)/(4*me))

        k = Fourier.waveNumbers(psi, x0, L)

        self.v_elements = [create_v_element(v_i) for v_i in V]
        self.t_elements = [create_t_element(k_i) for k_i in k]

    def use(self, psi, x0, L):
        psi_k = FT(psi, x0, L)
        
        psi_tk = numpy.multiply(psi_k, self.t_elements,  dtype=numpy.complex64)
        
        psi_x = IFT(psi_tk, x0, L)
        
        psi_vx = numpy.multiply(psi_x, self.v_elements, dtype=numpy.complex64)
        
        psi_k_new = FT(psi_vx, x0, L)
        
        psi_tk_new = numpy.multiply(psi_k_new, self.t_elements, dtype=numpy.complex64)
        
        psi_x_new = IFT(psi_tk_new, x0, L)
        
        return psi_x_new
