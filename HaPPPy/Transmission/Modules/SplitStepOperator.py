import numpy
import cmath
from scipy.constants import codata

class SplitStepMethod:
    dt = 10**-15
    me   = codata.value("electron mass energy equivalent in MeV") * 1e9 ;  # Convert to milli eV
    hbar = codata.value("Planck constant over 2 pi in eV s")      * 1e15;  # Convert to ps
    
    x =[1, 2, 3, 4, 5, 6] # example, getting it from Mo
    k = FT(x)
    v = [1, 2, 3, 4, 5, 6] #example, gettoing from mo
    
    
    
    def use(self, psi):
        
        

        def v_element_create(self, dt, me, hbar, vi, v):
       
        
            v_element = cmath.exp(-1j*dt*vi/hbar)
        
            return v_element
        
        v_elements = [v_element_create(vi) for vi in v]
        
        
        
        
        
        def t_element_create(self, dt, me, hbar, ki, k):
            
            t_element = cmath.exp(-1j*hbar*dt(ki**2)/(4*me))
            
            return t_element
        
        t_elements = [t_element_create(ki) for ki in k]
            
            
            
        psi_k = FT(psi)
        
        psi_tk = numpy.multiply(psi_k, t_elements,  dtype=numpy.complex64)
        
        psi_x = IFT(psi_tk)
        
        psi_vx = numpy.multiply(psi_x, v_elements, dtype=numpy.complex64)
        
        psi_k_new = FT(psi_vx)
        
        psi_tk_new = numpy.multiply(psi_k_new, t_elements, dtype=numpy.complex64)
        
        psi_x_new = IFT(psi_tk_new)
        
        return psi_x_new
    
