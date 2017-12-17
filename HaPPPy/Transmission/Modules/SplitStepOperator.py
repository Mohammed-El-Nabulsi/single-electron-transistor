import numpy
import cmath


class SplitStepMethod:
    dt = 1
    me = 1
    hbar = 1
    
    positions = [1]
    wave_numbers = [1]
    potential = [1]
    
    T = [t_element(wave_number) for wave_number in wave_numbers]
    
    V = [v_element(vi) for vi in potential]
    
    
    def t_element_create(self, dt, me, hbar):
    
        t_element = cmath.exp(-1j*hbar*dt(ki**2)/(4*me))
        
        return t_element
    
    
    
    def v_element_create(self, dt, me, hbar):
       
        
        v_element = cmath.exp(-1j*dt*vi/hbar)
        
        return v_element
    
    
    def use(self, psi, T, V):
        
        psi_k = FT(psi)
        
        psi_tk = numpy.multiply(psi_k, T,  dtype=numpy.complex64)
        
        psi_x = IFT(psi_tk)
        
        psi_vx = numpy.multiply(psi_x, V, dtype=numpy.complex64)
        
        psi_k_new = FT(psi_vx)
        
        psi_tk_new = numpy.multiply(psy_k_new, T, dtype=numpy.complex64)
        
        psi_x_new = IFT(psi_tk_new)
        
        return psi_x_new
    
