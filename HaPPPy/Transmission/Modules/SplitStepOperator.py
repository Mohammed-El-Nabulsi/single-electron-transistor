<<<<<<< HEAD
import numpy
import cmath
from scipy.constants import codata
<<<<<<< HEAD
from numpy.fft import fft, ifft

<<<<<<< HEAD
dt = 10**-15
me   = codata.value("electron mass energy equivalent in MeV") * 1e9 ;  # Convert to milli eV
hbar = codata.value("Planck constant over 2 pi in eV s")      * 1e15;  # Convert to p    
    
class SplitStepOperator:
    def __init__(self, V, x):
        def create_v_element(v_i):
            return cmath.exp(-1j*dt*v_i/hbar)
=======
=======
>>>>>>> 6640497... Update SplitStepOperator.py




dt = 10**-15
me   = codata.value("electron mass energy equivalent in MeV") * 1e9 ;  # Convert to milli eV
hbar = codata.value("Planck constant over 2 pi in eV s")      * 1e15;  # Convert to p    
x =[1, 2, 3, 4, 5, 6] # example, getting it from Mo
k = FT(x)
v = [1, 2, 3, 4, 5, 6] #example, gettoing from mo
    
    
    
    
class SplitStepMethod:

    
    
    def use(self, psi):
        
        

        def v_element_create(self, vi):
       
        
            v_element = cmath.exp(-1j*dt*vi/hbar)
        
<<<<<<< HEAD
        return v_element
    
    
    def use(self, psi, T, V):
>>>>>>> 4625e5b... Update SplitStepOperator.py
        
        def create_t_element(k_i):
            return cmath.exp(-1j*hbar*dt*(k_i**2)/(4*me))

        k = fft(x)

        self.v_elements = [create_v_element(v_i) for v_i in V]
        self.t_elements = [create_t_element(k_i) for k_i in k]


    def use(self, psi):
        psi_k = fft(psi, norm='ortho')
        
        psi_tk = numpy.multiply(psi_k, self.t_elements,  dtype=numpy.complex64)
=======
            return v_element
        
        v_elements = [v_element_create(vi) for vi in v]
        
        
        
        
        
        def t_element_create(self, ki):
            
            t_element = cmath.exp(-1j*hbar*dt(ki**2)/(4*me))
            
            return t_element
        
        t_elements = [t_element_create(ki) for ki in k]
            
            
            
        psi_k = FT(psi)
        
        psi_tk = numpy.multiply(psi_k, t_elements,  dtype=numpy.complex64)
>>>>>>> 6640497... Update SplitStepOperator.py
        
        psi_x = ifft(psi_tk, norm='ortho')
        
<<<<<<< HEAD
        psi_vx = numpy.multiply(psi_x, self.v_elements, dtype=numpy.complex64)
=======
        psi_vx = numpy.multiply(psi_x, v_elements, dtype=numpy.complex64)
>>>>>>> 6640497... Update SplitStepOperator.py
        
        psi_k_new = fft(psi_vx, norm='ortho')
        
<<<<<<< HEAD
        psi_tk_new = numpy.multiply(psi_k_new, self.t_elements, dtype=numpy.complex64)
=======
        psi_tk_new = numpy.multiply(psi_k_new, t_elements, dtype=numpy.complex64)
>>>>>>> 6640497... Update SplitStepOperator.py
        
        psi_x_new = ifft(psi_tk_new, norm='ortho')
        
        return psi_x_new
=======
#Split-Step-Algorithm repeated M times
import numpy
import cmath

#diagonalelemts of kinetic operator in k-space
def TElement(ki):
    
    Telement = cmath.exp(-1j*hbar*dt*(ki**2)/(4*me))
    
    return Telement

#diagonalelements of potential in position-space
def VElement(Vi):
    
    Velement = cmath.exp(-1j*dt*Vi/hbar)
    
    return Velement

#SplitStep-Operator
def splitStep(wave):
    
    #transform into k space
    kWave1 = FT(wave)
    
    #applying modified kinetic operator
    TkWave1 = numpy.multiply(kWave1,T, dtype=numpy.complex64)
    
    #transform into x space
    xWave = IFT(TkWave1)
    
    #applying potential operator
    VxWave = numpy.multiply(xWave1,V, dtype=numpy.complex64)
    
    return VxWave

dt = 1 #timestep
me = 1 #mass of electron
hbar = 1 #reduced plack constant

positions = [1]
waveNumbers = [1]

potential = [1] #potential

#diagonalelements of kinetic operator
T = [TElement(waveNumber) for waveNumber in waveNumbers]

#diagonalelements of potential operator
V = [VElement(vi) for vi in potential]
>>>>>>> d7fea88... Update SplitStepOperator.py
