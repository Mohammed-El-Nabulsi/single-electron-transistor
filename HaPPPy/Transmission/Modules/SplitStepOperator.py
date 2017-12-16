#Split-Step-Algorithm repeated M times
import numpy
import cmath

#diagonalelemts of kinetic operator in k-space
def t_element(ki):
    
    t_element = cmath.exp(-1j*hbar*dt*(ki**2)/(4*me))
    
    return t_element

#diagonalelements of potential in position-space
def v_element(Vi):
    
    v_element = cmath.exp(-1j*dt*vi/hbar)
    
    return v_element

#SplitStep-Operator
def splitStep(wave):
    
    #transform into k space
    k_wave = FT(wave)
    
    #applying modified kinetic operator
    tk_wave_new = numpy.multiply(k_wave,T, dtype=numpy.complex64)
    
    #transform into x space
    x_wave = IFT(tk_wave)
    
    #applying potential operator
    vx_wave = numpy.multiply(x_wave,V, dtype=numpy.complex64)
    
    #transform into k space
    k_wave_new = FT(vx_wave)
    
    #applying modified kinetic operator
    tk_wave_new_2 = numpy.multiply(k_wave_new,T, dtype=numpy.complex64)
    
    #transform into x space
    x_wave_new = IFT(tk_nave_new_2)
    
    return x_wave_new

dt = 1 #timestep
me = 1 #mass of electron
hbar = 1 #reduced plack constant

positions = [1]
wave_numbers = [1]

potential = [1] #potential

#diagonalelements of kinetic operator
T = [t_element(wave_number) for wave_number in wave_numbers]

#diagonalelements of potential operator
V = [v_element(vi) for vi in potential]
