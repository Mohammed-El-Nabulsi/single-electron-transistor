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

    #transform into k space
    kWave1 = FT(wave)
    
def split_step(wave):
    
    
    #applying modified kinetic operator
    TkWave1 = numpy.multiply(kWave1,T, dtype=numpy.complex64)
    
    #transform into x space
    xWave = IFT(TkWave1)
    
    #applying potential operator
    VxWave = numpy.multiply(xWave1,V, dtype=numpy.complex64)
    
    #applying modified kinetic operator
    TkWave1 = numpy.multiply(kWave1,T, dtype=numpy.complex64)
    
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
