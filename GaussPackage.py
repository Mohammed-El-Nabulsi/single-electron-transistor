# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 16:14:35 2017

@author: Frederik Weißler
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 18:56:40 2017

@author: Frederik Weißler
"""
import cmath #calculations with complex numbers
import numpy
import matplotlib.pyplot as plt


#create 1D-gaussian-package at starting time
def gauss(x):
    me = 0.5109989461*(10**(9))            #mass of electron in meV
    hbar = 6.582119514*(10**(-1))          #reduced plack constant in meV*ps
    
    k = cmath.sqrt(2*me*energie)/hbar
    norm = (2/(cmath.pi*(a**2)))**(1/4)
    
    psi = complex(norm*cmath.exp(1j*k*x-(x/a)**2)) #wavefunktion at position x
    
    return psi
########################################
#encomment for visualization
#calc the squared abs value of an complex number
def absSqare(x):
    z = x*x.conjugate()
    return z
#testing the function 'gauss' for an well known package
#########################################
    
#def test(x):
#    diffReal = abs(gauss(x).real - 0.32857347123882683881472488457541208643981)
#    diffImag = abs(gauss(x).imag -(-0.0046267725731805804274342494105112176810)
#    if diffReal < 10**(-10): #and diffImag < 10**(-10) :
#        return True
#    else:
#        return False
    
#main
energie = float(input('Energie of the package in meV?'))
a = float(input('Width of the package?')) #width of the package

##########################################
#Vorschlag: Breite des Pakets = Breite der Barriere/N mit natuerlichem N
# N durch Vergleich mit bekannter Transmission von z.B. PotTopf finden.
##########################################

#########################################
#visualizing the probability-density
x = numpy.arange(0,6,1) 
psiSquared = [0,1,2,3,4,5]

for i in x:
    psiSquared[i] = absSqare(gauss(i)).real

line1 = plt.plot(x, psiSquared)

plt.show()
print(gauss(1))
#########################################
#testing
#plot(test(x))
