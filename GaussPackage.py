
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 17:20:48 2017
@author: Frederik Weißler
"""
import cmath #calculations with complex numbers
#import numpy
#create gaussian package
def gauss(x):
    k = cmath.sqrt(2*me*energie)/hbar
    norm = (2/(cmath.pi*(a**2)))**(1/4)
    
    psi = complex(norm*cmath.exp(1j*k*x-(x/a)**2)) #wavefunktion at position x
    
    return psi
#compare returned gauss(x) with testPsi
def test(z):
    absReal = abs(z.real - testReal)
    absImag = abs(z.imag - testImag)
    uncertain = 10**(-10)
    if absReal < uncertain  and absImag < uncertain:
        return True
    else:
        return False
#----------------------------------------------------------------------------
#testPsi for energie = x = a = 1 from wolfram alpha:
#testReal = 0.328573471238826838814724884575412086439812622961237236
#testImag = -0.0046267725731805804274342494105112176810998628709738939
#----------------------------------------------------------------------------
 
#const source: http://pdg.lbl.gov/2017/reviews/rpp2016-rev-phys-constants.pdf     
me = 0.5109989461*(10**(9))            #mass of electron in meV
hbar = 6.582119514*(10**(-1))          #reduced plack constant in meV*ps
#userinput: energie an width of the package to calc psi at position x
energie = float(input('Energie of the package in meV?'))
a = float(input('Width of the package?')) #width of the package
x = float(input('position in pm?'))
#userinput: expected values for wavepackage at position x
testReal = float(input('Test: realpart?')) #Realpart of known wavepackage
testImag = float(input('Test: imaginarypart?')) #Imagpart of known wavepackage
#calc and print complex value of gaussian package at position x
psi = gauss(x)
print('psi(x=',x,',t=0)=',psi)
#compare calculated psi to expected values and print the results
testPsi = complex(testReal + testImag*1j)
print('######################################################################')
print('Testresult:',test(psi),'.')
print('Allowed uncertain: 10**(-10)')
print('Difference psi-testPsi =', psi-testPsi)
print('######################################################################')
