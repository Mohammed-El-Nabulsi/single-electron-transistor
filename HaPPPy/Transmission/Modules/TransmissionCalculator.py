import numpy
#Integration by addition of sample points.
def integrate(funktion,start): # funktion = array with the samplepoints from the funktion , start sets the startingpoint for the integration.
    n = start
    a = 0
    # Squares the absolut value and adds the samplepoints.
    while n < numpy.size(funktion) : 
        a = a + funktion[n]*funktion.conjugate()[n]
        n = n+1
    return a
    
    # Calculates the propability for tunneling with a given final wavefunktion.
def warhscheinlichkeit(funktion,start):#  Requires only the final wavefunktion in the form of an array and the startingpoint for the propability.
    x = integrate(funktion,start)/ integrate(funktion,0)
    return x
