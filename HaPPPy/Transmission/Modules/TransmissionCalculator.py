import numpy
#Integration über die Funktion anhand von Stützstellen.
def integrate(funktion,start): # funktion  = Funktionsarray, start ist der endpunkt des Potentials.
    n = start
    a = 0
    # Komplexkonjugiert die werte und addiert sie auf.
    while n < numpy.size(funktion) : 
        a = a + funktion[n]*funktion.conjugate()[n]
        n = n+1
    return a
    
    # berechnet die Wahrscheinlicheit anhand der addierten Stützstellen.
def warhscheinlichkeit(funktion,start):#  benötigt die momentane Funktion und den Endpunkt des Potentials
    x = integrate(funktion,start)/ integrate(funktion,0)
    return x
 
#def main():
#    funktion = numpy.array([1,1,1,1,1,1,1,1,1])
#    start = 4
#    print (warhscheinlichkeit(funktion,start) )
#    
#main()
#    
