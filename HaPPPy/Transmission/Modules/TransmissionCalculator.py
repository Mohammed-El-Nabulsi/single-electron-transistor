import numpy

class TransmissionCalculator():
    def integrate(wavefunktion,start): 
        """ Documentation"""
        n = start
        a = 0
       while n < numpy.size(funktion) : 
            a = a + wavefunktion[n]*wavefunktion.conjugate()[n]
            n = n+1
        return a

        
    def calculate(wavefunktion,start):
        """Documentation"""
        x = TransmissionCalculator.integrate(wavefunktion,start) / TransmissionCalculator.integrate(wavefunktion,0)
        return x.real
