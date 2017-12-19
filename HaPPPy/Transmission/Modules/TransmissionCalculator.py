import numpy

class TransmissionCalculator():
    def integrate(wavefunction,start): 
        """ Documentation"""
        n = start
        a = 0
        while n < numpy.size(wavefunction) : 
            a = a + wavefunction[n]*wavefunction.conjugate()[n]
            n = n+1
        return a

        
    def calculate(wavefunction,start):
        """Documentation"""
        x = TransmissionCalculator.integrate(wavefunction,start) / TransmissionCalculator.integrate(wavefunction,0)
        return x.real
