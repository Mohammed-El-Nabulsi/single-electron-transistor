import numpy

class Transmission():
    def integrate(self, wavefunction,start): 
        """ Documentation"""
        n = start
        a = 0
        while n < numpy.size(wavefunction) : 
            a = a + wavefunction[n]*numpy.array(wavefunction).conjugate()[n]
            n = n+1
        return a

        
    def calculate(self, wavefunction,start):
        """Documentation"""
        x = self.integrate(wavefunction,start) / self.integrate(wavefunction,0)
        return x.real
