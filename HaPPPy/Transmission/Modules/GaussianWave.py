import cmath # Calculations with complex numbers
from scipy.constants import codata

me   = codata.value("electron mass energy equivalent in MeV") * 1e9 ;  # Convert to milli eV
hbar = codata.value("Planck constant over 2 pi in eV s")      * 1e15;  # Convert to ps

class GaussianWave():
    def create_package(self, x, energy, a):
        k = cmath.sqrt(2*me*energy)/hbar
        norm = (2/(cmath.pi*(a**2)))**(1/4)
        
        psi = complex(norm*cmath.exp(1j*k*x-(x/a)**2))
        
        return psi
