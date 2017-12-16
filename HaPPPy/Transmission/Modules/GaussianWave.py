import cmath #calculations with complex numbers
from scipy.constants import codata

me = codata.value("electron mass energy equivalent in MeV");
hbar = codata.value("Planck constant over 2 pi in eV s") * 1e15; # Convert to ps

# Is this wrong? Electron mass in MeV is 1e-1 not 1e-9.
me = 0.5109989461*(10**(9))            #mass of electron in meV

class GaussianWave():
    def create_package(self, x, energy, a):
        k = cmath.sqrt(2*me*energy)/hbar
        norm = (2/(cmath.pi*(a**2)))**(1/4)
        
        psi = complex(norm*cmath.exp(1j*k*x-(x/a)**2))
        
        return psi
