import cmath # Calculations with complex numbers
from scipy.constants import codata

me   = codata.value("electron mass energy equivalent in MeV") * 1e9 ;  # Convert to milli eV
hbar = codata.value("Planck constant over 2 pi in eV s")      * 1e15;  # Convert to ps

class GaussianWave():
    def CreatePackagePoint(self, x, x0, energy, a):
        k = cmath.sqrt(2*me*energy)/hbar
        norm = (2/(cmath.pi*(a**2)))**(1/4)
        
        psi = complex(norm*cmath.exp(1j*k*(x-x0)-((x-x0)/a)**2))
        
        return psi
    
    def createPackage(positions, x0, energy, a):
        phi = [createPackagePoint(x, x0, energy, a) for x in positions]

        return phi
    
 '''
 Creates a gaussian package in position-space symmetrically around x0 at time zero.
 
 Parameters
 ----------
 x : float
    Position on which the gaussian wave is evaluated
 x0 : float
    symmetry point of the package
 positions : array
    A sequence of positions for which to create the package
 energy: float
    Energy of the package
 a : float
    Width of the package
 
 Returns
 -------
 wavePoint : complex
    Value of the gaussian wave at time zero an position x
 wave : array, shape(len(positions))
    Array containing the value of the gaussian package for each desired position in positions,
'''
