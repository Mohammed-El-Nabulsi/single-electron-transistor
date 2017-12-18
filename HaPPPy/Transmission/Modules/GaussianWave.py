import cmath
from scipy.constants import codata

me   = codata.value("electron mass energy equivalent in MeV") * 1e9 ;  # Convert to milli eV
hbar = codata.value("Planck constant over 2 pi in eV s")      * 1e15;  # Convert to ps

class GaussianWave():
    def create_package_at_point(self, x, x0, energy, a):
        k = cmath.sqrt(2*me*energy)/hbar
        norm = (2/(cmath.pi*(a**2)))**(1/4)
        
        psi_x = complex(norm*cmath.exp(1j*k*(x-x0))*cmath.exp(-((x-x0)/a)**2))
        
        
        return psi_x
    
    def create_package(self, positions, x0, energy, a):
        '''
        Creates a gaussian package in position-space symmetrically around x0 at time zero.

        Parameters
        ----------
        positions : array
           A sequence of positions for which to create the package
        x0 : float
           symmetry point of the package           
        energy : float
           Energy of the package
        a : float
           Width of the package
        
        Returns
        -------
        psi : array, shape(len(positions))
           Array containing the value of the gaussian package for each desired position in positions,
        '''
        psi = [self.create_package_at_point(x, x0, energy, a) for x in positions]
 
        return psi
    

