import numpy as np
from HaPPPy.Transmission.Modules.Fourier import Fourier

from scipy.constants import codata

dt = 10**(2)
me   = codata.value("electron mass energy equivalent in MeV") * 1e8 ;  # Convert to meV
hbar = codata.value("Planck constant over 2 pi in eV s")      * 1e19;  # Convert to 10 meV fm


class Split_Step_Operator():
    """
    Parameters
    ----------
    position_grid: Array
        samplepoints in position-space
    pot : Array, len(position_grid)
        potentialalues at positiongrid
        
    Attributes
    ----------
    k : Array, len(position_grid)
        samplepoints in wavenumber-space
    kinetic_operator : Array, len(position_grid)
        diagonalelements of the kinetic energie operator at wavenumber-grid
    potential_operator: Array, len(position_grid)
        diagonalelements of the potential energie operator at position-grid
    """
    
    def __init__(self, position_grid, pot):
        
        self.const = hbar/(2*me)
        
        self.fourier = Fourier(position_grid)
        self.k = self.fourier.k
        
        self.pot = pot
        
        def kinetic_diagonal_element(self, k): #to be used for steps execpt the first and last
            return np.exp(-1j*self.const*dt*k**2)
        
        def kinetic_diagonal_element_half(self, k): #to be used in the first and last step
            return np.exp(-1j*(self.const/2)*dt*k**2)

        def potential_diagonal_element(self, v):
            return np.exp(-1j*dt*v/hbar)
        
        def create_kinetic_operator(self):
            kinetic_operator = np.zeros(self.k.size, dtype=np.complex64)
            
            indices = np.arange(self.k.size)
            
            for a in indices:
                kinetic_operator[a]=kinetic_diagonal_element(self,
                                                             self.fourier.k[a])
            return kinetic_operator
        
        def create_kinetic_operator_half(self):
            kinetic_operator_half = np.zeros(self.k.size, dtype=np.complex64)
            
            indices = np.arange(self.k.size)
            
            for a in indices:
                kinetic_operator_half[a]=kinetic_diagonal_element_half(self,
                                                             self.fourier.k[a])
            return kinetic_operator_half
        
        self.kinetic_operator = create_kinetic_operator(self)
        self.kinetic_operator_half = create_kinetic_operator_half(self)
        
#        print("shapes of the kinetic operators")
#        print(self.kinetic_operator.shape)
#        print(self.kinetic_operator_half.shape)
#        print("types of the kinetic operators")
#        print(type(self.kinetic_operator))
#        print(type(self.kinetic_operator_half))

        def create_potential_operator(self):
            potential_operator = np.zeros(self.pot.size, dtype=np.complex64)
            
            indices = np.arange(self.pot.size)
            
            for b in indices:
                potential_operator[b]=potential_diagonal_element(self,
                                                             self.pot[b])
            return potential_operator
        
        self.potential_operator = create_potential_operator(self)
        
#        print("shape of the potenital operator")
#        print(self.potential_operator.shape)
#        print("type of the potential operator")
#        print(type(self.potential_operator))
    
    def first_step(self, x_wave_1):
        """
        First interation step.
        
        Parameters
        ----------
        xx_wave_1: Array
            Gaussian package
            
        Returns
        --------
        ft_pot_x_wave : Array, len(xx_wave_1)
            Wavefunction after applying the half kinetc- and full potentialoperator
        """
        k_wave_f = self.fourier.dft(x_wave_1)
        kin_k_wave = np.multiply(self.kinetic_operator_half, k_wave_f, dtype=np.complex64)
        return self.fourier.idft(kin_k_wave)
    
    def steps(self, wave):
        """
        Performs the Split Step iteration. This command combines sequentices fourier transformations.
        
        Parameter
        ---------
        wave : Array
            Wavefunction in wavenumberspace
        Returns
        -------
        v_k_wave_s2 : Array, len(wave)
            wavefunction after applying full kinetic-, full potential- and again full kineticoperator
        """
        pot_wave_s = np.multiply(self.potential_operator, wave, dtype=np.complex64)
        k_wave_s = self.fourier.dft(pot_wave_s)
        kin_wave_s = np.multiply(self.kinetic_operator, k_wave_s, dtype=np.complex64)
        return self.fourier.idft(kin_wave_s)

    
    def final_step(self, k_wave_fs):
        """
        Concludes the last full timestep of the splitstep-iteration.
        
        Parameter
        ---------
        k_wave_fs : Array
            Wavefunction after reacing the stopping criterion.
        Returns
        -------
        ift_kin_k_wave_fs : Array, len(k_wave_fs)
            Wavefuntion in positionspace after the split-step-algorithm
        """
        pot_wave_fs = np.multiply(self.potential_operator, k_wave_fs, dtype=np.complex64)
        k_wave_fs = self.fourier.dft(pot_wave_fs)
        kin_wave_fs = np.multiply(self.kinetic_operator_half, k_wave_fs, dtype=np.complex64)
        return self.fourier.idft(kin_wave_fs)
