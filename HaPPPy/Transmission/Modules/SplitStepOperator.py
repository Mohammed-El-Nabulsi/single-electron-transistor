#TODO: docu, test(kommt nach Silvester, habe nun endlich Maple...)
#In dieser Version ist berücksichtigt, dass im Algo idft und dft aufeinander folgen.
#auf diese weise wird eine vielzahl an Operationen eingespart.
#Ausfuehrliche Erklärungen folgen, docu hier unbrauchbar/zT falsch(Relikt aus Experimenten)

import numpy as np

from scipy.constants import codata

dt = 10**(-15)
me   = codata.value("electron mass energy equivalent in MeV") * 1e9 ;  # Convert to milli eV
hbar = codata.value("Planck constant over 2 pi in eV s")      * 1e15;  # Convert to p    

class SplitStepOperator():
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
        
        self.k = np.fft.fft(position_grid, norm='ortho')
        
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
                                                             self.k[a])
            return kinetic_operator
        
        def create_kinetic_operator_half(self):
            kinetic_operator_half = np.zeros(self.k.size, dtype=np.complex64)
            
            indices = np.arange(self.k.size)
            
            for a in indices:
                kinetic_operator_half[a]=kinetic_diagonal_element_half(self,
                                                                       self.k[a])
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
    
    def first_step(self, x_wave_is):
        """
        ACHTUNG: Rechne mit Operatoren, dies ist im Algorithmus somit der ERSTE Befehl der anzuwenden ist!
        """
        k_wave_f = np.fft.fft(x_wave_is, norm='ortho')
        kin_k_wave = np.multiply(self.kinetic_operator_half, k_wave_f, dtype=np.complex64)
        x_wave2 = np.fft.ifft(kin_k_wave, norm='ortho')
        pot_x_wave = np.multiply(self.potential_operator, x_wave2, dtype=np.complex64)
        return np.fft.fft(pot_x_wave)
    
    def step(self, wave):
        """
        Anwendung nach final step, bis Abbruchkriterium erfuellt. 
        """
        kin_k_wave_s = np.multiply(self.kinetic_operator, wave, dtype=np.complex64)
        x_wave_s = np.fft.ifft(kin_k_wave_s, norm='ortho')
        pot_x_wave_s = np.multiply(self.potential_operator, x_wave_s, dtype=np.complex64)
        k_wave_s2 = np.fft.fft(pot_x_wave_s, norm='ortho')
        return np.multiply(self.kinetic_operator, k_wave_s2, dtype=np.complex64)
    
    def final_step(self, k_wave_fs):
        """
        Anwendung nachdem das Abbruchkriterium erfuellt wurde.
        """
        x_wave_fs = np.fft.ifft(k_wave_fs, norm='ortho')
        pot_x_wave_fs = np.multiply(self.potential_operator, x_wave_fs, dtype=np.complex64)
        k_wave_fs2 = np.fft.fft(pot_x_wave_fs, norm='ortho')
        kin_k_wave_fs = np.multiply(self.kinetic_operator_half, k_wave_fs2, dtype=np.complex64)
        return np.fft.ifft(kin_k_wave_fs)
