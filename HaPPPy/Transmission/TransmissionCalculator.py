from scipy.constants import codata
import numpy as np
from scipy.fftpack import fft, ifft
import matplotlib.pyplot as plt
import string

eV = codata.value("electron volt")

me   = codata.value("electron mass") * 1e15 # convert to pico grams            
hbar = codata.value("Planck constant over 2 pi") * 1e27  # convert to pico-compatible values


class TransmissionCalculator():
    """
    Calculates the probability of an electron tunneling through an potential barrier.
    
    Parameters
    ----------
    E : float
        Energy of the electron in meV
    dx : float
        Distance between two neigboured sample points in pm
    barrier : Array
        Potenital barrier to be tunneled with values in meV
    """
    def __init__(self, step_callback = None, value_callback = None,  _me = None, _hbar = None, package_wdh = None, disable_electron_potential_validation = None):
        global me
        global hbar

        self.step_callback = step_callback
        self.value_callback = value_callback
        self.mock_package_width = package_wdh
        self.disable_electron_potential_validation = disable_electron_potential_validation
        
        if (_me):
            me = _me

        if (_hbar):
            hbar = _hbar

    def validate_input(self, E, barrier, dx):
        error_msg_tmpl = string.Template("Validation Error: $message Details: $details")

        if (not hasattr(barrier, "__iter__") or np.array(barrier).ndim != 1):
            message = "The potential must be an array of one dimension"
            details = f"Dimension: {np.array(barrier).ndim}. Instance: {barrier.__class__.__name__}"
            raise ValueError(error_msg_tmpl.substitute(message = message, details = details))

        max_V = max(barrier)
        min_V = min(barrier)

        if (E < 0):
            message = "Electron energy must be greater than 0."
            details = f"Electron energy given: {E}."
            raise ValueError(error_msg_tmpl.substitute(message = message, details = details))

        if (E > max_V and not self.disable_electron_potential_validation):
            message = "Electron energy cannot be bigger than max potential value."
            details = f"Electron energy given: {E}. Max Potential: {max_V}."
            raise ValueError(error_msg_tmpl.substitute(message = message, details = details))

        if (max_V == 0 and min_V == 0 and not self.disable_electron_potential_validation):
            message = "Potential must contain values other than 0."
            details = f"Empty."
            raise ValueError(error_msg_tmpl.substitute(message = message, details = details))

        if (dx <= 0):
            message = "dx must be greater than 0"
            details = f"dx: {dx}."
            raise ValueError(error_msg_tmpl.substitute(message = message, details = details))

    
    def calculate_transmission(self, E, barrier, dx):
        self.validate_input(E, barrier, dx)
        
        E = E * eV * 1e-3 * 1e12
        barrier = barrier * eV * 1e-3 * 1e12
        # dx = dx * 1e-12 

        #input parameters
        self.E = E
        self.dx = dx
        self.barrier = barrier
        
        #build potential and helper objects
        self.N = barrier.size
        self.V = self.get_potential()
        self.trans_index = self.get_first_index_behind_barrier()
        
        #build position grid and fix sizes with respect to V if needed. This could be nessesary if N is odd
        self.x = self.get_position_grid()
        while self.x.size<self.V.size:
            self.x = np.append(self.x,[self.x[-1]+self.dx])
        while self.x.size>self.V.size:
            self.V = np.append(self.V,[0])
        
        #build wavenumbergrid depending on the now known length len(self.x)
        self.M  = self.V.size
        self.k_min = -np.pi/self.dx
        self.dk = 2*np.pi/(self.x[-1]-self.x[0])
        self.k = self.get_wavenumber_grid()
        
        max_k = self.get_wavenumber_grid()[-1]
        self.k0 = np.sqrt(2*me*self.E)/hbar
        max_E =   max_k**2 * hbar**2 / (2 * me)

        if (self.k0 >= max_k):
            error_msg_tmpl = string.Template("Validation Error: $message Details: $details")

            message = "The energy provided results in too much intertia. The algorithm cannot provide results for such energies."
            details = f"Energy: {E}. Allowed energy range: (0, {max_E / eV} eV]."
            raise ValueError(error_msg_tmpl.substitute(message = message, details = details))

        #choose parameters for the initial state (gaussian package)
        self.sigma = self.dx*self.N / 10 if self.mock_package_width is None else self.mock_package_width #initial width of the package
        # self.sigma = hbar / np.sqrt(2*E*me) if self.mock_package_width is None else self.mock_package_width #initial width of the package
        self.x0 = self.x[self.trans_index-self.N-5*int(self.sigma/self.dx+1)] #symmetrypoint in positions space
        
        self.v = hbar / (me * self.sigma)
        
        self.t = (self.dx * self.x.size) / self.v

        #build initial state (gaussian package)
        self.psi0 = (2/(np.pi*self.sigma**2))**(1/4) * np.exp(1j*self.k0*(self.x-self.x0)) * np.exp(-((self.x-self.x0)/self.sigma)**2)
        # self.psi0 = (self.sigma*np.sqrt(np.pi))**(-1/2)*np.exp(1j*self.x*self.k0-1/2*((self.x-self.x0)/self.sigma)**2)

        #chosse time between to two neigbourt time sample points (to be used in the split step algorithm)
        # steps = 1000  #  150 #  20**3 #  0.25 * 1 / self.k0 # 10**10
        self.dt = self.t / 10000  #  0.05 # self.t / (400 * steps)  #  150 #  20**3 #  0.25 * 1 / self.k0 # 10**10

        self.steps = self.N * 0.8
        
        #build Operators for the split-step-algorithm
        self.V_Op_half = np.exp(-(1j*self.dt/2)*(self.V/hbar))
        self.T_Op = np.exp(-(1j*hbar*np.multiply(self.k,self.k)*self.dt/(2*me)))
        
        if (self.value_callback != None):
            self.value_callback(self)

        #performing z steps and calculate the probability of presence behind the barrier
        self.psi_after_z_steps = self.perform_z_split_steps(z = self.steps)
        
        #calculate density of presence
        self.density_after_z_steps = self.get_density_of_probability(self.psi_after_z_steps)
        
        #calculate transmission probability after 1500 steps
        self.T_1500 = self.get_transmission(self.density_after_z_steps)
            
        psi = self.psi_after_z_steps
        
        psi = self.dx*self.psi0*np.exp(-1j*self.k_min*self.x)/np.sqrt(2*np.pi) #fourier scaling
        
        delta = 1
        
        T_stop_crit = self.T_1500
        
#        while delta > 10**(-1):
#            psi = self.step(psi)
#            
#            #stop criterion
#            psi_stop_crit = psi*np.sqrt(2*np.pi)*np.exp(1j*self.k_min*self.x)/self.dx #undo the fourier scaling 
#            psi_stop_crit_dens = self.get_density_of_presens(psi_stop_crit)
#            T_stop_crit_n = self.get_transmission(psi_stop_crit_dens)
#            
#            delta = np.abs(T_stop_crit_n-T_stop_crit)
#            
#            T_stop_crit = T_stop_crit_n
#            
#        print(T_stop_crit)
#            
#        psi = psi*np.sqrt(2*np.pi)*np.exp(1j*self.k_min*self.x)/self.dx #undo the fourier scaling   
#
#        
#        psi_abs_squared = self.get_density_of_presens(psi)
#        T = self.get_transmission(psi_abs_squared)

        T = self.T_1500
        
        return T
        
    def get_potential(self):
            zeros = np.zeros(25*self.N)
            V_ref = np.append(zeros,self.barrier) #part of the xgrid where the reflectet wave is living
            return np.append(V_ref,zeros)
    
    def get_first_index_behind_barrier(self):
        zeros = np.zeros(25*self.N)
        V_ref = np.append(zeros,self.barrier) #part of the xgrid where the reflectet wave is living
        return V_ref.size
    
    def get_position_grid(self):
        x = np.arange(-int(self.V.size/2)*self.dx,int(self.V.size/2)*self.dx,self.dx)
        return x
    
    def get_wavenumber_grid(self):
        return self.k_min+self.dk*np.arange(0,self.M,1)
    
    def step(self,psi):
        psi1 = np.array(np.multiply(psi,self.V_Op_half))
        psi2 = np.array(fft(psi1))
        psi3 = np.array(np.multiply(psi2,self.T_Op))
        psi4 = np.array(ifft(psi3))
        return np.array(np.multiply(psi4,self.V_Op_half))
    
    def perform_z_split_steps(self,z):
        psi = self.dx*self.psi0*np.exp(-1j*self.k_min*self.x)/np.sqrt(2*np.pi) #fourier scaling
        
        for n in np.arange(0,z,1):
            psi = self.step(psi)

            if (self.step_callback != None):
                psi_plot = psi*np.sqrt(2*np.pi)*np.exp(1j*self.k_min*self.x)/self.dx
                psi_plot = np.multiply(psi_plot,psi_plot.conj()).real
                
                self.step_callback(self, psi, psi_plot, self.x, self.dt * n,  n, z)
            
            n+=1

        return psi*np.sqrt(2*np.pi)*np.exp(1j*self.k_min*self.x)/self.dx #undo the fourier scaling
    
    def get_density_of_probability(self,wave):
        return np.multiply(wave,wave.conj()).real
    
    def get_transmission(self,psi2):
        """
        input: psi squared after sucsessful split step iteration
        """
        p_tr = self.dx*np.sum(psi2[self.trans_index:])
        return p_tr

##constants
# hbar = 1
# me = 1
#
##input
# E0 = 20
# E = E0
# dx = 0.1
# barrier = np.array(20+np.zeros(3000))

# testinstance = TransmissionCalculator(_me = 1, _hbar = 1)
# #
# ##print(testinstance.T_600)
# print(testinstance.calculate_transmission(E, barrier, dx))
