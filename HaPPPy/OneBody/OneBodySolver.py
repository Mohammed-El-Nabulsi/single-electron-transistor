from scipy import constants
import matplotlib.pyplot as plt
import numpy as np
import h5py
from math import sqrt

class OneBodySolver:
    """ Solves the one body problem

    TODO: Add more documentation

    """
    # validating constants for unit calculation of hamilton matrix
    hbar = constants.hbar
    me = constants.m_e
    e = constants.e
    
    def __init__(self, l, n):
        """ The constructor.
        """
        self.l = l
        self.n = n
        print("Hello from the OneBody Solver, what can I do for you?")
        # creating array in user input length l and user input elements n and not changeable array for plot x-axis
        # creating transposed array a -> at
        a = np.linspace((-self.l / 2), (self.l / 2), self.n)
        aplot = np.linspace((-self.l / 2), (self.l / 2), self.n)
    

    


    def doCalculation(self):
        """ Gives the energy levels and energy eigen states
            If portenil is equal to 2 it uses A * x^2
        
        Returns:
            double.  The result

        """
        
        if(l <= 0):
            print("l not set. Abort!")
        
        exit()
    
    
    def calcualteHarmonocPotential(n, intersection, factor):
        
        
        # creating potential matrix with user input elements
        pot_mat = np.zeros((self.n, self.n))
        
        # function to build potential matrix out of user adjusted array
        def mat_build(a):
            i = 0
            while i < self.n:
                pot_mat[i, i] = a[i]
                i += 1
            return pot_mat
        
        for x in np.nditer(a, op_flags=['readwrite']):
            x[...] = factor * (x**2) + intersection
        mat_build(a)  # build pot_mat
        
        
        # creating body of the kinetic matrix with second derivate of the location
        kin_mat = np.zeros((n, n))

        i = 0
        while i < n:
            kin_mat[i, i] = 2
            i += 1

        i = 0
        while i < n-1:
            kin_mat[i, i+1] = kin_mat[i+1, i] = -1
            i += 1
        print(kin_mat)

        # unit system and calculation of final hamiltonian matrix ham_mat
        # factor 1000 for the unit system in order to reach meV
        # unit_pot = ((1/2)*me*(3.5174*(10**29))*((10**-9))**2)/e  # potential matrix in meV

        unit_kin = ((hbar**2)*1000)/(2*me*(10**-18)*e)
        print(unit_kin, "\n", unit_pot)  # control print for unit
        # dx for the derivate of the matrix
        dx = l/n
        # build the final hamilton matrix ham_mat
        ham_mat = unit_kin * kin_mat * (1/(dx*dx)) + unit_pot * pot_mat

        # calculate eigenvalues (stored in la) and eigenvectors (stored in v)
        la, v = np.linalg.eigh(ham_mat)

        

        
        return la, V
        
        
    def calcualteGaussPotential(absdfjksdf):


        ham = kin + pot
    
        V = eigh(V)

        return V
        
    

        return 42.5
