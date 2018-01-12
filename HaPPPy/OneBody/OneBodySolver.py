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
    a = ""
	#aplot = ""
    def __init__(self, l, n):
        """ The constructor.
        """
        self.l = l
        self.n = n
        print("Hello from the OneBody Solver, what can I do for you?")
        # creating array in user input length l and user input elements n and not changeable array for plot x-axis
        # creating transposed array a -> at
        self.a = np.linspace((-self.l / 2), (self.l / 2), self.n)
        #aplot = np.linspace((-self.l / 2), (self.l / 2), self.n)
        print(self.a)
        print("der geht")
		return self.a
    


    def doCalculation(self):
        """ Gives the energy levels and energy eigen states
            If portenil is equal to 2 it uses A * x^2
        
        Returns:
            double.  The result

        """
        
        if(l <= 0):
            print("l not set. Abort!")
        
        exit()
    
    # possible input factor but we defined factor already as unit_pot
    def calcualteHarmonocPotential(self, intersection):
        """ Gives the energy levels and energy eigenvalues and eigenvectors
            It uses (x^2 + m) as the potential
        
        Returns:
            double.  The result

        """
        print("Hallo Welt")
		
        # creating potential matrix with user input n elements
        pot_mat = np.zeros((self.n, self.n))
        
        # function to build potential matrix out of user adjusted array
        def mat_build(a):
            i = 0
            while i < self.n:
                pot_mat[i, i] = a[i]
                i += 1
            return pot_mat
        
        for x in np.nditer(self.a, op_flags=['readwrite']):
            x[...] =(x**2) + self.intersection
        mat_build(self.a)  # build pot_mat
        print(pot_mat)
        
        # creating body of kinetic matrix with second derivate of the location
        kin_mat = np.zeros((self.n, self.n))

        i = 0
        while i < self.n:
            kin_mat[i, i] = 2
            i += 1
        i = 0
        while i < n-1:
            kin_mat[i, i+1] = kin_mat[i+1, i] = -1
            i += 1
        print(kin_mat)

        # unit system and calculation of final hamiltonian matrix ham_mat
        # factor 1000 for the unit system in order to reach meV
        unit_pot = ((1/2)*me*(3.5174*(10**29))*((10**-9))**2)/e  # potential matrix in meV

        unit_kin = ((hbar**2)*1000)/(2*me*(10**-18)*e)
        print(unit_kin, "\n", unit_pot)  # control print for unit
        # dx for the derivate of the matrix
        dx = self.l/self.n
        # build the final hamilton matrix ham_mat
        ham_mat = unit_kin * kin_mat * (1/(dx*dx)) + unit_pot * pot_mat

        # calculate eigenvalues (stored in la) and eigenvectors (stored in v)
        la, v = np.linalg.eigh(ham_mat)

		# printing eingenvalues and eigenvectors
		# as option for debugging
		for i in range(10):
			print("Eigenvalue:\t\t ", la[i])
		
		# creat norm of the eigenvectors
		norm = 0.0
		for i in range(self.n):
			norm += (v[i,0]*v[i,0]) * dx

		sqrt_norm = sqrt(norm)
		print(sqrt_norm)

		v_norm = v / sqrt_norm # v is now norm matrix of eigenvectors

        # test if eigenvectors are normed here 25. 
		#norm = 0.0
		#for i in range(self.n):
		#	norm += (v[i,25]*v[i,25]) * dx
		#print(norm)
		
		# list of tuples for calculation information
		Info = np.array([("n-grids point" ,n),("l-lenght of potential",l),("HarmonocPotential",True), \
        ("GaussPotential",False),("BoxPotential",False)])
		
        return la, v_norm , Info
        
        
    def calcualteGaussPotential(absdfjksdf):


        ham = kin + pot
    
        V = eigh(V)

        return V
        
    

        return 42.5
	
	def exportData(self, self.la, self.v_norm):
		# export data as hdf5 file in two sets
		dataFile = h5py.File("data_group1.hdf5", "w")
		# dataSet_calcInfo = dataFile.create_dataset("Input Information", data=self.Info)
		dataSet_eigvalues = dataFile.create_dataset("eigenvalues_group1", data=self.la)  
		dataSet_eigvectors = dataFile.create_dataset("eigenvectors_group1", data=self.v_norm)
		dataFile.close()
