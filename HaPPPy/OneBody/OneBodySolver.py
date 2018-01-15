from scipy import constants
import matplotlib.pyplot as plt
import numpy as np
import h5py
from math import sqrt

class OneBodySolver:
    """ Solves the one body problem
    class has to receive via user input: l, n
        --> l, length of the potential in nm (best at 100 nm)
        --> n, number of grid points for numerical calc.(best at 1000 - 3000)

        Solves the one body problem with 3 different potentials as 3 functions:
            call functions with: .calculate*(), 
            give the parameters stated in parantheses behind the functions name
                --> BoxPotential(intersection)
                --> HarmonicPotential(intersection)
                --> GaussPotential(A, sigma)
                       intersection == y-position of potential
                       A == unit_kin in meV (best at 1 meV)
                       sigma == width of gauss function (best at 5)

        Data stored in Hdf5 format:
		    call function with : .exportData()
                --> exportData(la, v_norm, info)

        validating constants as class variables
            --> hbar, planck´s constant
            --> me, mass electron
            --> e, charge of electron
    """
    # validating constants for unit calculation of Hamilton matrix
    hbar = constants.hbar
    me = constants.m_e
    e = constants.e

    def __init__(self, l, n):
        """ The constructor receives l, n
        and creates an array from (-l/2) - (l/2) with n elements
        array is used as x-values for potential (Harmonic, Gauss, Box)
        """
        self.l = l
        self.n = n
        print("Hello from the OneBody Solver, what can I do for you?")
        
        # creating array in user input length l and user input elements n and not changeable array for plot x-axis
        self.a = np.linspace((-self.l / 2), (self.l / 2), self.n)

        
        


    def calculateBoxPotential(self, intersection):
        """ Gives energy eigenvalues and eigenvectors and Info
            needs intersection as input:
                --> uses y = intersection, as the potential
            eigenvalues:
                --> stored in la
            standardized eigenvectors: 
                --> stored in v_norm as ordered matrix, columns are eigenvector
            Info:
                --> list of list with user input information in each sub-list
                    choosen potential stored as True 
                    not choosen potential stored as False
        
        Returns:
            la, v_norm, Info

        """
        print("\nYou've choosen BoxPotential.")
		
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
            x[...] = intersection
        mat_build(self.a)  # build pot_mat
        #print(pot_mat)
        
        # creating body of kinetic matrix with second derivate of the location
        kin_mat = np.zeros((self.n, self.n))

        i = 0
        while i < self.n:
            kin_mat[i, i] = 2
            i += 1
        i = 0
        while i < self.n-1:
            kin_mat[i, i+1] = kin_mat[i+1, i] = -1
            i += 1
        #print(kin_mat)

        # unit system and calculation of final hamiltonian matrix ham_mat
        # factor 1000 for the unit system in order to reach meV
        #unit_pot = ((1/2)*self.me*(3.5174*(10**29))*((10**-9))**2)/self.e  # potential matrix in meV
        unit_pot = 1
		
        unit_kin = ((self.hbar**2)*1000)/(2*self.me*(10**-18)*self.e)
        #print(unit_kin, "\n", unit_pot)  # control print for unit
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
        
        # create norm of the eigenvectors
        norm = 0.0
        for i in range(self.n):
            norm += (v[i,0]*v[i,0]) * dx

        sqrt_norm = sqrt(norm)
        #print(sqrt_norm)

        v_norm = v / sqrt_norm # v is now norm matrix of eigenvectors

        # test if eigenvectors are normed here 25. 
        #norm = 0.0
        #for i in range(self.n):
        #    norm += (v[i,25]*v[i,25]) * dx
        #print(norm)
        
        # list of tuples for calculation information
        Info = np.array([["n-grids point" ,str(self.n)],["l-lenght of potential",str(self.l)],["HarmonocPotential",str(True)], \
        ["GaussPotential",str(False)],["BoxPotential",str(False)]])
        
        return la, v_norm , Info
			
    # possible input factor but we defined factor already as unit_pot
    def calcualteHarmonicPotential(self, intersection):
        """ Gives energy eigenvalues and eigenvectors and Info
            needs intersection as input:
                --> uses y = x^2 + intersection, as the potential
            eigenvalues:
                --> stored in la
            standardized eigenvectors:
                --> stored in v_norm as ordered matrix, columns are eigenvector
            Info:
                --> list of list with user input information in each sub-list
                    choosen potential stored as True 
                    not choosen potential stored as False
        
        Returns:
            la, v_norm, Info

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
            x[...] =(x**2) + intersection
        mat_build(self.a)  # build pot_mat
        #print(pot_mat)
        
        # creating body of kinetic matrix with second derivate of the location
        kin_mat = np.zeros((self.n, self.n))

        i = 0
        while i < self.n:
            kin_mat[i, i] = 2
            i += 1
        i = 0
        while i < self.n-1:
            kin_mat[i, i+1] = kin_mat[i+1, i] = -1
            i += 1
        #print(kin_mat)

        # unit system and calculation of final hamiltonian matrix ham_mat
        # factor 1000 for the unit system in order to reach meV
        unit_pot = ((1/2)*self.me*(3.5174*(10**29))*((10**-9))**2)/self.e  # potential matrix in meV

        unit_kin = ((self.hbar**2)*1000)/(2*self.me*(10**-18)*self.e)
        #print(unit_kin, "\n", unit_pot)  # control print for unit
        # dx for the derivate of the matrix
        dx = self.l/self.n
        # build the final hamilton matrix ham_mat
        ham_mat = unit_kin * kin_mat * (1/(dx*dx)) + unit_pot * pot_mat

        # calculate eigenvalues (stored in la) and eigenvectors (stored in v)
        la, v = np.linalg.eigh(ham_mat)

        # printing eingenvalues and eigenvectors
        # as option for debugging
        #for i in range(10):
        #    print("Eigenvalue:\t\t ", la[i])
        
        # create norm of the eigenvectors
        norm = 0.0
        for i in range(self.n):
            norm += (v[i,0]*v[i,0]) * dx

        sqrt_norm = sqrt(norm)
        print(sqrt_norm)

        v_norm = v / sqrt_norm # v is now norm matrix of eigenvectors

        # test if eigenvectors are normed here 25. 
        #norm = 0.0
        #for i in range(self.n):
        #    norm += (v[i,25]*v[i,25]) * dx
        #print(norm)
        
        # list of tuples for calculation information
        Info = np.array([["n-grids point" ,str(self.n)],["l-lenght of potential",str(self.l)],["HarmonocPotential",str(True)], \
        ["GaussPotential",str(False)],["BoxPotential",str(False)]])
        
        return la, v_norm , Info
        
        
    def calcualteGaussPotential(self, A, sigma):
        """ Gives energy eigenvalues and eigenvectors and Info
        needs intersection as input
            --> uses x^2 + intersection as the potential
        eigenvalues
            --> stored in la
        standardized eigenvectors 
            --> stored in v_norm as ordered matrix, columns are eigenvector
        Info
            --> list of list with user input information in each sub-list
                choosen potential stored as True 
                not choosen potential stored as False

        Returns:
            la, v_norm, Info

        """
        print("\nYou have chosen the Gauss potential!")
        

        # creating potential matrix with user input n elements
        pot_mat = np.zeros((self.n,self.n))
        
        # function to build potential matrix out of user adjusted array
        def mat_build(a):
            i = 0
            while i < self.n:
                pot_mat[i, i] = a[i]
                i += 1
            return pot_mat

        for x in np.nditer(self.a, op_flags=['readwrite']):
            x[...] = A*(-np.exp(-np.power((x) / sigma, 2.) / 2.))
        mat_build(self.a)  # build pot_mat


        # creating body of kinetic matrix with second derivate of the location
        kin_mat = np.zeros((self.n, self.n))

        i = 0
        while i < self.n:
            kin_mat[i, i] = 2
            i += 1
        i = 0
        while i < self.n-1:
            kin_mat[i, i+1] = kin_mat[i+1, i] = -1
            i += 1
        print(kin_mat)

        unit_pot = A
        unit_kin = ((self.hbar**2)*1000)/(2*self.me*(10**-18)*self.e)
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
        #    norm += (v[i,25]*v[i,25]) * dx
        #print(norm)
        
        # list of tuples for calculation information
        Info = np.array([["n-grids point" ,str(self.n)],["l-lenght of potential",str(self.l)],["HarmonocPotential",str(False)], \
        ["GaussPotential",str(True)],["BoxPotential",str(False)]])
        
        return la, v_norm, Info


 


    
    def exportData(self, la, v_norm, info):
        # export data as hdf5 file in three sets
        dataFile = h5py.File("data_group1.hdf5", "w")
        dataSet_calcInfo = dataFile.create_dataset("Input Information", data=info)
        dataSet_eigvalues = dataFile.create_dataset("eigenvalues_group1", data=la)  
        dataSet_eigvectors = dataFile.create_dataset("eigenvectors_group1", data=v_norm)
        dataFile.close()

