#Attention: This class uses dx = L/(N-1) with the total length of
#the potential and the number samplepoints N.
#It is different to the wikipedia formular dx=L/N.
#Argumentation: the length of the grid is L=(N-1)*dx

import numpy as np

class Fourier():
    def __init__(self, x):
        """
        Creates a grid in wavenumber-space.
        
        Transformation between position- and wavenumber-space.
        
        Parameters
        ----------
        func: Array
            function on discret grid
        x : Array
            grid in position space
            
        Attributes
        ----------
        self.k : Array
            grid of samplepoints in wavenumber-space
        """
        #set Parameters for 
        self.x = x
        self.dx = x[1]-x[0]
        self.L = x[-1]-x[0]
        self.N = x.size
        
        self.k0 = -np.pi/self.dx #smallest useable wavenumber
        self.dk = 2*np.pi/self.L
        self.k = np.arange(self.k0, self.k0 + (self.N)*self.dk, self.dk)

        self.norm = np.sqrt(2*np.pi*(self.N-1))/self.L
        
        #create trafo matrix
        def trafo_matrix_element(self, x, k):
            return np.exp(1j*k*x)/np.sqrt(self.N-1)
        
        def create_trafo_matrix(self):
            #create matrix with correct shape filled with zeros
            matrix = np.zeros(self.N**2).reshape((self.N,self.N))
            trafo_matrix = matrix + 1j*matrix
            
            #fill trafo_matrix with correct elements
            indices = np.arange(0, self.N, 1)
            
            for a in indices:
                for b in indices:
                    trafo_matrix[a][b] = trafo_matrix_element(self,
                                                              self.x[b],
                                                              self.k[a])
            print("shape of trafo matrix")
            print(trafo_matrix.shape)
            return trafo_matrix
        
        def get_trafo_matrix(self):
            return (np.asmatrix(create_trafo_matrix(self)))
        
        self.trafo_matrix = create_trafo_matrix(self)
        self.dagger_trafo_matrix = self.trafo_matrix.conjugate().T
        
    #define functions that perform the transformations
    def dft(self, x_func):
        
        var = np.dot(self.dagger_trafo_matrix, x_func)/self.norm

        return var.reshape(-1)
        
    def idft(self, k_func):
        return self.norm*np.dot(self.trafo_matrix, k_func).reshape(-1)
