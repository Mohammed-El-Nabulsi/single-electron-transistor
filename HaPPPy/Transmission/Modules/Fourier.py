import numpy as np
import cmath

class Fourier():
    #inverse transformationmatrixelement
    def ITMElement(self, arg,newArg, N):
    
        element = cmath.exp(1j*newArg*arg)/np.sqrt(N)
    
        return element
    
    #FT transform to new arguments
    def FT(self, xFunction, x0, L):
        '''
        Returns DFT of an given complex array in wavenumber-space, centured around the origin.
        
        Parameters
        ----------
        xFunction : Array
            A sequence of wavefunction values in position space 
        x0 : float
            lowest sample point argument in position space         
        L : float
            range of the sample points in position space
            
        Returns
        -------
        kFunction : array, shape(len(xFunction))
            Array containing the value of the DFT for each desired sample point in wavenumber space
        '''   
        
        N =  xFunction.size #Number of Arguments, needs to be odd for centering k's around 0
    
        dx = L/N
        dk = 2*cmath.pi/L
        k0 = -(int(N-1))/2*dk #lowest new arguments so that k's centred around k = 0
        
        #create transformationMatrices
        indices = np.arange(0, N, 1) # indices of transformation matrices
        
        arguments = np.array(np.arange(x0, x0 + L, dx)) #samplepoints of the function
        
        
        newArguments = np.array(np.arange(k0,k0+N*dk, dk))
        
        ITM = np.zeros((N,N), dtype=np.complex64) #InverseTransformationMatrix
        
        for n in indices: #rowindeces
            for m in indices: #columindices
                ITM[n][m] = self.ITMElement(arguments[n],newArguments[m], N)
        
        TM = np.transpose(ITM.conjugate()) #transformationsmatrix    
        
        kFunction = L*np.dot(TM,xFunction)/np.sqrt(2*np.pi*N)
        
        return kFunction
    
    #IFT:retransform to arguments
    def IFT(self, kFunction, x0, L):
        '''
        Returns IDFT of an given complex array, with lowest samplepoint x0.
        
        Parameters
        ----------
        kFunction : Array
            A sequence of wavefunction values in in wavenumber-space
        x0 : float
            lowest sample point argument in position space         
        L : float
            range of the sample points in position space
            
        Returns
        -------
        xFunction : array, shape(len(xFunction))
            Array containing the value of the IDFT for each desired sample point in position space
        '''   
        N =  kFunction.size #Number of Arguments, needs to be odd for centering k's around 0
    
        dx = L/N
        dk = 2*cmath.pi/L
        k0 = -(int(N-1))/2*dk #lowest new arguments so that k's centred around k = 0
        
        #create transformationMatrices
        indices = np.arange(0, N, 1) # indices of transformation matrices
        
        arguments = np.array(np.arange(x0, x0 + L, dx)) #samplepoints of the function
        
        
        newArguments = np.array(np.arange(k0,k0+N*dk, dk))
        
        ITM = np.zeros((N,N), dtype=np.complex64) #InverseTransformationMatrix
        
        for n in indices: #rowindeces
            for m in indices: #columindices
                ITM[n][m] = self.ITMElement(arguments[n],newArguments[m], N)
        
        xFunction = np.sqrt(2*np.pi*N)*np.dot(ITM,kFunction)/L
        
        return xFunction
        
    def waveNumbers(self, xfunction, x0, L):
         '''
        Calculates wavenumbers centured around 0.
        
        Parameters
        ----------
        xfunction : Array
            A sequence of wavefunction values in in position-space
        x0 : float
            lowest sample point argument in position space         
        L : float
            range of the sample points in position space
            
        Returns
        -------
        waveNumbers : Array, shape(len(xFunction))
        '''   
        N =  xfunction.size

        dk = 2*np.pi/L
        k0 = -(N-1)/2*dk #lowest new arguments so that k's centred around k = 0
        
        waveNumers = np.array(np.arange(k0,k0+N*dk, dk))
        return waveNumers
