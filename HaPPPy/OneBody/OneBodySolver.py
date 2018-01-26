from scipy import constants
import matplotlib.pyplot as plt
import numpy as np
import h5py
from math import sqrt
from HaPPPy.TwoBody.OneParticleLoader import SpectrumData

__docformat__ = 'reStructuredText'

class OneBodySolver:
    """ 
    Global class constants are planck's constant :math:`\hbar`, electron mass :math:`m_e`, and electron charge :math:`e`.
    
    :param l: length of the quantum dot's potential in :math:`nm`.    
    :type l: int or float
    
    :param n: number of grid points used for the numerical approximation.
    :type n: int
    
    :param m: number of eigenvalues and corresponding eigenvectors which are saved to a hdf5 file
    :type m: int
    
    Two equal global numpy arrays (``a`` and ``a_axis``) are created from :math:`-\\frac{l}{2}` to :math:`\\frac{l}{2}` with *n* grid points. 
            
        - ``a_axis`` serves as holder for the x values
        - ``a`` serves as holder for the (in each potential to calculate) y values
            
    .. figure::  _static/linspace.jpg
        :align:   center
    
    
    """
    # validating constants for unit calculation of Hamilton matrix
    hbar = constants.hbar
    me = constants.m_e
    e = constants.e

    def __init__(self, l, n, m):

        self.l = l
        self.n = n
        self.m = m
#        print("Hello from the OneBody Solver, what can I do for you?")

        # creating array in user input length l and user input elements n and not changeable array for plot x-axis
        
        self.a = np.linspace((-self.l / 2), (self.l / 2), self.n)
        self.a_axis = np.linspace((-self.l / 2), (self.l / 2), self.n)
    
    # possible input factor but we defined factor already as unit_pot
    def calcualteHarmonicPotential(self, intersection, omega=1.875537359*10**13):
        """ 
        Calculates the OneBody problem for harmonic potential.
        
        :param intersection: y position of the harmonic potential in :math:`meV`.
        :type intersection: int or float
        
        :param omega: oscillation frequency in :math:`\\frac{1}{s}`. Default value is :math:`1,875537359*10^{13}`.
        :type omega: float
        
        :return: la, v_norm , Info, kin_mat, pot_mat, self.a_axis, la_l, v_norm_l
        :rtype: np.array,matrix, np.array, matrix, matrix, np.array, np.array, matrix


        **Potential term:**
  
                :math:`\\widehat{V} =` ``unit_pot`` :math:`*\\begin{pmatrix} V(x_0)&0&0&0\\\ 0 &V(x_2)&0&0\\\ 0 &0&\\ddots&0\\\ 0 &0&0&V(x_{n-1})\\end{pmatrix}`
                
                :math:`\\widehat{V} = \\frac{m_e\omega^2}{2}*\\underline{V}(x)`
        
            Since the unit of the potential is currently Joule, ``unit_pot`` has to be multiplied by :math:`\\frac{1000}{e}` in order to reach :math:`meV`. 
            The program work best with a ``unit_pot`` close to 1.0. Therefore the default value of omega is initialized. 
        
                ``unit_pot`` = :math:`\\frac{m_e\omega^2}{2}*\\frac{1000}{e}`
                
                :math:`V(x) = x^2+intersection`
            
            Values for x are taken from ``a``. The vector is then overwritten and put on the diagonal of the before initialized ``pot_mat`` 
            as shown in the code example.
            
                .. figure::  _static/pot_mat.jpg
                    :align:   center
                    
            A plot of ``a`` vs. ``a_axis`` is shown in figure 2.

                .. figure::  _static/harmonic_potential.jpg
                    :align:   center
            
            **Figure 2:** Plot of harmonic potential in arbitrary units.
            
        **Kinetic term:**
            
                :math:`\\widehat{K} =` ``unit_kin`` :math:`* \\begin{pmatrix} -2&1&0&0&0\\\ 1 &-2&1&0&0\\\ 0 &1&-2&1&0\\\ 0 &0&\\ddots&\\ddots&\\ddots\\\ 0 &0&0&1&-2 \\end{pmatrix}`
                
                :math:`\\widehat{K} = -\\frac{\\hbar^2}{2m} * \\frac{1}{\\Delta x^2} * \\underline{K}`
                
                :math:`\\underline{K}` is saved as the matrix ``kin_mat``. Since the unit of the potential is currently Joule, ``unit_kin`` has 
                to be multiplied by :math:`\\frac{1000}{e}` in order to reach :math:`meV`.
                
                
                ``unit_kin`` :math:`= -\\frac{\\hbar^2}{2m} * \\frac{1}{\\Delta x^2}*\\frac{1000}{e}`
                
        
        **Hamilton matrix:**
        
                :math:`\\widehat{H} = \\widehat{K} + \\widehat{V}`
                
                ``ham_mat`` = ``unit_kin`` * ``kin_mat`` + ``unit_pot`` * ``pot_mat``
                
        The eigenvalues and eigenvectors are calculated with the numpy library ``linalg.eigh()``. Eigenvalues are stored as a vector
        in ``la_l`` and eigenvectors are stored as a matrix in ``v``. The matrix ``v_norm_l`` is a matrix with the normalized eigenvectors. "_l" stands for "long".
        To the next group only ``m`` eigenvalues and -vectors as ``la``, ``v_norm`` are given.  The code is shown below:
        
                .. figure::  _static/linalg_eigh.jpg
                    :align:   center
        
        The probability density of the first 3 wave functions is shown in figure 3. The absolute square of the first 3 eigenvectors vs. ``a_axis`` is plotted.
        
            .. figure::  _static/harmonic_wavefkt.jpg
                :align:   center

        **Figure 3:** Probability density of the first 3 wave functions to corresponding energy eigenvalues (``l`` = 100, ``n`` = 1000). The density is given arbitrary units.
        
        ``Info`` is an array with the used module parameter saved as strings. It saves the grid points ``n``, length of potential ``l``, number of eigenvalues -vectors ``m``,
        not used potentials as ``False`` and the used potential as ``True``.
             
        """
#        print("\nYou've chosen HarmonicPotential.")
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

        # unit system and calculation of final Hamiltonian matrix ham_mat
        # factor 1000 for the unit system in order to reach meV
        unit_pot = ((1/2)*self.me*1000*omega**2*((10**-9))**2)/self.e  # potential matrix in meV

        unit_kin = ((self.hbar**2)*1000)/(2*self.me*(10**-18)*self.e)
        #print(unit_kin, "\n", unit_pot)  # control print for unit
        # dx for the derivate of the matrix
        dx = self.l/self.n
        # build the final Hamilton matrix ham_mat
        ham_mat = unit_kin * kin_mat * (1/(dx*dx)) + unit_pot * pot_mat

        # calculate eigenvalues (stored in la) and eigenvectors (stored in v)
        la_l, v = np.linalg.eigh(ham_mat)
        
        # normalize the eigenvectors
        norm = 0.0
        for i in range(self.n):
            norm += (v[i,0]*v[i,0]) * dx
        sqrt_norm = sqrt(norm)

        v_norm_l = v / sqrt_norm # v_norm_l is now norm matrix of eigenvectors
        
        # return m eigenvalues/eigenvector
        la = la_l[:self.m]
        v_norm = v_norm_l[:,:self.m]

        # printing eigenvalues and eigenvectors
        # as option for debugging
        #for i in range(10):
        #    print("Eigenvalue:\t\t ", la[i])
        
        # create norm of the eigenvectors
        
        # test if eigenvectors are normed here 25. 
        #norm = 0.0
        #for i in range(self.n):
        #    norm += (v[i,25]*v[i,25]) * dx
        #print(norm)
        
        # list of tuples for calculation information
        Info = np.array([["n-grids point" ,str(self.n)],["l-length of potential",str(self.l)],
                ["number of eigenvalues -vectors",str(self.m)],["HarmonicPotential",str(True)], \
        ["GaussPotential",str(False)],["BoxPotential",str(False)]]).astype('S9')

        pot_mat = unit_pot * pot_mat

        return la, v_norm , Info, kin_mat, pot_mat, self.a_axis, la_l, v_norm_l


    def calculateBoxPotential(self, intersection):
        """ 
        Calculates the OneBody problem for box potential. The calculation is done according to the description for the **HarmonicPotential**.
        
        :param intersection: Ground level of potential.
        :type intersection: int or float
        
        :return: la, v_norm , Info, kin_mat, pot_mat, la_l, v_norm_l
        :rtype: np.array,matrix, np.array, matrix, matrix, np.array, matrix
        
        The potential :math:`V(x)` has the unit ``unit_pot`` = 1 in :math:`meV`. The potential is plotted in figure 4.
        
            :math:`V(x) = intersection`
        
            .. figure::  _static/box.jpg
                :align:   center

        **Figure 4:** Plot of the box potential in arbitrary units.
        
        The probability density of the first 3 wave functions are shown in figure 5. The absolute square of the first 3 eigenvectors vs. ``a_axis`` is plotted.

            .. figure::  _static/BoxPot_wavefkt.jpg
                :align:   center

        **Figure 5:** Probability density of the first 3 wave functions to corresponding energy eigenvalues (``l`` = 100, ``n`` = 1000). The density is given arbitrary units.
        
        """
#        print("\nYou've chosen BoxPotential.")

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

        # unit system and calculation of final Hamiltonian matrix ham_mat
        # factor 1000 for the unit system in order to reach meV
        #unit_pot = ((1/2)*self.me*(3.5174*(10**29))*((10**-9))**2)/self.e  # potential matrix in meV
        unit_pot = 1

        unit_kin = ((self.hbar**2)*1000)/(2*self.me*(10**-18)*self.e)
        #print(unit_kin, "\n", unit_pot)  # control print for unit
        # dx for the derivate of the matrix
        dx = self.l/self.n
        # build the final Hamilton matrix ham_mat
        ham_mat = unit_kin * kin_mat * (1/(dx*dx)) + unit_pot * pot_mat

        # calculate eigenvalues (stored in la) and eigenvectors (stored in v)
        la_l, v = np.linalg.eigh(ham_mat)

        # create norm of the eigenvectors
        norm = 0.0
        for i in range(self.n):
            norm += (v[i,0]*v[i,0]) * dx

        sqrt_norm = sqrt(norm)

        v_norm_l = v / sqrt_norm # v is now norm matrix of eigenvectors

        # test if eigenvectors are normed here 25. 
        #norm = 0.0
        #for i in range(self.n):
        #    norm += (v[i,25]*v[i,25]) * dx
        #print(norm)
        
        # list of tuples for calculation information
        Info = np.array([["n-grids point" ,str(self.n)],["l-length of potential",str(self.l)],
                ["number of eigenvalues -vectors",str(self.m)],["HarmonicPotential",str(True)], \
        ["GaussPotential",str(False)],["BoxPotential",str(False)]]).astype('S9')
        
        # return m eigenvalues/eigenvector
        la = la_l[:self.m]
        v_norm = v_norm_l[:,:self.m]
        return la, v_norm , Info, kin_mat, pot_mat, la_l, v_norm_l

    def calcualteGaussPotential(self, A, sigma):
        """ 
        Calculates the OneBody problem for gauss potential. The calculation is done according to the description for the **HarmonicPotential**.
        
        :param A: Unit of the potential in :math:`meV`.
        :type A: int or float
        
        :param sigma: Width of the gauss curve.
        :type A: int or float
        
        :return: la, v_norm, Info, kin_mat, pot_mat, la_l, v_norm_l
        :rtype: np.array,matrix, np.array, matrix, matrix, np.array, matrix
        
        The potential :math:`V(x)` has the unit ``unit_pot`` = ``A`` in :math:`meV`. The potential is plotted in figure 6.
        
            :math:`V(x) = -A*\\exp{({-\\frac{x}{2 \\sigma}})^2}`
        
            .. figure::  _static/gauss.jpg
                :align:   center

        **Figure 6:** Plot of the gauss potential in arbitrary units (``A`` = 1, ``sigma`` = 5).
        
        The probability density of the first 3 wave functions are shown in figure 7. The absolute square of the first 3 eigenvectors vs. ``a_axis`` is plotted.

            .. figure::  _static/gauss_wavefkt.jpg
                :align:   center

        **Figure 5:** Probability density of the first 3 wave functions to corresponding energy eigenvalues (``l`` = 100, ``n`` = 1000). The density is given arbitrary units.
        
        """
        #        print("\nYou have chosen GaussPotential!")


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
            x[...] = -A*(np.exp(-np.power((x) / sigma, 2.) / 2.))
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
        #        print(kin_mat)

        unit_pot = A
        unit_kin = ((self.hbar**2)*1000)/(2*self.me*(10**-18)*self.e)
        #        print(unit_kin, "\n", unit_pot)  # control print for unit
        # dx for the derivate of the matrix
        dx = self.l/self.n
        # build the final Hamilton matrix ham_mat
        ham_mat = unit_kin * kin_mat * (1/(dx*dx)) + unit_pot * pot_mat
        
        # calculate eigenvalues (stored in la) and eigenvectors (stored in v)
        la_l, v = np.linalg.eigh(ham_mat)

        # printing eigenvalues and eigenvectors
        # as option for debugging
        #for i in range(10):
        #    print("Eigenvalue:\t\t ", la[i])
        
        # create norm of the eigenvectors
        norm = 0.0
        for i in range(self.n):
            norm += (v[i,0]*v[i,0]) * dx

        sqrt_norm = sqrt(norm)
        #        print(sqrt_norm)

        v_norm_l = v / sqrt_norm # v is now norm matrix of eigenvectors

        # test if eigenvectors are normed here 25. 
        #norm = 0.0
        #for i in range(self.n):
        #    norm += (v[i,25]*v[i,25]) * dx
        #print(norm)
        
        # list of tuples for calculation information
        Info = np.array([["n-grids point" ,str(self.n)],["l-length of potential",str(self.l)],
                        ["number of eigenvalues -vectors",str(self.m)], ["HarmonicPotential",str(False)], \
        ["GaussPotential",str(True)],["BoxPotential",str(False)]]).astype('S9')
        # return m eigenvalues/eigenvector
        la = la_l[:self.m]
        v_norm = v_norm_l[:,:self.m]
        return la, v_norm, Info, kin_mat, pot_mat, la_l, v_norm_l

        
    def exportData(self, la, v_norm, info, path="data_group1"):
        """
        Saves the calculated eigenvalues and -vectors as well as the information array Info to a hdf5 format file.
        
        :param la: Array holding the eigenvalues.
        :type la: numpy.array
        
        :param v_norm: Matrix holding the normalized eigenvectors.
        :type v_norm: matrix
        
        :param info: Array holding information about the chosen parameters.
        :type info: numpy.array holding strings
        
        :param path: Name of the hdf5 file with default value "data_group1".
        :type path: string
        
        To secure the correct save and load procedure, the hdf5 file save system of the group 2 is used. Each value holding parameter is saved in its own data set.
        For any further information, please read the **Two Body Module** documentation.
        """
        
        Data = SpectrumData()
        Data.init(path, len(la), len(la), info=info)
        Data.energies = la
        Data.waves = v_norm
        Data.close()





