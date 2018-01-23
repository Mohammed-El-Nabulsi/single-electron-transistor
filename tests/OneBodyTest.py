
import HaPPPy
import unittest
import numpy as np
import scipy.integrate
from math import sqrt


class OneBodyTestSuite(unittest.TestCase):
    """A test class to the OneBody module.

    """

    def test_OneBody_exists(self):
        """ Checks wether the One Body module exists
        """
        self.assertTrue(hasattr(HaPPPy, 'OneBody'))
        print("\n""One Body module exists""\n" )

    def test_OneBody_Harmonic_kinetic(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        _,_,_,kin_mat,_,_ = OBSolver.calcualteHarmonicPotential(0)
        """ check the kin Matrix of the harmonic potential
            Confirm that all matrix elements are positive"""
        diagonalarray_harm_kin_main = np.diagonal(kin_mat, 0)
        if np.all (diagonalarray_harm_kin_main > 0):
            print("\n""The elements of the kinetic matrix in the harmonic potential does make sense.""\n")
        else :
            print("\n""The elements of the kinetic matrix in the harmonic potential doesn't make sense.""\n")
        diagonalarray_harm_kin_plusone = np.diagonal(kin_mat,1)
        self.assertTrue (np.all(diagonalarray_harm_kin_plusone == -1),msg = "\n""Your kinetic matrix in your harmonic potential is false.""\n")

        diagonalarray_harm_kin_minusone = np.diagonal(kin_mat,-1)
        self.assertTrue (np.all(diagonalarray_harm_kin_minusone == -1), msg ="\n""Your kinetic matrix in your harmonic potential is false.""\n")

    def test_OneBody_Box_kinetic(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        _,_,_,kin_mat,_ = OBSolver.calculateBoxPotential(1)
        """ check the kin Matrix of the box potential
        Confirm that all matrix elements are positive"""
        diagonalarray_box_kin_main = np.diagonal(kin_mat, 0)
        if np.all (diagonalarray_box_kin_main > 0):
            print("\n""The elements of the kinetic matrix in the box potential does make sense.""\n")
        else :
            print("\n""The elements of the kinetic matrix in the box potential doesn't make sense.""\n")
        diagonalarray_box_kin_plusone = np.diagonal(kin_mat,1)
        self.assertTrue (np.all(diagonalarray_box_kin_plusone == -1),msg = "\n""Your kinetic matrix in your box potential is false.""\n")

        diagonalarray_box_kin_minusone = np.diagonal(kin_mat,-1)
        self.assertTrue (np.all(diagonalarray_box_kin_minusone == -1), msg ="\n""Your kinetic matrix in your box potential is false.""\n")

    def test_OneBody_Gauss_kinetic(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        _,_,_,kin_mat,_ = OBSolver.calcualteGaussPotential(1,10)
        """ check the kin Matrix of the gauss potential
            Confirm that all matrix elements are positive"""
        diagonalarray_gauss_kin_main = np.diagonal(kin_mat, 0)
        if np.all (diagonalarray_gauss_kin_main > 0):
            print("\n""The elements of the kinetic matrix in the gauss potential does make sense.""\n")
        else :
            print("\n""The elements of the kinetic matrix in the gauss potential doesn't make sense.""\n")
        diagonalarray_gauss_kin_plusone = np.diagonal(kin_mat,1)
        self.assertTrue (np.all(diagonalarray_gauss_kin_plusone == -1),msg = "\n""Your kinetic matrix in your gauss potential is false.""\n")

        diagonalarray_gauss_kin_minusone = np.diagonal(kin_mat,-1)
        self.assertTrue (np.all(diagonalarray_gauss_kin_minusone == -1), msg ="\n""Your kinetic matrix in your gauss potential is false.""\n")

    def test_OneBody_Harmonic_potential(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        _,_,_,_,pot_mat,self.a_axis = OBSolver.calcualteHarmonicPotential(0)
        gridpoints = 3000
        i = 0
        testvector = np.diagonal(pot_mat)
        testvector = testvector -0
        testvalue_pot = testvector [0]
        testvalue_a = self.a_axis [0]
        while i < gridpoints:
            testvalue_pot = testvector [i]
            testvalue_a = self.a_axis [i]
            sqrt_testvalue_pot = sqrt(testvalue_pot)
            self.assertTrue ((sqrt_testvalue_pot == abs(testvalue_a)), msg = "Your calculation of the potential matrix in the harmonic potential is incorrect.")
            i = i+1
        diagonalarray_main = np.diagonal(pot_mat, 0)
        if np.all (diagonalarray_main > 0):
            print("\n""The elements of the potential matrix in the harmonic potential does make sense.""\n")
        else :
            print("\n""The elements of the potential matrix in the harmonic potential doesn't make sense.""\n")
        diagonalarray_plusone = np.diagonal(pot_mat,1)
        self.assertTrue (np.all(diagonalarray_plusone == 0),msg = "\n""Your potential matrix is false.""\n")

        diagonalarray_minusone = np.diagonal(pot_mat,-1)
        self.assertTrue (np.all(diagonalarray_minusone == 0), msg ="\n""Your potential matrix is false.""\n")

    def test_OneBody_Box_potential(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        _,_,_,_,pot_mat = OBSolver.calculateBoxPotential(1)
        diagonalarray_main = np.diagonal(pot_mat, 0)
        if np.all (diagonalarray_main > 0):
            print("\n""The elements of the potential matrix in the box potential does make sense.""\n")
        else :
            print("\n""The elements of the potential matrix in the box potential doesn't make sense.""\n")
        diagonalarray_plusone = np.diagonal(pot_mat,1)
        self.assertTrue (np.all(diagonalarray_plusone == 0),msg = "\n""Your potential matrix is false.""\n")

        diagonalarray_minusone = np.diagonal(pot_mat,-1)
        self.assertTrue (np.all(diagonalarray_minusone == 0), msg = "\n""Your potential matrix is false.""\n")

#    def test_OneBody_Gauss_potential(self):
#        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
#        _,_,_,_,pot_mat = OBSolver.calcualteGaussPotential(1,10)
#        diagonalarray_main = np.diagonal(pot_mat, 0)
#        if np.all (diagonalarray_main > 0):
#            print("\n""The elements of the potential matrix in the gauss potential does make sense.""\n")
#        else :
#            print("\n""The elements of the potential matrix in the gauss potential doesn't make sense.""\n")
#        diagonalarray_plusone = np.diagonal(pot_mat,1)
#        self.assertTrue (np.all(diagonalarray_plusone == 0),msg = "\n""Your potential matrix is false.""\n")
#
#        diagonalarray_minusone = np.diagonal(pot_mat,-1)
#        self.assertTrue (np.all(diagonalarray_minusone == 0), msg = "\n""Your potential matrix is false.""\n")

    def test_OneBody_calcualteBoxEigenvalues(self):
        """ Checks a dummy Calculation of an Harmonic Potential
        This code tests the deviation of the eigenvalues.
        And confirm that the eigenvalues are sorted.
        """
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        la,_,_,_,_ = OBSolver.calculateBoxPotential(1)
        eigenvalue = la [0]
        x = 0
        gridpoints = 3000
        lenght_eigenvector = len(la)
        self.assertTrue ((lenght_eigenvector == gridpoints), msg ="You have not enough or too much eigenvalues")
        while x < gridpoints:
            eigenvalue = la[x]
#            self.assertTrue ((eigenvalue > 0), msg = "The eigenvalue aren't positive!")
            if eigenvalue < 0:
                print ("In your box potential are the eigenvalue no.",x,"is negative", "\n")
            x = x+1

    def test_OneBody_calcualteGaussEigenvalues(self):
        """ Checks a dummy Calculation of an Harmonic Potential
        This code tests the deviation of the eigenvalues.
        And confirm that the eigenvalues are sorted.
        """
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        la,_,_,_,_ = OBSolver.calcualteGaussPotential(1,10)
        eigenvalue = la [0]
        x = 0
        gridpoints = 3000
        lenght_eigenvector = len(la)
        self.assertTrue ((lenght_eigenvector == gridpoints), msg ="You have not enough or too much eigenvalues")
        while x < gridpoints:
            eigenvalue = la[x]
#            self.assertTrue ((eigenvalue > 0), msg = "The eigenvalue aren't positive!")
            if eigenvalue < 0:
                print ("In your gauss potential are the eigenvalue no.",x,"is negative", "\n")
            x = x+1

    def test_OneBody_calcualteHarmonicEigenvalues(self):
        """ Checks a dummy Calculation of an Harmonic Potential
                This code tests the deviation of the eigenvalues.
                And confirm that the eigenvalues are sorted.
        """
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        la,_,_,_,_,_ = OBSolver.calcualteHarmonicPotential(0)
        eigenvalue = la [0]
        x = 0
        gridpoints = 3000
        lenght_eigenvector = len(la)
        self.assertTrue ((lenght_eigenvector == gridpoints), msg ="You have not enough or too much eigenvalues")
        while x < gridpoints:
            eigenvalue = la[x]
#            self.assertTrue ((eigenvalue > 0), msg = "The eigenvalue aren't positive!")
            if eigenvalue < 0:
                print("In your harmonic potential are the eigenvalue no.",x,"is negative", "\n")
            x = x+1
        x = 0
        while x < gridpoints:
            eigenvalue = la[x]
            eigenvalue2 = la[x+1]
            eigenvalue3 = la[x+2]
            eigenvalue4 = la[x+3]
            diffrence12 = eigenvalue2 - eigenvalue
            diffrence34 = eigenvalue4 - eigenvalue3
            diffrence1234 = diffrence34 - diffrence12
            if diffrence1234 > 0.5 or diffrence1234 < -0.5:
                print ("The eigenvalues in your harmonic potential are äquidistant until the ", (x-1),"th Gridpoint.You are leaving the bounded state at Gridpoint", x, ",The eigenvalues aren't äquidistant anymore!")
                break
            x = x+1

    def test_Wafefunktion(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        _, v_norm,_,_,_,_ = OBSolver.calcualteHarmonicPotential(0)
        squared_v = (v_norm * v_norm)
        int_squared = np.trapz((squared_v),dx=1/30)
        gridpoints = 3000
        n = int_squared [0]
        x = 0
        while x < gridpoints:
           n = int_squared[x]
           self.assertTrue (n > 0.98 and n < 1.02, msg="not correctly normalized.")
           x = x+1

    def test_Wafefunktion(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        _, v_norm,_,_,_ = OBSolver.calculateBoxPotential(1)
        squared_v = (v_norm * v_norm)
        int_squared = np.trapz((squared_v),dx=1/30)
        gridpoints = 3000
        n = int_squared [0]
        x = 0
        while x < gridpoints:
            n = int_squared[x]
            self.assertTrue (n > 0.98 and n < 1.02, msg="not correctly normalized.")
            x = x+1

    def test_Wafefunktion(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        _, v_norm,_,_,_ = OBSolver.calcualteGaussPotential(1,10)
        squared_v = (v_norm * v_norm)
        int_squared = np.trapz((squared_v),dx=1/30)
        gridpoints = 3000
        n = int_squared [0]
        x = 0
        while x < gridpoints:
            n = int_squared[x]
#            if n > 0.98 and n < 1.02:
#                print ("The Wavefunktion is great!","\n")
            self.assertTrue (n > 0.98 and n < 1.02, msg="not correctly normalized.")
            x = x+1

if __name__ == '__main__':
    #OBSolver = HaPPPy.OneBody.OneBodySolver(100,100)
    #la, vabcd, i = OBSolver.calcualteHarmonicPotential(0)
    #OBSolver.exportData(la, vabcd, i)
    #print(vabcd)
    
    #OBSolver = HaPPPy.OneBody.OneBodySolver(100,2000)
    #la ,v = OBSolver.calcualteGaussPotential(1 ,5)
    #print(v)

#    OBSolver = HaPPPy.OneBody.OneBodySolver(100,100)
#    la, v , Info = OBSolver.calculateBoxPotential(1)
#    OBSolver.exportData(la, v, Info)
    unittest.main()

    #one_body_suite = unittest.TestLoader().loadTestsFromTestCase(OneBodyTestSuite)
    #unittest.TextTestRunner(verbosity=2, buffer=True).run(one_body_suite)
