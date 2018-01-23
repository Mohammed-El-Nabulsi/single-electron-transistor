
import HaPPPy
import unittest
import numpy as np
import scipy.integrate
from math import sqrt
from scipy import constants

hbar = constants.hbar
e = constants.e
l = 100
n = 3000
m = 10
intersection = 0
sigma = 10
dx = l/n
unit_gauss = 1

class OneBodyTestSuite(unittest.TestCase):
    """A test class to the OneBody module.

    """

    def test_OneBody_exists(self):
        """ Checks wether the One Body module exists
        """
        self.assertTrue(hasattr(HaPPPy, 'OneBody'))
        print("\n""One Body module exists""\n" )

    def test_OneBody_Harmonic_kinetic(self):
        """
            abc
            """
        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
        _,_,_,kin_mat,_,_,_,_ = OBSolver.calcualteHarmonicPotential(intersection)
        """ check the kin Matrix of the harmonic potential
            Confirm that all matrix elements are positive"""
        diagonalarray_harm_kin_main = np.diagonal(kin_mat, 0)
        self.assertTrue(np.all(diagonalarray_harm_kin_main > 0), msg = "\n""The elements of the kinetic matrix in the harmonic potential doesn't make sense.""\n")
        diagonalarray_harm_kin_plusone = np.diagonal(kin_mat,1)
        self.assertTrue(np.all(diagonalarray_harm_kin_plusone == -1),msg = "\n""Your kinetic matrix in your harmonic potential is false.""\n")

        diagonalarray_harm_kin_minusone = np.diagonal(kin_mat,-1)
        self.assertTrue(np.all(diagonalarray_harm_kin_minusone == -1), msg ="\n""Your kinetic matrix in your harmonic potential is false.""\n")

    def test_OneBody_Box_kinetic(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
        _,_,_,kin_mat,_,_,_ = OBSolver.calculateBoxPotential(intersection)
        """ check the kin Matrix of the box potential
        Confirm that all matrix elements are positive"""
        diagonalarray_box_kin_main = np.diagonal(kin_mat, 0)
        self.assertTrue(np.all(diagonalarray_box_kin_main > 0), msg = "\n""The elements of the kinetic matrix in the box potential doesn't make sense.""\n")
        diagonalarray_box_kin_plusone = np.diagonal(kin_mat,1)
        self.assertTrue (np.all(diagonalarray_box_kin_plusone == -1),msg = "\n""Your kinetic matrix in your box potential is false.""\n")

        diagonalarray_box_kin_minusone = np.diagonal(kin_mat,-1)
        self.assertTrue (np.all(diagonalarray_box_kin_minusone == -1), msg ="\n""Your kinetic matrix in your box potential is false.""\n")

    def test_OneBody_Gauss_kinetic(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
        _,_,_,kin_mat,_,_,_ = OBSolver.calcualteGaussPotential(unit_gauss,sigma)
        """ check the kin Matrix of the gauss potential
            Confirm that all matrix elements are positive"""
        diagonalarray_gauss_kin_main = np.diagonal(kin_mat, 0)
        self.assertTrue(np.all(diagonalarray_gauss_kin_main > 0), msg = "\n""The elements of the kinetic matrix in the gauss potential doesn't make sense.""\n")
        diagonalarray_gauss_kin_plusone = np.diagonal(kin_mat,1)
        self.assertTrue (np.all(diagonalarray_gauss_kin_plusone == -1),msg = "\n""Your kinetic matrix in your gauss potential is false.""\n")

        diagonalarray_gauss_kin_minusone = np.diagonal(kin_mat,-1)
        self.assertTrue (np.all(diagonalarray_gauss_kin_minusone == -1), msg ="\n""Your kinetic matrix in your gauss potential is false.""\n")

    def test_OneBody_Harmonic_potential(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
        _,_,_,_,pot_mat,self.a_axis,_,_ = OBSolver.calcualteHarmonicPotential(intersection)
        i = 0
        testvector = np.diagonal(pot_mat)
        testvalue_pot = testvector [0]
        testvalue_a = self.a_axis [0]
        while i < n:
            testvalue_pot = testvector [i]
            testvalue_poti = testvalue_pot - intersection
            testvalue_a = self.a_axis [i]
            sqrt_testvalue_pot = sqrt(testvalue_poti)
            self.assertTrue((round((sqrt_testvalue_pot),5) == round(abs(testvalue_a),5)), msg = "Your calculation of the potential matrix in the harmonic potential is incorrect.")
            i = i+1
        diagonalarray_main = np.diagonal(pot_mat, 0)
        self.assertTrue(np.all(diagonalarray_main > 0), msg = "\n""The elements of the potential matrix in the harmonic potential doesn't make sense.""\n")
        diagonalarray_plusone = np.diagonal(pot_mat,1)
        self.assertTrue (np.all(diagonalarray_plusone == 0),msg = "\n""Your potential matrix is false.""\n")

        diagonalarray_minusone = np.diagonal(pot_mat,-1)
        self.assertTrue (np.all(diagonalarray_minusone == 0), msg ="\n""Your potential matrix is false.""\n")

    def test_OneBody_Box_potential(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
        _,_,_,_,pot_mat,_,_ = OBSolver.calculateBoxPotential(intersection)
        diagonalarray_main = np.diagonal(pot_mat, 0)
        self.assertTrue(np.all(diagonalarray_main >= 0), msg = "\n""The elements of the potential matrix in the box potential doesn't make sense.""\n")
        diagonalarray_plusone = np.diagonal(pot_mat,1)
        self.assertTrue (np.all(diagonalarray_plusone == 0),msg = "\n""The elements of the potential matrix in the box potential doesn't make sense.""\n")
        diagonalarray_minusone = np.diagonal(pot_mat,-1)
        self.assertTrue (np.all(diagonalarray_minusone == 0), msg = "\n""The elements of the potential matrix in the box potential doesn't make sense.""\n")

#    def test_OneBody_Gauss_potential(self):
#        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
#        _,_,_,_,pot_mat,_,_ = OBSolver.calcualteGaussPotential(unit_gauss,sigma)
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
        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
        _,_,_,_,_,la_l,_ = OBSolver.calculateBoxPotential(intersection)
        eigenvalue = la_l [0]
        x = 0
        lenght_eigenvector = len(la_l)
        self.assertTrue ((lenght_eigenvector == n), msg ="The ammount of eigenvalues are incorrect!")
        while x < n:
            eigenvalue = la_l[x]
#            self.assertTrue ((eigenvalue > 0), msg = "The eigenvalue aren't positive!")
            if eigenvalue < 0:
                print ("In your box potential are the eigenvalue no.",x,"is negative", "\n")
            x = x+1

#        for i in range(2000):
#            bla = HaPPPy.OneBody.OneBodySolver(100, i, 10):
#            bla.10
#
#            assertTrue(len(reuslt) == i)

    def test_OneBody_calcualteGaussEigenvalues(self):
        """ Checks a dummy Calculation of an Harmonic Potential
        This code tests the deviation of the eigenvalues.
        And confirm that the eigenvalues are sorted.
        """
        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
        _,_,_,_,_,la_l,_ = OBSolver.calcualteGaussPotential(unit_gauss,sigma)
        eigenvalue = la_l [0]
        x = 0
        lenght_eigenvector = len(la_l)
        self.assertTrue ((lenght_eigenvector == n), msg ="The ammount of eigenvalues are incorrect!")
        while x < n:
            eigenvalue = la_l[x]
#            self.assertTrue ((eigenvalue > 0), msg = "The eigenvalue aren't positive!")
            if eigenvalue < 0:
                print ("In your gauss potential are the eigenvalue no.",x,"is negative", "\n")
            x = x+1

    def test_OneBody_calcualteHarmonicEigenvalues(self):
        """ Checks a dummy Calculation of an Harmonic Potential
                This code tests the deviation of the eigenvalues.
                And confirm that the eigenvalues are sorted.
        """
        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
        _,_,_,_,_,_,la_l,_ = OBSolver.calcualteHarmonicPotential(intersection)
        eigenvalue = la_l [0]
        x = 0
        lenght_eigenvector = len(la_l)
        self.assertTrue ((lenght_eigenvector == n), msg ="You have not enough or too much eigenvalues")
        while x < n:
            eigenvalue = la_l[x]
#            self.assertTrue ((eigenvalue > 0), msg = "The eigenvalue aren't positive!")
            if eigenvalue < 0:
                print("In your harmonic potential are the eigenvalue no.",x,"is negative", "\n")
            x = x+1
        x = 0
        while x < n:
            eigenvalue = la_l[x]
            eigenvalue2 = la_l[x+1]
            eigenvalue3 = la_l[x+2]
            eigenvalue4 = la_l[x+3]
            diffrence12 = eigenvalue2 - eigenvalue
            diffrence34 = eigenvalue4 - eigenvalue3
            diffrence1234 = diffrence34 - diffrence12
            if diffrence1234 > 0.05 or diffrence1234 < -0.05:
                print ("\n","The eigenvalues in your harmonic potential are äquidistant until the ", (x-1),"th Gridpoint.")
                print  ("You hitting the wall of your potential at Gridpoint", x, ",The eigenvalues aren't äquidistant anymore!","\n")
                break
            x = x+1

    def test_Wafefunktion(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
        _,_,_,_,_,_,_,v_norm_l = OBSolver.calcualteHarmonicPotential(intersection)
        squared_v = (v_norm_l * v_norm_l)
        int_squared = np.trapz((squared_v),dx)
        i = int_squared [0]
        x = 0
        while x < n:
           i = int_squared[x]
           self.assertTrue (i > 0.98 and i < 1.02, msg="not correctly normalized.")
           x = x+1

    def test_Wafefunktion(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
        _,_,_,_,_,_,v_norm_l = OBSolver.calculateBoxPotential(intersection)
        squared_v = (v_norm_l * v_norm_l)
        int_squared = np.trapz((squared_v),dx)
        u = int_squared [0]
        x = 0
        while x < n:
            u = int_squared[x]
            self.assertTrue (u > 0.98 and u < 1.02, msg="not correctly normalized.")
            x = x+1

    def test_Wafefunktion(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
        _,_,_,_,_,_,v_norm_l = OBSolver.calcualteGaussPotential(unit_gauss,sigma)
        squared_v = (v_norm_l * v_norm_l)
        int_squared = np.trapz((squared_v),dx= dx)
        i = int_squared [0]
        x = 0
        while x < n:
            i = int_squared[x]
#            if n > 0.98 and n < 1.02:
#                print ("The Wavefunktion is great!","\n")
            self.assertTrue (i > 0.98 and i < 1.02, msg="not correctly normalized.")
            x = x+1

    def test_output_eigenvalue(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
        la, v_norm,_,_,_,la_l, v_norm_l = OBSolver.calculateBoxPotential(intersection)
        cut_la_l = la_l [0]
        x = 0
        while x < m:
            cut_la_l = la_l[x]
            check_la = la[x]
            self.assertTrue((cut_la_l == check_la), msg ="Your choosen eigenvalues do not correlate with the original")
            x = x+1
    def test_output_eigenvalue(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
        la, v_norm,_,_,_,_,la_l, v_norm_l = OBSolver.calcualteHarmonicPotential(intersection)
        cut_la_l = la_l [0]
        x = 0
        while x < m:
            cut_la_l = la_l[x]
            check_la = la[x]
            self.assertTrue((cut_la_l == check_la), msg ="Your choosen eigenvalues do not correlate with the original")
            x = x+1
    def test_output_eigenvalue(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
        la, v_norm,_,_,_,la_l, v_norm_l = OBSolver.calcualteGaussPotential(unit_gauss,sigma)
        cut_la_l = la_l [0]
        x = 0
        while x < m:
            cut_la_l = la_l[x]
            check_la = la[x]
            self.assertTrue((cut_la_l == check_la), msg ="Your choosen eigenvalues do not correlate with the original")
            x = x+1

#    def test_eigenvalues_gauss_compared(self):
#        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)


    def test_output_eigenvalue_harmonic_potential(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver (l,n,m)
        la,_,_,_,_,_,_,_ = OBSolver.calcualteHarmonicPotential(intersection)
        f = 1.876077526e13
        energie_first_eigenvalue = (((1/2)*hbar*f)/e)*1000
        energie_first_eigenvalue_calc = la [0]
        self.assertTrue((round(energie_first_eigenvalue,5) == round(energie_first_eigenvalue_calc,5)), msg ="The calculated eigenvalues aren't equal too the numeric calculated eigenvalues")

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
