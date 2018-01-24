
import HaPPPy
import unittest
import numpy as np
import scipy.integrate
import math
from math import sqrt
from scipy import constants

hbar = constants.hbar
e = constants.e
l = 100
n = 3000
m = 10
intersection = 1
sigma = 1
dx = l/n
unit_gauss = 1
omega = 1.875537349*(10**13)
lenght_gauss = int(round((n/4),0))
lenght_gauss2 = int(round((3*lenght_gauss),0))

class OneBodyTestSuite(unittest.TestCase):
    """A test class to the OneBody module.

    """

    def test_OneBody_exists(self):
        """ Checks wether the One Body module exists
        """
        self.assertTrue(hasattr(HaPPPy, 'OneBody'),msg ="One Body Module doesn't exist")
    

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
        self.assertTrue(np.all(diagonalarray_main == intersection), msg = "\n""The elements of the potential matrix in the box potential doesn't make sense.""\n")
        diagonalarray_plusone = np.diagonal(pot_mat,1)
        self.assertTrue (np.all(diagonalarray_plusone == 0),msg = "\n""The elements of the potential matrix in the box potential doesn't make sense.""\n")
        diagonalarray_minusone = np.diagonal(pot_mat,-1)
        self.assertTrue (np.all(diagonalarray_minusone == 0), msg = "\n""The elements of the potential matrix in the box potential doesn't make sense.""\n")

    def test_OneBody_Gauss_potential(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
        _,_,_,_,pot_mat,_,_ = OBSolver.calcualteGaussPotential(unit_gauss,sigma)
        _,_,_,_,_,self.a_axis,_,_ = OBSolver.calcualteHarmonicPotential(intersection)
        i = lenght_gauss
        testvector = np.diagonal(pot_mat)
        testvalue_pot = testvector [0]
        testvalue_a = self.a_axis [0]
        while i < lenght_gauss2:
            testvalue_pot = testvector [i]
            testvalue_poti = testvalue_pot
            testvalue_a = self.a_axis [i]
            calc_testvalue_poti = (sigma)*sqrt(-2*np.log(-testvalue_poti/unit_gauss))
            self.assertTrue((round(calc_testvalue_poti,5) == round(abs(testvalue_a),5)), msg = "Your calculation of the potential matrix in the gauss potential is incorrect.")
            i = i+1

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
            self.assertTrue ((eigenvalue > 0), msg = "The eigenvalue aren't positive!")
#            if eigenvalue < 0:
#                print ("In your box potential are the eigenvalue no.",x,"is negative", "\n")
            x = x+1

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
#            if eigenvalue < 0:
#                print ("In your gauss potential are the eigenvalue no.",x,"is negative", "\n")
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

    def test_eigenvalues_gauss_compared(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
        la_gauss,_,_,_,_,_,_ = OBSolver.calcualteGaussPotential(1,1)
        la_harmonic,_,_,_,_,_,_,_ = OBSolver.calcualteHarmonicPotential(intersection)
        first_eigenvalue_harm = la_harmonic [1] - intersection
        first_eigenvalue_gauss = la_gauss [1]
        self.assertTrue((round(first_eigenvalue_harm,3) == round(first_eigenvalue_gauss,3)), msg ="Your eigenvalues from the gauss potential are wrong calculated")

    def test_output_eigenvalue_harmonic_potential(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver (l,n,m)
        la,_,_,_,_,_,_,_ = OBSolver.calcualteHarmonicPotential(intersection)
        energie_first_eigenvalue = (((1/2)*hbar*omega)/e)*1000
        energie_first_eigenvalue_calc = la [0] - intersection
        p = 0
        while p < 5:
            energie_eigenvalue = (((p+(1/2))*hbar*omega)/e)*1000
            energie_eigenvalue_calc = la [p] - intersection
            self.assertTrue((round(energie_eigenvalue,1)==round(energie_eigenvalue_calc,1)), msg= "The eigenvalues form the harmonic potential are incorrect!")
            p = p+1

    def test_gridpoints_potential_matrix_box(self):
        n = 1
        while n < 200:
            OBSolver = HaPPPy.OneBody.OneBodySolver (l,n,m)
            _,_,_,_,pot_mat,_,_ = OBSolver.calculateBoxPotential(intersection)
            diagonalarray = np.diagonal(pot_mat)
            lenght_diagonalarray = len(diagonalarray)
            self.assertTrue((lenght_diagonalarray == n), msg = "The potenialmatrix in your box potential has a wrong ammount of values.")
            n = n+1

    def test_gridpoints_potential_matrix_harmonic(self):
        n = 1
        while n < 200:
            OBSolver = HaPPPy.OneBody.OneBodySolver (l,n,m)
            _,_,_,_,pot_mat,_,_,_ = OBSolver.calcualteHarmonicPotential(intersection)
            diagonalarray = np.diagonal(pot_mat)
            lenght_diagonalarray = len(diagonalarray)
            self.assertTrue((lenght_diagonalarray == n), msg = "The potenialmatrix in your harmonic potential has a wrong ammount of values.")
            n = n+1

    def test_gridpoints_potential_matrix_gauss(self):
        n = 1
        while n < 200:
            OBSolver = HaPPPy.OneBody.OneBodySolver (l,n,m)
            _,_,_,_,pot_mat,_,_ = OBSolver.calcualteGaussPotential(unit_gauss,sigma)
            diagonalarray = np.diagonal(pot_mat)
            lenght_diagonalarray = len(diagonalarray)
            self.assertTrue((lenght_diagonalarray == n), msg = "The potenialmatrix in your gauss potential has a wrong ammount of values.")
            n = n+1

    def test_gridpoints_kinetic_matrix_gauss(self):
        n = 1
        while n < 200:
            OBSolver = HaPPPy.OneBody.OneBodySolver (l,n,m)
            _,_,_,kin_mat,_,_,_ = OBSolver.calcualteGaussPotential(unit_gauss,sigma)
            diagonalarray = np.diagonal(kin_mat)
            lenght_diagonalarray = len(diagonalarray)
            self.assertTrue((lenght_diagonalarray == n), msg = "The kineticmatrix in your gauss potential has a wrong ammount of values.")
            n = n+1

    def test_gridpoints_kinetic_matrix_harmonic(self):
        n = 1
        while n < 200:
            OBSolver = HaPPPy.OneBody.OneBodySolver (l,n,m)
            _,_,_,kin_mat,_,_,_,_ = OBSolver.calcualteHarmonicPotential(intersection)
            diagonalarray = np.diagonal(kin_mat)
            lenght_diagonalarray = len(diagonalarray)
            self.assertTrue((lenght_diagonalarray == n), msg = "The kineticmatrix in your harmonic potential has a wrong ammount of values.")
            n = n+1

    def test_gridpoints_kinetic_matrix_box(self):
        n = 1
        while n < 200:
            OBSolver = HaPPPy.OneBody.OneBodySolver (l,n,m)
            _,_,_,kin_mat,_,_,_ = OBSolver.calculateBoxPotential(intersection)
            diagonalarray = np.diagonal(kin_mat)
            lenght_diagonalarray = len(diagonalarray)
            self.assertTrue((lenght_diagonalarray == n), msg = "The kineticmatrix in your box potential has a wrong ammount of values.")
            n = n+1





if __name__ == '__main__':
    
    unittest.main()

