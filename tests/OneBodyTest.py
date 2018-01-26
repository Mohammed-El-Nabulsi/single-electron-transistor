
import HaPPPy
import unittest
import numpy as np
import scipy.integrate
import math
from math import sqrt
from scipy import constants

hbar = constants.hbar
e = constants.e
m = 3
intersection = 0
sigma = 1
n_i = 200
l_i = 50
target_n = 1000
target_l = 200
step_n = 200
step_l = 50

unit_gauss = 1
omega = 1.875537349*(10**13)
class OneBodyTestSuite(unittest.TestCase):
    """A test class to the OneBody module.
        
        First of all, i would like to test the import and if the constants are choosen.
        These are the Main target of my first few tests.
            """
    def test_input(self):
        """ This test should confirm, that the ammount of gridpoints are positive."""
        self.assertTrue((n_i > 0),msg ="You have to choose a positive ammount of gridpoints.")
        
    def test_OneBody_exists(self):
        """ This test, checks wether the One Body module exists.
        """
        self.assertTrue(hasattr(HaPPPy, 'OneBody'),msg ="One Body Module doesn't exist")
        
    def test_ammount_cutted_eigenvalues(self):
        """ This confirm that the user didn't want to get more cutted eigenvalues then eigenvalues exists."""
        self.assertTrue((m <= n_i),msg =" You want to get more cutted eigenvalues, then eigenvalues.")
        
    """ Now the test start to check the kinetic matrices from the chosen potential."""
    def test_OneBody_Harmonic_kinetic(self):
        """ This test check the kinetic Matrix in the harmonic potential by building diagonalarrays."""
        for n in range (n_i,target_n,step_n):
            for l in range (l_i,target_l,step_l):
                OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
                _,_,_,kin_mat,_,_,_,_ = OBSolver.calcualteHarmonicPotential(intersection)
                diagonalarray_harm_kin_main = np.diagonal(kin_mat, 0)
                """ all elements in the diagonalarray of the main "diagonale" should be positive."""
                self.assertTrue(np.all(diagonalarray_harm_kin_main > 0), msg = "\n""The elements of the kinetic matrix in the harmonic potential doesn't make sense.""\n")
                """ all elements of the both "nebendiagonalen" should be -1."""
                diagonalarray_harm_kin_plusone = np.diagonal(kin_mat,1)
                self.assertTrue(np.all(diagonalarray_harm_kin_plusone == -1),msg = "\n""Your kinetic matrix in your harmonic potential is false.""\n")
                diagonalarray_harm_kin_minusone = np.diagonal(kin_mat,-1)
                self.assertTrue(np.all(diagonalarray_harm_kin_minusone == -1), msg ="\n""Your kinetic matrix in your harmonic potential is false.""\n")



    def test_OneBody_Box_kinetic(self):
        for n in range (n_i,target_n,step_n):
            for l in range (l_i,target_l,step_l):
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
        for n in range (n_i,target_n,step_n):
            for l in range (l_i,target_l,step_l):
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

    def test_gridpoints_kinetic_matrix_gauss(self):
        for l in range (l_i,target_l,step_l):
            n = 1
            while n < 200:
                OBSolver = HaPPPy.OneBody.OneBodySolver (l,n,m)
                _,_,_,kin_mat,_,_,_ = OBSolver.calcualteGaussPotential(unit_gauss,sigma)
                diagonalarray = np.diagonal(kin_mat)
                lenght_diagonalarray = len(diagonalarray)
                self.assertTrue((lenght_diagonalarray == n), msg = "The kineticmatrix in your gauss potential has a wrong ammount of values.")
                n = n+1

    def test_gridpoints_kinetic_matrix_harmonic(self):
        for l in range (l_i,target_l,step_l):
            n = 1
            while n < 200:
                OBSolver = HaPPPy.OneBody.OneBodySolver (l,n,m)
                _,_,_,kin_mat,_,_,_,_ = OBSolver.calcualteHarmonicPotential(intersection)
                diagonalarray = np.diagonal(kin_mat)
                lenght_diagonalarray = len(diagonalarray)
                self.assertTrue((lenght_diagonalarray == n), msg = "The kineticmatrix in your harmonic potential has a wrong ammount of values.")
                n = n+1
            
    def test_gridpoints_kinetic_matrix_box(self):
        for l in range (l_i,target_l,step_l):
            n = 1
            while n < 200:
                OBSolver = HaPPPy.OneBody.OneBodySolver (l,n,m)
                _,_,_,kin_mat,_,_,_ = OBSolver.calculateBoxPotential(intersection)
                diagonalarray = np.diagonal(kin_mat)
                lenght_diagonalarray = len(diagonalarray)
                self.assertTrue((lenght_diagonalarray == n), msg = "The kineticmatrix in your box potential has a wrong ammount of values.")
                n = n+1


        """ Now the test start to check the potential matrices from the chosen potential."""
            
    def test_OneBody_Harmonic_potential(self):
        for n in range (n_i,target_n,step_n):
            for l in range (l_i,target_l,step_l):
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
                    self.assertTrue(((abs(sqrt_testvalue_pot) - abs(testvalue_a)) < 0.2), msg = "Your calculation of the potential matrix in the harmonic potential is incorrect.")
                    i = i+1
                diagonalarray_main = np.diagonal(pot_mat, 0)
                if intersection >= 0:
                    self.assertTrue(np.all(diagonalarray_main > 0), msg = "\n""The elements of the potential matrix in the harmonic potential doesn't make sense.""\n")
                    diagonalarray_plusone = np.diagonal(pot_mat,1)
                    self.assertTrue (np.all(diagonalarray_plusone == 0),msg = "\n""Your potential matrix is false.""\n")
                    diagonalarray_minusone = np.diagonal(pot_mat,-1)
                    self.assertTrue (np.all(diagonalarray_minusone == 0), msg ="\n""Your potential matrix is false.""\n")

    def test_OneBody_Box_potential(self):
        for n in range (n_i,target_n,step_n):
            for l in range (l_i,target_l,step_l):
                OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
                _,_,_,_,pot_mat,_,_ = OBSolver.calculateBoxPotential(intersection)
                diagonalarray_main = np.diagonal(pot_mat, 0)
                self.assertTrue(np.all(diagonalarray_main == intersection), msg = "\n""The elements of the potential matrix in the box potential doesn't make sense.""\n")
                diagonalarray_plusone = np.diagonal(pot_mat,1)
                self.assertTrue (np.all(diagonalarray_plusone == 0),msg = "\n""The elements of the potential matrix in the box potential doesn't make sense.""\n")
                diagonalarray_minusone = np.diagonal(pot_mat,-1)
                self.assertTrue (np.all(diagonalarray_minusone == 0), msg = "\n""The elements of the potential matrix in the box potential doesn't make sense.""\n")

    def test_OneBody_Gauss_potential(self):
        for n in range (n_i,target_n,step_n):
            for l in range (l_i,target_l,step_l):
                OBSolver = HaPPPy.OneBody.OneBodySolver(100,n,m)
                _,_,_,_,pot_mat,_,_ = OBSolver.calcualteGaussPotential(unit_gauss,sigma)
                _,_,_,_,_,self.a_axis,_,_ = OBSolver.calcualteHarmonicPotential(intersection)
                lenght_gauss = int(round((n/4),0))
                lenght_gauss2 = int(round((3*lenght_gauss),0))
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
            
        ####### HANDLUNGSBEDARF

    def test_gridpoints_potential_matrix_box(self):
        for l in range (l_i,target_l,step_l):
            n = 1
            while n < 200:
                OBSolver = HaPPPy.OneBody.OneBodySolver (l,n,m)
                _,_,_,_,pot_mat,_,_ = OBSolver.calculateBoxPotential(intersection)
                diagonalarray = np.diagonal(pot_mat)
                lenght_diagonalarray = len(diagonalarray)
                self.assertTrue((lenght_diagonalarray == n), msg = "The potenialmatrix in your box potential has a wrong ammount of values.")
                n = n+1

    def test_gridpoints_potential_matrix_harmonic(self):
        for l in range (l_i,target_l,step_l):
            n = 1
            while n < 200:
                OBSolver = HaPPPy.OneBody.OneBodySolver (l,n,m)
                _,_,_,_,pot_mat,_,_,_ = OBSolver.calcualteHarmonicPotential(intersection)
                diagonalarray = np.diagonal(pot_mat)
                lenght_diagonalarray = len(diagonalarray)
                self.assertTrue((lenght_diagonalarray == n), msg = "The potenialmatrix in your harmonic potential has a wrong ammount of values.")
                n = n+1

    def test_gridpoints_potential_matrix_gauss(self):
        for l in range (l_i,target_l,step_l):
            n = 1
            while n < 200:
                OBSolver = HaPPPy.OneBody.OneBodySolver (l,n,m)
                _,_,_,_,pot_mat,_,_ = OBSolver.calcualteGaussPotential(unit_gauss,sigma)
                diagonalarray = np.diagonal(pot_mat)
                lenght_diagonalarray = len(diagonalarray)
                self.assertTrue((lenght_diagonalarray == n), msg = "The potenialmatrix in your gauss potential has a wrong ammount of values.")
                n = n+1

    """ Now the test start to check the eigenvalues the chosen potential."""

    def test_OneBody_calcualteBoxEigenvalues(self):
        """ Checks a dummy Calculation of an Harmonic Potential
        This code tests the deviation of the eigenvalues.
        And confirm that the eigenvalues are sorted.
            """
        for n in range (n_i,target_n,step_n):
            for l in range (l_i,target_l,step_l):
                OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
                _,_,_,_,_,la_l,_ = OBSolver.calculateBoxPotential(intersection)
                eigenvalue = la_l [0]
                x = 0
                lenght_eigenvector = len(la_l)
                self.assertTrue ((lenght_eigenvector == n), msg ="The ammount of eigenvalues are incorrect!")
                if intersection >= 0:
                    while x < n:
                        eigenvalue = la_l[x]
                        self.assertTrue ((eigenvalue > 0), msg = "The eigenvalue aren't positive!")
                        x = x+1
                    x = 0
                    while x < n-3:
                        first_eigenvalue_box = la_l[x]
                        second_eigenvalue_box = la_l[x+1]
                        third_eigenvalue_box = la_l[x+2]
                        forth_eigenvalue_box = la_l[x+3]
                        diffrence12 = second_eigenvalue_box - first_eigenvalue_box
                        diffrence23 = third_eigenvalue_box - second_eigenvalue_box
                        diffrence1234 = diffrence23-diffrence12
                        if diffrence1234 > 0.5 or diffrence1234 < -0.5:
                            if x > 0:
                                print ("\n","The eigenvalues in your box potential are äquidistant until the ",(x-1),"th Gridpoint.")
                                break
                            else:
                                print ("\n","The eigenvalues in your box potential aren't äquidistant.")
                                break
                        x = x+1

    def test_OneBody_calcualteGaussEigenvalues(self):
        """ Checks a dummy Calculation of an Harmonic Potential
        This code tests the deviation of the eigenvalues.
        And confirm that the eigenvalues are sorted.
        """
        for n in range (n_i,target_n,step_n):
            for l in range (l_i,target_l,step_l):
                OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
                _,_,_,_,_,la_l,_ = OBSolver.calcualteGaussPotential(unit_gauss,sigma)
                eigenvalue = la_l [0]
                x = 0
                lenght_eigenvector = len(la_l)
                self.assertTrue ((lenght_eigenvector == n), msg ="The ammount of eigenvalues are incorrect!")
                while x < n:
                    eigenvalue = la_l[x]
        #          self.assertTrue ((eigenvalue > 0), msg = "The eigenvalue aren't positive!")
        #           if eigenvalue < 0:
        #               print ("In your gauss potential are the eigenvalue no.",x,"is negative", "\n")
                    x = x+1

    def test_OneBody_calcualteHarmonicEigenvalues(self):
        """ Checks a dummy Calculation of an Harmonic Potential
                This code tests the deviation of the eigenvalues.
                And confirm that the eigenvalues are sorted.
        """
        for n in range (n_i,target_n,step_n):
            for l in range (l_i,target_l,step_l):
                OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
                _,_,_,_,_,_,la_l,_ = OBSolver.calcualteHarmonicPotential(intersection)
                eigenvalue = la_l [0]
                x = 0
                lenght_eigenvector = len(la_l)
                self.assertTrue ((lenght_eigenvector == n), msg ="You have not enough or too much eigenvalues")
                while x < n:
                        eigenvalue = la_l[x]
                        self.assertTrue ((eigenvalue > 0), msg = "The eigenvalue aren't positive!")
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
                        if x > 0:
                            print ("\n","The eigenvalues in your harmonic potential are äquidistant until the ",(x-1),"th Gridpoint.")
                            print  ("You hitting the wall of your potential at Gridpoint", x, ",The eigenvalues aren't äquidistant anymore!","\n")
                            break
                        else:
                            print ("\n","The eigenvalues in your harmonic potential aren't äquidistant.")
                            break
                    x = x+1

    def test_output_eigenvalue(self):
        for n in range (n_i,target_n,step_n):
            for l in range (l_i,target_l,step_l):
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
        for n in range (n_i,target_n,step_n):
            for l in range (l_i,target_l,step_l):
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
        for n in range (n_i,target_n,step_n):
            for l in range (l_i,target_l,step_l):
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
        for n in range (n_i,target_n,step_n):
            OBSolver = HaPPPy.OneBody.OneBodySolver(100,2000,m)
            la_gauss,_,_,_,_,_,_ = OBSolver.calcualteGaussPotential(1,1)
            intersection = 0
            la_harmonic,_,_,_,_,_,_,_ = OBSolver.calcualteHarmonicPotential(intersection)
            print(la_harmonic[0])
            intersection = 0
            la_harmonic,_,_,_,_,_,_,_ = OBSolver.calcualteHarmonicPotential(intersection)
            print(la_harmonic[0])
            intersection = 0
            la_harmonic,_,_,_,_,_,_,_ = OBSolver.calcualteHarmonicPotential(intersection)
            print(la_harmonic[0])
            first_eigenvalue_harm = la_harmonic [1] - intersection
            first_eigenvalue_gauss = la_gauss [1]           #ABHÄNGIGKEIT VON L
            self.assertTrue((round(first_eigenvalue_harm,3) == round(first_eigenvalue_gauss,3)), msg ="Your eigenvalues from the gauss potential are wrong calculated")

    def test_output_eigenvalue_harmonic_potential(self):
        for n in range (n_i,target_n,step_n):
            for l in range (l_i,target_l,step_l):
                OBSolver = HaPPPy.OneBody.OneBodySolver (l,3000,m)
                la,_,_,_,_,_,_,_ = OBSolver.calcualteHarmonicPotential(intersection)
                energie_first_eigenvalue = (((1/2)*hbar*omega)/e)*1000
                print("abc", energie_first_eigenvalue)
                print("abc", la[0])
                energie_first_eigenvalue_calc = la [0] - intersection
                p = 0
                if m >= 5:
                    while p < 5:
                        energie_eigenvalue = (((p+(1/2))*hbar*omega)/e)*1000
                        energie_eigenvalue_calc = la [p] - intersection
                        self.assertTrue((round(energie_eigenvalue,1)==round(energie_eigenvalue_calc,1)), msg= "The eigenvalues form the harmonic potential are incorrect!")
                        p = p+1
                else:
                    while p < m:
                        energie_eigenvalue = (((p+(1/2))*hbar*omega)/e)*1000
                        energie_eigenvalue_calc = la [p] - intersection
                        self.assertTrue((round(energie_eigenvalue,1)==round(energie_eigenvalue_calc,1)), msg= "The eigenvalues form the harmonic potential are incorrect!")
                        p = p+1


    """ Now the test start to check the output from the chosen potential."""

    def test_Wafefunktion(self):
        for n in range (n_i,target_n,step_n):
            for l in range (l_i,target_l,step_l):
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
        for n in range (n_i,target_n,step_n):
            for l in range (l_i,target_l,step_l):
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
        for n in range (n_i,target_n,step_n):
            for l in range (l_i,target_l,step_l):
                OBSolver = HaPPPy.OneBody.OneBodySolver(l,n,m)
                _,_,_,_,_,_,v_norm_l = OBSolver.calcualteGaussPotential(unit_gauss,sigma)
                squared_v = (v_norm_l * v_norm_l)
                dx = l/n
                int_squared = np.trapz((squared_v),dx= dx)
                i = int_squared [0]
                x = 0
                while x < n:
                    i = int_squared[x]
                    self.assertTrue (i > 0.85 and i < 1.15, msg="not correctly normalized.")
                    x = x+1



if __name__ == '__main__':
        unittest.main()

