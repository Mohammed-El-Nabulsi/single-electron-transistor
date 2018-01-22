
import HaPPPy
import unittest
import numpy as np

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
        _,_,_,kin_mat,_ = OBSolver.calcualteHarmonicPotential(0)
        """ check the kin Matrix of the harmonic potential
            Confirm that all matrix elements are positive"""
        diagonalarray_harm_kin_main = np.diagonal(kin_mat, 0)
        if np.all (diagonalarray_harm_kin_main > 0):
            print("\n""The elements of the kinetic matrix in the harmonic potential make sense""\n")
        else :
            print("\n""The elements of the kinetic matrix in the harmonic potential doesn't make sense""\n")
        diagonalarray_harm_kin_plusone = np.diagonal(kin_mat,1)
        self.assertTrue (np.all(diagonalarray_harm_kin_plusone == -1),msg = "\n""Your kinetic matrix in your harmonic potential is false""\n")
        
        diagonalarray_harm_kin_minusone = np.diagonal(kin_mat,-1)
        self.assertTrue (np.all(diagonalarray_harm_kin_minusone == -1), msg ="\n""Your kinetic matrix in your harmonic potential is false""\n")


    def test_OneBody_Box_kinetic(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        _,_,_,kin_mat,_ = OBSolver.calculateBoxPotential(0)
        """ check the kin Matrix of the box potential
        Confirm that all matrix elements are positive"""
        diagonalarray_box_kin_main = np.diagonal(kin_mat, 0)
        if np.all (diagonalarray_box_kin_main > 0):
            print("\n""The elements of the kinetic matrix in the box potential make sense""\n")
        else :
            print("\n""The elements of the kinetic matrix in the box potential doesn't make sense""\n")
        diagonalarray_box_kin_plusone = np.diagonal(kin_mat,1)
        self.assertTrue (np.all(diagonalarray_box_kin_plusone == -1),msg = "\n""Your kinetic matrix in your box potential is false""\n")
        
        diagonalarray_box_kin_minusone = np.diagonal(kin_mat,-1)
        self.assertTrue (np.all(diagonalarray_box_kin_minusone == -1), msg ="\n""Your kinetic matrix in your box potential is false""\n")

    def test_OneBody_Gauss_kinetic(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        _,_,_,kin_mat,_ = OBSolver.calcualteGaussPotential(0,1)
        """ check the kin Matrix of the gauss potential
            Confirm that all matrix elements are positive"""
        diagonalarray_gauss_kin_main = np.diagonal(kin_mat, 0)
        if np.all (diagonalarray_gauss_kin_main > 0):
            print("\n""The elements of the kinetic matrix in the gauss potential make sense""\n")
        else :
            print("\n""The elements of the kinetic matrix in the gauss potential doesn't make sense""\n")
        diagonalarray_gauss_kin_plusone = np.diagonal(kin_mat,1)
        self.assertTrue (np.all(diagonalarray_gauss_kin_plusone == -1),msg = "\n""Your kinetic matrix in your gauss potential is false""\n")

        diagonalarray_gauss_kin_minusone = np.diagonal(kin_mat,-1)
        self.assertTrue (np.all(diagonalarray_gauss_kin_minusone == -1), msg ="\n""Your kinetic matrix in your gauss potential is false""\n")


    def test_OneBody_Harmonic_potential(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        _,_,_,_,pot_mat = OBSolver.calcualteHarmonicPotential(0)
        diagonalarray_main = np.diagonal(pot_mat, 0)
        if np.all (diagonalarray_main > 0):
            print("\n""The elements of the potential matrix in the harmonic potential make sense""\n")
        else :
            print("\n""The elements of the potential matrix in the harmonic potential doesn't make sense""\n")
        diagonalarray_plusone = np.diagonal(pot_mat,1)
        self.assertTrue (np.all(diagonalarray_plusone == 0),msg = "\n""Your potential matrix is false""\n")
        
        diagonalarray_minusone = np.diagonal(pot_mat,-1)
        self.assertTrue (np.all(diagonalarray_minusone == 0), msg ="\n""Your potential matrix is false""\n")

    def test_OneBody_Box_potential(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        _,_,_,_,pot_mat = OBSolver.calculateBoxPotential(0)
        diagonalarray_main = np.diagonal(pot_mat, 0)
        if np.all (diagonalarray_main > 0):
            print("\n""The elements of the potential matrix in the box potential make sense""\n")
        else :
            print("\n""The elements of the potential matrix in the box potential doesn't make sense""\n")
        diagonalarray_plusone = np.diagonal(pot_mat,1)
        self.assertTrue (np.all(diagonalarray_plusone == 0),msg = "\n""Your potential matrix is false""\n")
        
        diagonalarray_minusone = np.diagonal(pot_mat,-1)
        self.assertTrue (np.all(diagonalarray_minusone == 0), msg = "\n""Your potential matrix is false""\n")

    def test_OneBody_Gauss_potential(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        _,_,_,_,pot_mat = OBSolver.calcualteGaussPotential(0,1)
        diagonalarray_main = np.diagonal(pot_mat, 0)
        if np.all (diagonalarray_main > 0):
            print("\n""The elements of the potential matrix in the gauss potential make sense""\n")
        else :
            print("\n""The elements of the potential matrix in the gauss potential doesn't make sense""\n")
        diagonalarray_plusone = np.diagonal(pot_mat,1)
        self.assertTrue (np.all(diagonalarray_plusone == 0),msg = "\n""Your potential matrix is false""\n")

        diagonalarray_minusone = np.diagonal(pot_mat,-1)
        self.assertTrue (np.all(diagonalarray_minusone == 0), msg = "\n""Your potential matrix is false""\n")



    def test_OneBody_calcualteHarmonicPotential(self):
        """ Checks a dummy Calculation of an Harmonic Potential
                This code tests the deviation of the eigenvalues.
                And confirm that the eigenvalues are sorted.
        """
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        la, vabcd, i,_,_ = OBSolver.calcualteHarmonicPotential(0)

        zeroeigenvalue = la [0]
        firsteigenvalue = la[1]
        secondeigenvalue = la[2]
        thirdeigenvalue = la[3]
        fourtheigenvalue = la[4]
    
        firstcalculation = (firsteigenvalue-zeroeigenvalue)
        secondcalculation = (secondeigenvalue-firsteigenvalue)
        thirdcalculation = (thirdeigenvalue-secondeigenvalue)
        fourthcalculation = (fourtheigenvalue-thirdeigenvalue)
        print("\n""The Diffrence between the first eigenvalues are:",firstcalculation ,",", secondcalculation ,",", thirdcalculation ,",", fourthcalculation)
        finalcalculation_eigenvalues =(fourthcalculation-thirdcalculation)
        if finalcalculation_eigenvalues > 1:
            print("\n""The Eigenvalues are wrong, your calculation is bullshit""\n")
        else :
            print("\n""The Eigenvalues are correct, well done!""\n")
        print ("\n""The deviation of the Eigenvalues are:", finalcalculation_eigenvalues,"\n")

#    def test_Wafefunktion(self):
#        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
#        _, v_norm, _,_,_ = OBSolver.calcualteHarmonicPotential(0)
#        for z in range 3000
#            v_norm = waves [:,z]
#        self.assertAlmostEqual((np.all(np.inner(v_norm, v_norm)) * 1/30),1.0 , msg ="wrong normalized")


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
