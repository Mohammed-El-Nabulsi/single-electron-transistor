
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
        print("One Body module exists" )

    def test_OneBody_Harmonic_kinetic(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        kin_mat, self.n, i = OBSolver.calcualteHarmonicPotential(0)
        """ check the kin Matrix of the harmonic potential
            Confirm that all matrix elements are positive"""
        lowest_value_harmonic = np.min(kin_mat)
        print ("The Element with the lowest value in your kinetic matrix in the box potential are:", lowest_value_harmonic)
        if lowest_value_harmonic > 0:
            print("The kinetic Matrix in your harmonic potential is correct.")
        else:
            print ("Your make a mistake in your calculation of your kinetic matrix in your harmonic potential")

    def test_OneBody_Box_kinetic(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        kin_mat, self.n, i = OBSolver.calculateBoxPotential(0)
        """ check the kin Matrix of the box potential
        Confirm that all matrix elements are positive"""
        lowest_value_box = np.min(kin_mat)
        print ("The Element with the lowest value in your kinetic matrix in the box potential are:", lowest_value_box)
        if lowest_value_box > 0:
            print("The kinetic Matrix in your box potential is correct.")
        else:
            print ("Your make a mistake in your calculation of your kinetic matrix in your box potential")


    def test_OneBody_Gauss_kinetic(self):
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        kin_mat, self.n, i = OBSolver.calcualteGaussPotential(0,1)
        """ check the kin Matrix of the gauss potential
            Confirm that all matrix elements are positive"""
        lowest_value_gauss = np.min(kin_mat)
        print ("The Element with the lowest value in your kinetic matrix in the gauss potential are:", lowest_value_gauss)
        if lowest_value_gauss > 0:
            print("The kinetic Matrix in your gauss potential is correct.")
        else:
            print ("Your make a mistake in your calculation of your kinetic matrix in your gauss potential")


    def test_OneBody_calcualteHarmonicPotential(self):
        """ Checks a dummy Calculation of an Harmonic Potential
                This code tests the deviation of the eigenvalues.
        """
        OBSolver = HaPPPy.OneBody.OneBodySolver(100, 3000)
        la, vabcd, i = OBSolver.calcualteHarmonicPotential(0)

        zeroeigenvalue = la [0]
        firsteigenvalue = la[1]
        secondeigenvalue = la[2]
        thirdeigenvalue = la[3]
        fourtheigenvalue = la[4]
    
        firstcalculation = (firsteigenvalue-zeroeigenvalue)
        secondcalculation = (secondeigenvalue-firsteigenvalue)
        thirdcalculation = (thirdeigenvalue-secondeigenvalue)
        fourthcalculation = (fourtheigenvalue-thirdeigenvalue)
        print("The Diffrence between the first eigenvalues are:",firstcalculation ,",", secondcalculation ,",", thirdcalculation ,",", fourthcalculation)
        finalcalculation_eigenvalues =(fourthcalculation-thirdcalculation)
        if finalcalculation_eigenvalues > 1:
            print("The Eigenvalues are wrong, your calculation is bullshit")
        else :
            print("The Eigenvalues are correct, well done!")
        print ("The deviation of the Eigenvalues are:", finalcalculation_eigenvalues)


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
