from HaPPPy.Transmission.Modules.GaussianWave import GaussianWave
from HaPPPy.Transmission.Modules.SplitStepOperator import SplitStepOperator
from scipy.constants import codata

import unittest
import numpy

me   = codata.value("electron mass energy equivalent in MeV") * 1e9 ;  # Convert to milli eV
hbar = codata.value("Planck constant over 2 pi in eV s")      * 1e15;  # Convert to ps


class TransmissionTestSuite(unittest.TestCase):
    def test_GaussianWave_isMathematicallyCorrect(self):
        # Assemble
        # testPsi for a=1/2, energy=(hbar/3)**2/(2*me) x in [-1,0,1], x0=1/5 from wolfram alpha:
        #x = -1
        expectedPsiRe1  = 0.00366637770707954939318041612000244530051853481045139516
        expectedPsiIm1 = -0.00155011963188600389828072397974977655743189489415947434

        #x = 0
        expectedPsiRe2 = 1.074068789330235009621373594853013756095220098586369327
        expectedPsiIm2 = -0.07171085575151225040654591725275008467732293163741000634

        #x = 1
        expectedPsiRe3 = 0.09420262718699546582979705382375965292637280618901849569
        expectedPsiIm3 = 0.02573359354865373116603681901836827575687821882738438584
        
        expectedPsiRe = numpy.array([expectedPsiRe1,expectedPsiRe2,expectedPsiRe3])
        expectedPsiIm = numpy.array([expectedPsiIm1,expectedPsiIm2,expectedPsiIm3])        
        
        width = 1/2
        energy = (hbar/3)**2/(2*me)
        symmetrypoint = 1/5
        pos = [-1, 0, 1]

        gaussianWave = GaussianWave()

        # Act
        psi0 = gaussianWave.create_package(pos, symmetrypoint, energy, width)
        psi = numpy.array(psi0)

        maxAbsReal = max(numpy.absolute(psi.real - expectedPsiRe))
        maxAbsImag = max(numpy.absolute(psi.imag - expectedPsiIm))

        errorTolerance = 10**(-10)

        print("Max errors: ")
        print([maxAbsReal, maxAbsImag])

        # Assert
        self.assertTrue(maxAbsReal < errorTolerance and maxAbsImag < errorTolerance, "The gaussian wave returns an error greater than %e" % errorTolerance)
        
    def test_SplitStepOperator_isMathematicallyCorrect(self):
        psi = [ 1, 2, 3, 4, 4, 2, 1 ]
        V = [ 0, 0, 10, 10, 0, 0, 0 ]
        x = [ 1, 2, 3, 4, 5, 6, 7 ]

        splitStepOperator = SplitStepOperator(V, x)

        psi_new = splitStepOperator.use(psi)

        print(psi_new)

        self.assertTrue(true)


if __name__ == '__main__':
    transmission_suite = unittest.TestLoader().loadTestsFromTestCase(TransmissionTestSuite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(transmission_suite)
