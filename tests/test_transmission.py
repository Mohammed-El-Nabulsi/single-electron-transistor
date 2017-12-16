from HaPPPy.Transmission.Modules.GaussianWave import GaussianWave

import unittest

class TransmissionTestSuite(unittest.TestCase):
    def test_GaussianWave_isMathematicallyCorrect(self):
        # Assemble
        # testPsi for energie = x = a = 1 from wolfram alpha:
        expectedPsiReal = 0.328573471238826838814724884575412086439812622961237236
        expectedPsiImag = -0.0046267725731805804274342494105112176810998628709738939
        
        energy = 1
        x = 1
        a = 1

        gaussianWave = GaussianWave()

        # Act
        psi = gaussianWave.CreatePackage(x, energy, a)

        absReal = abs(psi.real - expectedPsiReal)
        absImag = abs(psi.imag - expectedPsiImag)

        errorTolerance = 10**(-10)

        # Assert
        self.assertTrue(absReal < errorTolerance and absImag < errorTolerance)
        
if __name__ == '__main__':
    transmission_suite = unittest.TestLoader().loadTestsFromTestCase(TransmissionTestSuite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(transmission_suite)