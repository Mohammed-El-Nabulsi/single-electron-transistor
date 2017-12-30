from HaPPPy.Transmission.Modules.GaussianWave import GaussianWave
from HaPPPy.Transmission.Modules.SplitStepOperator import SplitStepOperator
from HaPPPy.Transmission.Modules.TransmissionCalculator import TransmissionCalculator
from HaPPPy.Transmission.Modules.Fourier import Fourier
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
        
    def test_Transmissioncalculator_isMathematicallyCorrect(self):
        psi = [ -2, 1, 2, 1j ]
        start = 2
        expectedx = 0,5
        x = TaransmissionCalculator.calculate(psi,start)

        self.assertTrue(x == expectedx)


if __name__ == '__main__':
    transmission_suite = unittest.TestLoader().loadTestsFromTestCase(TransmissionTestSuite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(transmission_suite)
    
  def test_Fourier_isMathematicallyCorrect(self):
# Assemble
# test x = np.arange(-1,-1+5/2,1/2), function = np.array([np.sin(0.123*x)+np.sin(0.234*x) for x in positions])
        
    positions = numpy.arange(-1,-1+5/2,1/2)

    function = numpy.array([numpy.sin(0.123*x)+numpy.sin(0.234*x) for x in positions])
        
    #x[0]
    expected_idft_real1 = 0
    expected_idft_imag1 = 0
    #x[1]
    expected_idft_real2 = 0
    expected_idft_imag2 = -0.44666733802575149239
    #x[2]
    expected_idft_real3 = 0
    expected_idft_imag3 = 0
    #x[3]
    expected_idft_real4 = 0
    expected_idft_imag4 = 0.44666733802575149239
    #x[4]
    expected_idft_real5 = 0
    expected_idft_imag5 = 0
    
    #create expected array
    expected_idft = numpy.array(
            [expected_idft_real1+1*1j*expected_idft_imag1,
             expected_idft_real2+1*1j*expected_idft_imag2,
             expected_idft_real3+1*1j*expected_idft_imag3,
             expected_idft_real4+1*1j*expected_idft_imag4,
             expected_idft_real5+1*1j*expected_idft_imag5
             ]
            )
    
    expected_idft_real = expected_idft.real
    expected_idft_imag = expected_idft.imag
    
    
    # Act
    fourier = Fourier(positions)
    
    maxAbsReal = max(numpy.absolute(fourier.idft(function).real - expected_idft_real))
    maxAbsImag = max(numpy.absolute(fourier.idft(function).imag - expected_idft_imag))
    
    errorTolerance = 10**(-15)
    
    print("Max errors: ")
    print([maxAbsReal, maxAbsImag])
    
    # Assert
    self.assertTrue(maxAbsReal < errorTolerance and maxAbsImag < errorTolerance, "The Fourier-Method returns an error greater than %e" % errorTolerance)
