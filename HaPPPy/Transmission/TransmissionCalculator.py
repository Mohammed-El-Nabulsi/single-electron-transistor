

class TransmissionCalculator:
    """ Calculates the transmission coefficients

    TODO: Add more documentation

    """


    def __init__(self):
        """ The constructor.
        """

        print("Hello from the TransmissionCalculator, what can I do for you?")

    def doCalculation(self, E, V, x):
        """ Some dummy calculation
        
        Returns:
            double. The result

        """
        
        a = 1;

        psi_0 = GaussianWave().CreatePackage(a, E, x) 

        delta = infinity;
        margin = 10e-3
 
        while (delta > margin):
            if (Init):
                psi_n = psi0

            psi_n1 = SplitStepOperator().Use(psi_n)

            delta = GetLargestDiff(psi_n, psi_n1)

            psi_n = psi_n1

        return CalculateTransmission().Calculate(psi_n, V)

        

        return 1.0+1.0
