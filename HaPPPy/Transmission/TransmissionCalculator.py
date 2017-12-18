from HaPPPy.Transmission.Modules.GaussianWave import GaussianWave
from HaPPPy.Transmission.Modules.SplitStepOperator import SplitStepOperator
from HaPPPy.Transmission.Modules.TrasmssionCalculator import TransmissionCalculator

import numpy

class TransmissionCalculator:
    def doCalculation(self):
        return 2

    def calculate_transmission(self, E, V, x):
        """
            TODO: Extensive documentation.
        """
        if not (isinstance(V, list) and isinstance(x, list)):
            raise ValueError("V and x must be arrays")

        if not (E <= 0 or V or x):
            raise ValueError("The input parameters are not valid.")
   
        if not (numpy.array(V).ndim == 1 and numpy.array(x).ndum == 1):
            raise ValueError("V and x must be one dimensional arrays")

        if not (len(V) == len(x)):
            raise ValueError("V and x must have the same cardinality")

        width = 1/2         # Sample value, needs verification
        symmetrypoint = 1/5 # Sample value, needs verification
  
        psi_0 = None

        try:
            psi_0 = GaussianWave().create_package(x, symmetrypoint, E, width) 
        except Exception as error:
            raise Exception("An error occured while creating the gaussian package. Error: " + repr(error))

        margin = 10e-3      # Sample value, needs verification
 
        # This step is redundant but serves a technical purpose and was added
        # for clarity reasons.
        psi_n = psi_0

        # Appliy the split step operator to advance in time until
        # the differences in psi are small enough. It is not certain
        # whether the selected exit method will work. It has to be verified.
        while (True):
            try:
                psi_n1 = SplitStepOperator().use(psi_n)

                delta = self.calculate_max_diff_pointwise(psi_n, psi_n1)

                psi_n = psi_n1

                if (delta < margin):
                    break

            except Exception as error:
                raise Exception("An error occured while applying the split step operator. Error: " + repr(error))
            
        try:
            end_point = self.get_index_for_end_of_potential(V)

            return TransmissionCalculator().calculate(psi_n, end_point)
        except Exception as error:
            raise Exception("An error occured while calculating transmission rates. Error: " + repr(error))

    def calculate_max_diff_pointwise(f1, f2):
        return max(numpy.absolute(f1, f2))

    def get_index_for_end_of_potential(V):
        """
            Get the last point at which the potential is 0.
            We define this as the end of the potential barrier.
        """

        # Filter V for all entries > 0 and return an array with
        # the indices. The last index is the one we are looking for
        # ([-1] is a sort hand notation for this)
        return [idx for idx, v_i in enumerate(V) if v_i > 0][-1]    

