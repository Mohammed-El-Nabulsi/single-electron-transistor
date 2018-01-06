from HaPPPy.Transmission.Modules.GaussianWave import GaussianWave
from HaPPPy.Transmission.Modules.SplitStepOperator import SplitStepOperator
from HaPPPy.Transmission.Modules.Transmission import Transmission
from HaPPPy.Transmission.Modules.PotentialUtils import PotentialUtils

import numpy
import math

class TransmissionCalculator:
    def doCalculation(self):
        return 2

    def calculate_transmission(self, E, V):
        """
            TODO: Extensive documentation.
        """
        if not (isinstance(V, list)):
            raise ValueError("V and x must be arrays")

        if not (E <= 0 or V):
            raise ValueError("The input parameters are not valid.")
   
        if not (numpy.array(V).ndim == 1):
            raise ValueError("V must be one dimensional arrays")

        potential_utils = PotentialUtils(V, 1)

        potential = potential_utils.potential
        symmetry_point = potential_utils.gauss_symmetry_point
        width = potential_utils.gauss_width

        x = potential_utils.position_grid
        
        psi_0 = None

        try:
            psi_0 = GaussianWave(width, symmetry_point, E, x).create_gauss_x_package()
        except Exception as error:
            raise Exception("An error occured while creating the gaussian package. Error: " + repr(error))

        margin = 10e-3      # Sample value, needs verification
 
        splitStepOperator = SplitStepOperator(x, potential)

        psi_1 = splitStepOperator.first_step(psi_0)

        # This step is redundant but serves a technical purpose and was added
        # for clarity reasons.
        psi_n = psi_1

        # Appliy the split step operator to advance in time until
        # the differences in psi are small enough. It is not certain
        # whether the selected exit method will work. It has to be verified.
        while (True):
            try:
                psi_n1 = splitStepOperator.step(psi_n[:])

                delta = self.calculate_max_diff_pointwise(psi_n, psi_n1)

                psi_n = psi_n1[:]

                # BUG: Delta never changes!
                print(delta)

                if (delta < margin):
                    psi_n = splitStepOperator.final_step(psi_n)
                    break

            except Exception as error:
                raise Exception("An error occured while applying the split step operator. Error: " + repr(error))
            
        try:
            end_point = self.get_index_for_end_of_potential(V)

            return Transmission().calculate(psi_n, end_point)
        except Exception as error:
            raise Exception("An error occured while calculating transmission rates. Error: " + repr(error))

    def calculate_max_diff_pointwise(self, f1, f2):
        return numpy.absolute(max(numpy.absolute(f1, f2)))

    def get_index_for_end_of_potential(self, V):
        """
            Get the last point at which the potential is 0.
            We define this as the end of the potential barrier.
        """

        # Filter V for all entries > 0 and return an array with
        # the indices. The last index is the one we are looking for
        # ([-1] is a sort hand notation for this)
        return [idx for idx, v_i in enumerate(V) if v_i > 0][-1]    

# V = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# E = 1.8

# transmissionCalculator = TransmissionCalculator()

# rate = transmissionCalculator.calculate_transmission(E, V)

# print(rate)
