from HaPPPy.Transmission.Modules.GaussianWave import GaussianWave
from HaPPPy.Transmission.Modules.SplitStepOperator import SplitStepOperator
from HaPPPy.Transmission.Modules.Transmission import Transmission
from HaPPPy.Transmission.Modules.PotentialUtils import PotentialUtils

import numpy
import math

class TransmissionCalculator:
    """
    Calculates the transmission for a particle with energy E
    moving in within a potential V.
    """
    def doCalculation(self):
        return 2

    def calculate_transmission(self, E, V):
        """
        Performs the calculation for a particle with energy E
        moving in within a potential V to return the trasmission rate.

        Internally it creates a gaussian wave package with energy E and uses
        the split step method to iterate the gaussian package through the potential
        using small time steps.

        The transmission is calculated by letting enough time pass so that the wave package
        completely traveles through the potential V and calculating the propabiliy behind the barrier. 
        
        Parameters
        ----------
        E : float
            The energy of the particle in milli eV
        V : Array
            An array of potential values for every point in space x (x is extracted using the indices of V)
            
        Returns
        -------
        rate : float
            The rate of transmission for the particle within the Potential V (number between 0 and 1)
        """
        if not (isinstance(V, Array)):
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
            psi_0 = GaussianWave(width, symmetry_point, E, x, []).create_gauss_x_package()
        except Exception as error:
            raise Exception("An error occured while creating the gaussian package. Error: " + repr(error))

        margin = 10e-3      # Sample value, needs verification
 
        splitStepOperator = SplitStepOperator(x, potential)
        transmission = Transmission()

        psi_1 = splitStepOperator.first_step(psi_0)

        # This step is redundant but serves a technical purpose and was added
        # for clarity reasons.
        psi_n = psi_1
 
        end_point = self.get_index_for_end_of_potential(V)

        # Appliy the split step operator to advance in time until
        # the differences in psi are small enough. It is not certain
        # whether the selected exit method will work. It has to be verified.
        while (True):
            try:
                psi_n1 = splitStepOperator.step(psi_n[:])

                rate_n1 = transmission.calculate(psi_n1, end_point)
                rate_n = transmission.calculate(psi_n, end_point)

                delta = numpy.absolute(rate_n1 - rate_n) # self.calculate_max_diff_pointwise(psi_n, psi_n1)

                psi_n = psi_n1[:]

                # BUG: Delta never changes!
                print(delta)

                if (delta < margin):
                    psi_n = splitStepOperator.final_step(psi_n)
                    break

            except Exception as error:
                raise Exception("An error occured while applying the split step operator. Error: " + repr(error))
            
        try:
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

        if (max(V) == 0):
            return 0

        # Filter V for all entries > 0 and return an array with
        # the indices. The last index is the one we are looking for
        # ([-1] is a sort hand notation for this)
        return [idx for idx, v_i in enumerate(V) if v_i > 0][-1]    
