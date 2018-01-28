from HaPPPy.Transmission.TransmissionCalculator import TransmissionCalculator
from scipy.constants import codata

import unittest
import pytest
import numpy as np
import matplotlib.pyplot as plt

class TransmissionTestSuite(unittest.TestCase):
    def test_validation_does_not_allow_negative_electron_energies(self):
        # Assemble
        E = -1
        dx = 0.1
        barrier = np.array(20+np.zeros(3000))

        transmission_calculator = TransmissionCalculator()
        
        # Act
        with pytest.raises(ValueError) as exception_results:
            transmission_calculator.calculate_transmission(E, barrier, dx)

        # Assert
        self.assertTrue("Electron energy must be greater than 0." in str(exception_results.value))

    def test_validation_does_not_allow_electron_energies_bigger_than_potential(self):
        # Assemble
        E = 25 
        dx = 0.1
        barrier = np.array(20+np.zeros(3000))

        transmission_calculator = TransmissionCalculator()
        
        # Act
        with pytest.raises(ValueError) as exception_results:
            transmission_calculator.calculate_transmission(E, barrier, dx)

        # Assert
        self.assertTrue("Electron energy cannot be bigger than max potential value." in str(exception_results.value))
        
    def test_validation_does_not_allow_invalid_potential(self):
        # Assemble
        E = 25 
        dx = 0.1
        barrier = [np.array(20+np.zeros(3000)), np.array(20+np.zeros(3000))]
  
        transmission_calculator = TransmissionCalculator()
        
        # Act
        with pytest.raises(ValueError) as exception_results1:
            transmission_calculator.calculate_transmission(E, barrier, dx)

        with pytest.raises(ValueError) as exception_results2:
            transmission_calculator.calculate_transmission(E, "some string", dx)

        # Assert
        self.assertTrue("The potential must be an array of one dimension" in str(exception_results1.value))
        self.assertTrue("The potential must be an array of one dimension" in str(exception_results2.value))

        
    def test_validation_does_not_allow_invalid_dx(self):
        # Assemble
        E = 10 
        barrier = np.array(20+np.zeros(3000))
  
        transmission_calculator = TransmissionCalculator()
        
        # Act
        with pytest.raises(ValueError) as exception_results1:
            transmission_calculator.calculate_transmission(E, barrier, 0)

        with pytest.raises(ValueError) as exception_results2:
            transmission_calculator.calculate_transmission(E, barrier, -1)

        # Assert
        self.assertTrue("dx must be greater than 0" in str(exception_results1.value))
        self.assertTrue("dx must be greater than 0" in str(exception_results2.value))

    def test_width_of_free_gaussian_package_grows_correctly(self):
        hbar = 1
        me = 1

        E = 0.5
        dx = 0.1
        barrier = np.zeros(300)

        initial_package_width = 1
        error_tolerance = 1 
  
        package_widths = []
        package_expected_widths = []
 
        # ref: https://stackoverflow.com/a/16489955
        def find_width_at_half_max_height(X,Y):
            half_max = max(Y) / 2.
            #find when function crosses line half_max (when sign of diff flips)
            #take the 'derivative' of signum(half_max - Y[])
            d = np.sign(half_max - np.array(Y[0:-1])) - np.sign(half_max - np.array(Y[1:]))
            #plot(X,d) #if you are interested
            #find the left and right most indexes
            left_idx = np.argwhere(d > 0)[0]
            right_idx = np.argwhere(d < 0)[-1]

            return (X[right_idx] - X[left_idx])[0] #return the difference (full width)

        def _step_callback(self, psi, psi_squared, x, t, b, z):
            package_widths.append(find_width_at_half_max_height(x, psi_squared))
            width = (1 / initial_package_width * np.sqrt(initial_package_width ** 4 + (hbar / me * t) ** 2)) * 1.7
            package_expected_widths.append(width)

        transmission_calculator = TransmissionCalculator(
            disable_electron_potential_validation = True,
            _hbar = hbar,
            _me = me,
            package_wdh = initial_package_width,
            step_callback = _step_callback
        )
        
        # Act
        transmission_calculator.calculate_transmission(E, barrier, dx)
       
        error = max(np.absolute(np.array(package_widths) - np.array(package_expected_widths)))
     
        # Assert
        self.assertTrue(error < error_tolerance)

    def test_propability_density_is_1(self):
        E = 500 * codata.value("electron volt")
        V0 = 600 * codata.value("electron volt")

        dx = 0.1
        barrier = np.array(V0 + np.zeros(250))

        prob_dens = [] 
        error_tolerance = 0.1
  
        def _step_callback(self, psi, psi_squared, x, t, b, z):
            prob =  np.multiply(psi, psi.conj()).real
            prob_dens.append(dx*np.sum(psi_squared))

        transmission_calculator = TransmissionCalculator(
            step_callback = _step_callback
        )
        
        # Act
        transmission_calculator.calculate_transmission(E, barrier, dx)

        error = max(np.absolute(np.array(prob_dens) - 1))
        
        # Assert
        self.assertTrue(error < error_tolerance)

    def test_potential_to_energy_ratio(self):
        # Assemble
        V_0 = 20
        dx = 0.1
        barrier = np.array(V_0 + np.zeros(500))

        transmission_calculator = TransmissionCalculator(_me = 1, _hbar = 1, disable_electron_potential_validation = True)

        V_over_E = []
        transmissions = []

        for E in np.arange(0, 1.5 * V_0, 1):
            # Act
            transmission = transmission_calculator.calculate_transmission(E, barrier, dx)

            transmissions.append(transmission**2)

            V_over_E.append(V_0 / E)
 
        # print(transmissions)

        # plt.xlabel('V/E')
        # plt.ylabel('T^2')
        # plt.plot(V_over_E, transmissions)
        # plt.show()

        self.assertTrue(False)

        
if __name__ == '__main__':
    transmission_suite = unittest.TestLoader().loadTestsFromTestCase(TransmissionTestSuite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(transmission_suite)
