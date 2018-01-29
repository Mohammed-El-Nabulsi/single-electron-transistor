from HaPPPy.Transmission.TransmissionCalculator import TransmissionCalculator

from threading import Thread
from scipy.constants import codata

import threading
import unittest
import pytest
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation

_x = []
_psi_plot = []
_V = []

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

    def xtest_width_of_free_gaussian_package_grows_correctly(self):
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
            half_max = max(Y) / 2
            #find when function crosses line half_max (when sign of diff flips)
            #take the 'derivative' of signum(half_max - Y[])
            d = np.sign(half_max - np.array(Y[0:-1])) - np.sign(half_max - np.array(Y[1:]))
            #plot(X,d) #if you are interested
            #find the left and right most indexes
            left_idx = np.argwhere(d > 0)[0]
            right_idx = np.argwhere(d < 0)[-1]

            return (X[right_idx] - X[left_idx])[0] / (2 * np.sqrt(2 * np.log(2))) #return the difference (full width)

        def _step_callback(self, psi, psi_plot, x, n, finished):
            t = self.dt * n
            package_widths.append(find_width_at_half_max_height(x, psi_plot))
            width = initial_package_width/2 * np.sqrt(1 + (4*(hbar**2)*(t**2))/((me**2)*(initial_package_width**4)))
            package_expected_widths.append(width)

        def _step_exit(self, n):
            if (n == 1000):
                return True

            return False

        transmission_calculator = TransmissionCalculator(
            disable_electron_potential_validation = True,
            _hbar = hbar,
            _me = me,
            package_wdh = initial_package_width,
            step_callback = _step_callback,
            step_exit = _step_exit
        )


        # Act
        transmission_calculator.calculate_transmission(E, barrier, dx)
       
        error = max(np.absolute(np.array(package_widths) - np.array(package_expected_widths)))

        # Assert
        self.assertTrue(error < error_tolerance)

    def test_propability_density_is_1(self):
        E = 500 * codata.value("electron volt") * 1e-3
        V0 = 600 * codata.value("electron volt") * 1e-3

        dx = 0.1
        barrier = np.array(V0 + np.zeros(250))

        prob_dens = [] 
        error_tolerance = 0.1
  
        def _step_callback(self, psi, psi_squared, x, n, finished):
            if (finished == True):
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

    def xtest_transmission_to_energy_ratio(self):
        # Assemble
        E = 100 * codata.value("electron volt") * 1e-3
        V_0 = 400 * codata.value("electron volt") * 1e-3

        # E = 10 * codata.value("electron volt") * 1e-3
        # V_0 = 40 * codata.value("electron volt") * 1e-3

        dx = 0.1

        barrier = np.array(V_0 + np.zeros(50))

        V_over_E = []
        transmissions = []

        def _step_callback(self, psi, psi_plot, x, n, finished):
            if (finished == True):
               plt.xlabel('x in pm')
               plt.ylabel('$|\Psi(x)|^2$')
               plt.plot(
                       self.x,
                       psi_plot,
                       self.x,
                       max(psi_plot)*(self.V)/max(self.V))
               plt.show()

        transmission_calculator = TransmissionCalculator(disable_electron_potential_validation = True,
            # step_callback = _step_callback
        )

        for E in np.arange(0, 4 * V_0, V_0 / 5):
            # Act
            if (E == 0):
                E = 1 * codata.value("electron volt") * 1e-3

            transmission = transmission_calculator.calculate_transmission(E, barrier, dx)


            transmissions.append(transmission**2)

            V_over_E.append(E / V_0)
 
        # print(transmissions)

        plt.xlabel('E/V')
        plt.ylabel('T^2')
        plt.plot(V_over_E, transmissions)
        plt.show()

        self.assertTrue(True)

    def xtest_transmission_returns_realistic_values(self):
        E = 500 * codata.value("electron volt") * 1e-3
        V0 = 600 * codata.value("electron volt") * 1e-3

        dx = 0.05
        barrier = np.array(V0 + np.zeros(250))

        # E = 10 * codata.value("electron volt") * 1e-3
        # V0 = 40 * codata.value("electron volt") * 1e-3

        # dx = 0.1
        # barrier = np.array(V0 + np.zeros(3000))

        prob_dens = [] 
        error_tolerance = 0.1

        def _step_callback(self, psi, psi_plot, x, n, finished):
            if (finished == True):
               plt.xlabel('x in pm')
               plt.ylabel('$|\Psi(x)|^2$')
               plt.plot(
                       self.x,
                       psi_plot,
                       self.x,
                       max(psi_plot)*(self.V)/max(self.V)
               )

               plt.show()

        transmission_calculator = TransmissionCalculator(
            # step_callback = _step_callback,
        )
        
        # Act
        transmission = transmission_calculator.calculate_transmission(E, barrier, dx)
        
        # Assert
        self.assertTrue(np.abs(1 - transmission / 0.175) < error_tolerance) 
        
    def xtest_animate_progress(self):
        E = 550 * codata.value("electron volt") * 1e-3
        V0 = 600 * codata.value("electron volt") * 1e-3

        dx = 0.05
        barrier = np.array(V0 + np.zeros(200))

        prob_dens = [] 
        error_tolerance = 0.1

        fig, ax = plt.subplots()
        line = ax.plot([], [])

        plt.xlabel('x')
        plt.ylabel('$|\Psi|^2$')

        def update(data):
            plt.gca().cla() 
            plt.plot(_x, 0.6  * (_V / max(_V)))
            ax.plot(_x, _psi_plot),

            ax.set_xlim(-100, 100)
            ax.set_ylim(0, 0.8)

        def data_gen():
            global _psi_plot

            yield _psi_plot

        def _step_callback(self, psi, psi_plot, x, n, finished):
            global _psi_plot
            global _x
            global _V

            _psi_plot = psi_plot
            _x = self.x
            _V = self.V

        transmission_calculator = TransmissionCalculator(
            step_callback = _step_callback,
        )

        ani = animation.FuncAnimation(fig, update, data_gen, interval=30)

        # Act
        Thread(target = transmission_calculator.calculate_transmission, args = (E, barrier, dx)).start()
        
        plt.show()

        self.assertTrue(True)

if __name__ == '__main__':
    transmission_suite = unittest.TestLoader().loadTestsFromTestCase(TransmissionTestSuite)
    unittest.TextTestRunner(verbosity=2, buffer=True).run(transmission_suite)
