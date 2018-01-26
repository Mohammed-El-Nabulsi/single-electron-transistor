"""
.. module:: HaPPPy.MasterEquation
   :platform: Unix, Windows
   :synopsis: Tackles the master equation.

A module to simulate the time development of propabilities :math:`\\vec{P}(t)`
or to find stationary solutions :math:`\\vec{P}_{stat}` for an :math:`n`-state
quantum dot transistor with two reservoirs with transition rates
:math:`\\Gamma^{L}` and :math:`\\Gamma^{R}` including the netto current
:math:`I(t)` or the sationary netto current :math:`I_{stat}` respectively.

This module contains two classes. MasterEquationSolver does all mentioned
calculations and Simulation is a special return type for simulations produced by
MasterEquationSolver to provide a simple yet convenient object to operate with.
"""

from .MasterEquationSolver import MasterEquationSolver
from .MasterEquationSolver import Simulation
