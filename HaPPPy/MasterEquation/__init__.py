"""
.. module:: HaPPPy.MasterEquation
   :platform: Unix, Windows
   :synopsis: Tackles the master equation.

A module to simulate the time development of propabilities :math:`\\vec{P}(t)`
for an :math:`n`-state system with transition rates :math:`\\Gamma`
including the netto current :math:`I(t)`.

See documentation of MasterEquationSolver.doCalculation for details.
"""

from .MasterEquationSolver import MasterEquationSolver
from .MasterEquationSolver import Simulation
