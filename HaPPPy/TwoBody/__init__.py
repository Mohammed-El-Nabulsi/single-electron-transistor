"""
.. module:: HaPPPy.TwoBody
   :platform: Unix, Windows
   :synopsis: Tackles the two body program.

This finds eigenstates to the two-electron Hamiltonian :math:`H_{e^-e^-}`. 

.. math:: \hat{H}_{e^-e^-} = \hat{H}_{e^-}\otimes\\text{id}+\\text{id}\otimes \hat{H}_{e^-} + \hat{V}_\\text{int}
    
with

.. math:: \hat{V}_\\text{int}\Psi_{e^-e^-}(\\vec{r}_1,\\vec{r}_2) = \\frac{\\alpha_\\text{Coulomb}}{|\\vec{r}_1-\\vec{r}_2|} \Psi_{e^-e^-}(\\vec{r}_1,\\vec{r}_2)

For any basis of one-electron wave functions :math:`{\phi_i\}`, :math:`{\phi_i \otimes \phi_j\}` is a basis of the two-electron wave functions. Here :math:`{\phi_i}` was chosen as the basis of one-electron eigenfunctions. Hence, each two-electron eigenfunction can be written as

.. math:: \Psi^i = \sum_{j,k} Q^i_{j,k}\phi_i\otimes\phi_j

with some coefficient matrix :math:`Q^i`. These matrices are, given the one-electron eigenfunctions, a complete description of the two-electron eigenfunctions.
"""

from .TwoBodySolver import TwoBodySolver
from .OneParticleLoader import SpectrumData
from .TwoParticleLoader import TwoBodySpectrumData
