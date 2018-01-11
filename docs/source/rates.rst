Rates Module
============

we need:

from group  1: array of energy-eigenvalues of the single-particle-eigenstates

from group  2: array of coefficients-matrices, where each matrix represents a two-particle-state in the product-basis,
array of energy-eigenvalues of two-particle-states in the same order as the array of coefficient-matrices
tunnel-potential in the form of an array

from group  3: function that gives us the transmission-coefficients

from main: temperature, chemical potential L&R

our output:

array of rate-matrices sorted in the following order: (Gamma_01L,Gamma_01R,Gamma_10L,Gamma_10R,Gamma_12L,Gamma_12R,Gamma_21L,Gamma_21R)


.. automodule:: HaPPPy.Rates
.. autoclass:: RateCalculator
