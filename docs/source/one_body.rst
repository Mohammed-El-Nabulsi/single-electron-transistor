One Body Module
===============

Introduction
------------

The *OneBodyModule* is written by:

    - Gustaf Beier
    - Lars Klemeyer
    - Sani Harouna-Mayer
    - Ludwig Hendl

For any further information contact the developer team via ludwig@4xhendl.de

The One Body Module is written to calculate single bound electron states for the three
potentials listed below.

    1. Harmonic Potential

    2. Box Potential

    3. Gauss Potential

Choosing one of the potentials the user yields the one particle energy levels as well
as the corresponding wave functions. Plots shown in the documentation are either the chosen
potential or the probability density of a wave function (most likely of the first 3 wave functions).

Enjoy the following documentation.


Physics
-------

To calculate the energy states of a given potential, the Schrödinger quotation is needed.

        .. math::
            \\Eta =
            \\begin{pmatrix}
                \\vec{P}_1 & \\vec{P}_2 & \\dots & \\vec{P}_n
            \\end{pmatrix}
			
"schrödi"



https://pythonhosted.org/an_example_pypi_project/sphinx.html#headers
	
.. automodule:: HaPPPy.OneBody
.. autoclass:: OneBodySolver
