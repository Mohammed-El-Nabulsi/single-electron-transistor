One Body Module
===============
The *OneBodyModule* is written by:

    - Gustaf Beier
    - Lars Klemeyer
    - Sani Harouna-Mayer
    - Ludwig Hendl

For any further information contact the developer team via: ludwig@4xhendl.de

Introduction
------------

The One Body Module is the first part of the HaPPPy software project. For further information see **The HaPPPy program**. 
This module is written to calculate the single bound electron states in a quantum dot described by a given potential.

    .. figure::  _static/pic1.jpg
       :align:   center

**Figure 1:**    Drain and source as well as the quantum dot with single electron states.

The energy states for the quantum dot can be calculated for a *harmonic potential*, *box potential* and a *gauss potential*.

Choosing one of the potentials, the user yields the one particle energy levels as well
as the corresponding wave functions. Plots shown in the documentation are either the chosen
potential or the probability density of a wave function (most likely of the first 3 wave functions).

Enjoy the following documentation.


Physics
-------

To calculate the energy states of a given potential, the Schrödinger quotation is needed.

.. math::
    \widehat{H}\psi = E\psi

*H* is the Hamilton operator which is the sum of the kinetic and the potential term.

.. math::
    \widehat{H} = \widehat{K} + \widehat{V}
.. math::
    \widehat{H} = -\frac{\hbar^2}{2m} \frac{d^2}{dx^2} + V(x)
	
A analytic solution to this problem is only in a few cases possible. Therefore the described problem should be calculated numerical.
In order to do that the Hamilton operator must be translated to a matrix representation.  
A potential with the length *l* describes the quantum dot. To calculate the potential numerical from -*l*/2 to *l*/2 
the length is divided into n grid points. As more grid points are used, the approximation of the result becomes more and more accurate. 
The delta *x* between two grid points is the quotient of *l* and *n*.

.. math::
    \Delta x = \frac{l}{n} 

Potential term
++++++++++++++

The potential term generates to each *x* value of the quantum dots length *l* a *y* value of a function describing the potential.
Those values are then stated on the diagonal of the potential matrix.

.. math::
    \widehat{V} = \begin{pmatrix} V(x_0)&0&0&0\\0&V(x_2)&0&0\\0&0&\ddots&0\\0&0&0&V(x_{n-1})\end{pmatrix}

Kinetic term
++++++++++++

The kinetic term generates the second derivate of the wave function and multiplies it with constants. To yield the matrix representation, 
the second derivate of the differential operator must be calculated as shown below (the function *f* is just an example).

.. math::
    \frac{df}{dx} = \frac{f(x_{n+1}) - f(x_n)}{\Delta x}
.. math::
    \frac{d^2f}{dx^2} = \frac{1}{\Delta x}\biggl\{\frac{df}{\Delta x}(x_{n+1})-\frac{df}{\Delta x}(x_n)\biggl\}
.. math::
    \frac{d^2f}{dx^2} = \frac{1}{\Delta x^2}\big\{\big(f(x_{n+2})-f(x_{n+1})\big)-\big(f(x_{n+1})-f(x_n)\big)\big\}
.. math::
    \frac{d^2f}{dx^2} = \frac{1}{\Delta x^2}\big\{f(x_{n+2})-2f(x_{n+1})+f(x_n)\big\}

With this results the kinetic matrix can be build as:

.. math::
    \widehat{K} = -\frac{\hbar^2}{2m} \begin{pmatrix} -2&1&0&0&0\\ 1&-2&1&0&0\\0&1&-2&1&0\\0&0&\ddots&\ddots&\ddots\\0&0&0&1&-2 \end{pmatrix} \frac{1}{\Delta x^2}

Hamilton matrix
+++++++++++++++
	
The complete hamilton matrix is then:

.. math::
    \widehat{H} = -\frac{\hbar^2}{2m} \begin{pmatrix} -2&1&0&0&0\\ 1&-2&1&0&0\\0&1&-2&1&0\\0&0&\ddots&\ddots&\ddots\\0&0&0&1&-2 \end{pmatrix} \frac{1}{\Delta x^2} + \begin{pmatrix} V(x_0)&0&0&0&0\\0&V(x_2)&0&0&0\\0&0&V(x_3)&0&0\\0&0&0&\ddots&0\\0&0&0&0&V(x_{n-1})\end{pmatrix}
	
With the matrix representation of the Hamilton operator the Schrödinger equation becomes a eigenvalue problem. The calculations 
yields the eigenvalues (energy of the states)and the corresponding eigenvectors (wave function of the energy state). A plot of 
the squared wave function of a eigenvalue shows the probability density.
    - energy *E* --> eigenvalues
	- wf psi     --> eigenvectors

Unit convention
---------------

The length of the potential must be given in nano meters:

.. math::
    [l] = nm

The energy of the single bound electron states are needed in meV:

.. math::
    [E] = meV

Joule is the unit of the calculated energy levels. In order to reach meV, the unit of the kinetic matrix as well as the unit 
of the potential matrix are devided by the electron charge *e* and multiplied by 1000. 
	
https://pythonhosted.org/an_example_pypi_project/sphinx.html#headers
	
	
Code realization
----------------

.. automodule:: HaPPPy.OneBody
.. autoclass:: OneBodySolver
