

import numpy as np
from scipy.fftpack import fft, ifft
import cmath as cm

#define parameter
h = 1.054571800 * 10^-34
m =  9,109 383 56 * 10^-31
dt = 12

#Example
#getting gaussian wave as array
psi = np.array([1.0, 2.0, 1.0, -1.0, 1.5])

#changing into k-space
psi_k = fft(psi)





#importing pulse-eigenvalues
p_op_EWs = np.array([2.0, 2.0, 2.0, 2.0, 2.0])

#creating the first half of pulse operator
p_op = np.array(cm.exp(-dt/2 * h / (2*m) * (p_op_EWs)^2 j)
#use the p-op one our array
psi_k_neu = psi_k * p_op




#return to space
psi_x = ifft(psi_k_neu)

#define the potewntials
Vi = list([2.0, 2.0, 2.0, 2.0, 2.0])

#creating the space operator
x_op = np.array(cm.exp(-dt -Vi*dt/h j))

#use the x-op on the wave funktion
Psi_x_neu = psi_x * x_op




#changeing back to k-space
Psi_k2 = fft(Psi_x_neu)

#use the p operator again
psi_k_neu2 = psi_k * p_op

#remaining problems: operators?, how to loop, counter/break
