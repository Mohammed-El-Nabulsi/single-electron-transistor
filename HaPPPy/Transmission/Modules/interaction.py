#Input: energy in meV. internal calculations in meV and 10 femtoseconds.
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


from scipy.constants import codata

from unsure_A_potential002 import Potential
from B_Gauss00 import GaussianWave
from Fourier00 import Fourier
from D_Split_Step00 import Split_Step_Operator

me   = codata.value("electron mass energy equivalent in MeV") * 1e8 ;  # Convert to milli eV/c^2
hbar = codata.value("Planck constant over 2 pi in eV s")      * 1e19;  # Convert to 10femtoseconds*millielectronvolts

a = np.sqrt(me)/hbar

#new units
#sqrt_two_me_over_hbar = np.sqrt(2*me)/hbar
#print(sqrt_two_me_over_hbar)

#Main

#set testparameters
barriere = np.array([2,2,2,2,2])
dy = 1/2 #Abstand zwischen zwei sample points
energy = 1/3

#testing A_Potential00

potential = Potential(barriere, dy)


#print(potential.dx)
#print(potential.barrier)
#print(potential.gauss_index_width)
gauss_width = potential.gauss_width
N = potential.potential.size
#print(N)
#print(N)
#print()
pot = potential.potential
#print()
positions = potential.position_grid
#print(potential.pos_grid_width)
print(positions.size)
#print(potential.gauss_symmetry_index)
gauss_symmetry_point = potential.gauss_symmerey_point
#print(gauss_symmetry_point)


#finished, everything works fine here.

#testing interaction between gauss and potential class

gauss = GaussianWave(positions,gauss_symmetry_point,gauss_width,energy)

#print(gauss.width)
#print(gauss.symm)
print(gauss.k0)
#print(gauss.x_grid.size)
#print(gauss.x_probability)
#gauss.plot_x_package()
#plt.savefig(GaussPlot)
#gauss_package = gauss.x_package
#print(gauss_package)

#interaction test finished, everything works fine here.

#testing Fourier
fourier = Fourier(positions)
#print(fourier.get_wavenumbers().size)
#dft = fourier.dft(gauss.x_package)
#a = fourier.dft(gauss.x_package)
#fourier.plot_dft(gauss.x_package)
#fourier.plot_idft(a)

#print(abs(gauss.create_gauss_at_x(potential.gauss_symmerey_point)))

#dft2 = fourier.plot_dft(idft)
#print(dft.size)
#idft = fourier.plot_idft(dft)
#fourier.plot_dft(dft)
#k_pack = fourier.dft(gauss.x_package)
#fourier.plot_idft(k_pack)


#Anm.: Fourier Trafo des Gausspaktes ist anscheiend wieder Gauss, aber k-werte und amplithuden stimmen keinesfalls!

#Split Step Algorithm testing
splitstep = Split_Step_Operator(positions,pot)
#a = splitstep.first_step(gauss.x_package)
#
#indices = np.arange(0,200,1)
#
#psi = gauss.x_package
#for i in indices:
#    
#    i = i+1
#    psi = splitstep.first_step(psi)
##    psi21 =  np.multiply(psi,psi.conj()).real
##    plt.plot(fourier.k,psi21)#    plt.show()
#    
#    psi = splitstep.steps(psi)
#    psi = splitstep.steps(psi)
#    psi = splitstep.steps(psi)
#    psi = splitstep.steps(psi)
#    
#    psi = splitstep.final_step(psi)
#
#    psi2 = np.multiply(psi,psi.conj()).real
#    plt.plot(positions, psi2,positions,0.1*pot)
#    plt.xlabel('x')
#    plt.ylabel('$|\Psi|^2$')
#    
#    plt.show()
    #plt.savefig('i.png', bbox_inches='tight')

    
#done. everything wirks fine here!
    


print(positions.size)
#crrating an dynamic plot
   
fig, ax = plt.subplots()
line, = ax.plot(positions,gauss.x_package) #np.random.rand(10)
ax.set_ylim(0, 0.35)
#ax.set_xlim([0,positions[-1]])

#        plt.plot(positions, psi2,positions,0.1*pot)
plt.xlabel('x')
plt.ylabel('$|\Psi|^2$')

def update(data):
    line.set_ydata(data)
    return line,


def data_gen():
    indices = np.arange(0,1000,1)
    psi = gauss.x_package
    for i in indices:
        
        i = i+1
        psi = splitstep.first_step(psi)

        psi = splitstep.steps(psi)
        
        psi = splitstep.final_step(psi)
    
        psi2 = np.multiply(psi,psi.conj()).real

        yield psi2

ani = animation.FuncAnimation(fig, update, data_gen, interval=10)
plt.plot(positions,0.1*pot)
plt.show()


#
#"""
#A simple example of an animated plot
#"""
##import numpy as np
##import matplotlib.pyplot as plt
##import matplotlib.animation as animation
#
#fig, ax = plt.subplots()
#
#x = positions
#line, = ax.plot(x, np.sin(x))
#
#
#def animate(i):
#    
#    
#    line.set_ydata(np.sin(x + i/10.0))  # update the data
#    return line,
#
#
## Init only required for blitting to give a clean slate.
#def init():
#    line.set_ydata(np.ma.array(x, mask=True))
#    return line,
#
#ani = animation.FuncAnimation(fig, animate, np.arange(1, 200), init_func=init,
#                              interval=25, blit=True)
#plt.show()
