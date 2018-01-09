from scipy import constants
import matplotlib.pyplot as plt
import numpy as np
import h5py
from math import sqrt

n = 1000  # n-gridpoints via input at later stage; number of elements the lenght is divided trough user input
l = 100  # l-lenght of potential pot  length of the potential user input in nm runs best 10nm-99nm

# creating array in user input length l and user input elements n and not changeable array for plot x-axis
# creating transposed array a -> at
a = np.linspace((-l / 2), (l / 2), n)
aplot = np.linspace((-l / 2), (l / 2), n)


# validating constants for unit calculation of hamilton matrix
hbar = constants.hbar
me = constants.m_e
e = constants.e

print("hbar: ", hbar)
print("me: ",me)
print("e: ",e)

# creating potential matrix with user input elements
pot_mat = np.zeros((n, n))


# function to build potential matrix out of user adjusted array
def mat_build(a):
    i = 0
    while i < n:
        pot_mat[i, i] = a[i]
        i += 1
    return pot_mat


# function to plot potential matrix pot_mat
def plot_build(a):
    pot_plot = plt.plot(aplot, a)
    plt.title("Potential Matrix")
    plt.ylabel("energy [meV]")
    plt.xlabel("x [A]")
    plt.show()
    return pot_plot


# function to build gauss potential
def gaussian(x, A, sigma):
    return A*(-np.exp(-np.power((x) / sigma, 2.) / 2.))

# user potential decision and adjustment of the array
# which is then translated to potential matrix pot_mat
user = 0
unit_pot = 1.0
while user != 1 or user != 2 or user != 3 or user != 4:
    print("horizontal \t ==> 1\nlinear \t\t ==> 2\nsquare \t\t ==> 3\ngauss \t\t ==> 4")
    user = int(input("Press 1, 2, 3 or 4 for the listed potentials   "))
    if user == 1:
        print("\nYou've choosen a horizontal potential. \nnice choice!")
        hori = float(input("How big is your potential? [meV]  "))
        for x in np.nditer(a, op_flags=['readwrite']):
            x[...] = hori
        mat_build(a)  # build pot_mat
        plot_build(a)  # plot the function
        break
    elif user == 2:
        print("\nYou've choosen a linear potential. \nnice choice!")
        intersection = float(input("Where does your potential cross zero? [A]  "))
        slope = float(input("With which slope does your potential rise? [meV/A]  "))
        for x in np.nditer(a, op_flags=['readwrite']):
            x[...] = slope*x + intersection
        mat_build(a)  # build pot_mat
        plot_build(a)  # plot the function
        break
    elif user == 3:
        print("\nYou've choosen a square potential. \nnice choice!")
        intersection = float(input("Where does your potential cross zero? [A]  "))
        square = float(input("Which exponent for potential? [meV/A]   "))
        slope = float(input("With which slope does your potential rise? [meV/A]  "))
        for x in np.nditer(a, op_flags=['readwrite']):
            x[...] = slope * (x**square) + intersection
        mat_build(a)  # build pot_mat
        plot_build(a)  # plot the function
        # unit system and calculation of final hamiltonian matrix ham_mat
        # factor 1000 for the unit system in order to reach meV
        unit_pot = ((1/2)*me*(3.5174*(10**29))*((10**-9))**2)/e  # potential matrix in meV
        break
    elif user == 4:
        print("\nYou have chosen the Gauss potential!")
        A = float(input("Please enter A in meV: "))
        sigma = float(input("Please enter sigma: "))
        for x in np.nditer(a, op_flags=['readwrite']):
            x[...] = gaussian(x, A, sigma)
        mat_build(a)  # build pot_mat
        plot_build(a)  # plot the function
        unit_pot = A
        break
    else:
        print("=" * 30 + "\nPlease enter valid number!\n" + "=" * 30)
        user = True

# creating body of the kinetic matrix with second derivate of the location
kin_mat = np.zeros((n, n))

i = 0
while i < n:
    kin_mat[i, i] = 2
    i += 1

i = 0
while i < n-1:
    kin_mat[i, i+1] = kin_mat[i+1, i] = -1
    i += 1
print(kin_mat)

# unit system and calculation of final hamiltonian matrix ham_mat
# factor 1000 for the unit system in order to reach meV
# unit_pot = ((1/2)*me*(3.5174*(10**29))*((10**-9))**2)/e  # potential matrix in meV

unit_kin = ((hbar**2)*1000)/(2*me*(10**-18)*e)
print(unit_kin, "\n", unit_pot)  # control print for unit
# dx for the derivate of the matrix
dx = l/n
# build the final hamilton matrix ham_mat
ham_mat = unit_kin * kin_mat * (1/(dx*dx)) + unit_pot * pot_mat

# calculate eigenvalues (stored in la) and eigenvectors (stored in v)
la, v = np.linalg.eigh(ham_mat)

# printing eingenvalues and eigenvectors
# as option for debugging
for i in range(10):
    print("Eigenvalue:\t\t ", la[i])  # , "\ncorresponding eigenvector: ", v[:,i] )

# n_plot is x-axis for plotting the eigenvectors
n_plot = np.linspace(0, n, n)

# n_plot vs. square of eigenvectors to archieve the probabilty
# range 0,2 gives first two wavefunctions
# give option for debugging
for i in range(0,3):
    test = plt.plot(n_plot, (v[:,i]*v[:,i]))
    plt.title("Wavefunction")
    plt.ylabel("density")
    plt.xlabel(" ")
plt.show()


norm = 0.0
for i in range(n):
    norm += (v[i,0]*v[i,0]) * dx

sqrt_norm = sqrt(norm)
print(sqrt_norm)

v = v / sqrt_norm # normierte matrix von eigenvektoren die soll übergeben werden

norm = 0.0
for i in range(n):
    norm += (v[i,25]*v[i,25]) * dx
print(norm)


# export data as hdf5 file
data = h5py.File("data_group1.hdf5", "w")
dsetla = data.create_dataset("eigenvalues_group1", (n,), dtype='f')  # wenn er wegen n palaber macht mal mit ner zahl anstatt n versuchen, zb 100
dsetv = data.create_dataset("eigenvectors_group1", (n,), dtype='f')
# shape of dataset: n; datatype ist f für float ( bin mir nicht sicher ob f für float richtig ist, für integer nehmen die im inet aber 'i')
dsetla[...] = la
dsetv[...] = v
# la als dataset abspeichern, müssen das noch besser benennen und abchecken was gruppe 2 und 3 alles braucht
# hab hier kein scipy und hdf5 installiert desswegen kann ich nciht abchecken obs funktioniert
# obacht: import h5py ist rausgehachtagt oben !

# x = np.sort(la)

# print(x)

#with h5py.File('name-of-file.h5', 'w') as hf:
#    hf.create_dataset("name-of-dataset",  data=v[:])