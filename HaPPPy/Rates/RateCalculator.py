#the following units are being used throughout the HaPPPy-project: [energy]=meV, [length]=nm, [time]=ps, [temperature]=K

#We must import all modules of the other parts of the HaPPPy-project we use.
import math
import numpy as np
from HaPPPy.OneBody import OneBodySolver
from HaPPPy.TwoBody import TwoBodySolver
from HaPPPy.Transmission import TransmissionCalculator
from HaPPPy.MasterEquation import MasterEquationSolver

#V=Vinput[:half]
#E2 = Energie_Zweiteilchen 
#C= Koeffizienten_Zweiteilchen #TBSolver.doCalculation()
#muL, muR,T, Vinput aus main importieren 
#Transimissionsmatrix wird in def Gamma importiert

class RateCalculator:
    """This class solves the rate equation for every possible transmission.
    """
    
    def __init__(self):
        """ The constructor.
        """

        print("Hello from the RateCalculator, what can I do for you?")
     
    def doCalculation(self, E1, E2, muL, muR, T, pot, C, TCalc, Density, E0):

        """This class is responsible for calculating all rates of transmissions of electrons into the dot and from the dot.
        It needs the following inputs:
        E1 (np.array): array of eigenvalues of the single-particle-states on the dot (calculated in the One Body Module)
        E2 (np.array): array of eigenvalues of the two-particle-states on the dot (calculated in the Two Body Module)
        muL (float): chemical potential of the source (main input)
        muR(float): chemical potential of the drain (main input)
        T (float): temperature (main input)
        V (np.array): tunnel-barrier-potential in form of an array for input into TransmissionCalculator (main input)
        C (np.array): array of coefficients-matrices (calculated in the Two Body Module), where each matrix represents a two-particle-state in the single-particle-state-product-basis, that must be sorted according to the order of the E2 array (E2[i] must be the energy of the state represented by C[i])
        TCalc (class): class in the Transmission Module that contains the function "calculate_transmission(E,V)" which calculates the transmission coefficient for a given energy difference E and tunnel-barrier-potential V
        Density (function): function that gives us the density of states on source and drain (main input)

        And gives the following output:
        (Gamma_L, Gamma_R): Gamma_L and Gamma_R are two nxn-matrices (n=size(E1)+size(E2)+1) that contain the rates through the left barrier on the "source-side"(Gamma_L) and through the right barrier on the "drain-side"(Gamma_R), ready to be read out by the Master Equation Module.
        Both matrices have the same format and only differ in their values, due to the different chemical potentials in source and drain.

        Each matrix has following format:
        Gamma=[[0,Gamma_01,zeros],
               [Gamma_10,zeros,Gamma_12]
               [zeros, Gamma_21,zeros]] 
        """
        NEcut = len(E1)
        VG=np.diag(pot)
        E= int(0.5*np.size(VG))
        V = VG[0:E]
        print(V)
        kB=0.08629 #Boltzmann constant in meV/K
        #The fermi-function is called by the Gamma_ij-equations (Gamma_12,Gamma_21,Gamma_01,Gamma_10) that calculate our rates.

      
        def fermi(E,mu,T):
            f=1/(math.exp((E-mu)/(kB*T) )+1)
            return(f)
          

	#This function is called by the Gamma_ij-equations and includes the transmission-coefficient for each tunnelling-event and the density of state-function of the source and drain. 


        def Gamma(Ea,Eb,V):
             print(Ea)
             print(V)
             return (np.absolute(TCalc.calculate_transmission(Ea,V,0.1)*Density.calculate_DensityofStates(np.absolute(Ea-Eb))))
                     
        #These next four functions are used to calculate the rates of transmissions.
        #Each function represents a specific kind of equation that deals with a certain kind of transmission;
        #We distinguish between transmissions, where the number of electrons on the dot changes from one to two(Gamma_12) and reverse(Gamma_21).
        #And between transmissions where the number of electrons on the dot change from zero to one(Gamma_01) and reverse(Gamma_10).

        def Gamma_12(Ea,Eb,mu,T): 
            summe=0
            j=0
            Cb=C[np.where(E2==Eb)[0][0]]
            while j< NEcut:
                summe=Cb[np.where(E1==Ea)[0][0]][j]+summe
                j=j+1
            return(Gamma(Ea,Eb,V)*(np.absolute(summe))**2*fermi((Eb-Ea),mu,T))


        def Gamma_01(Eb,mu,T): 
            return(Gamma(0,Eb,V)*fermi((Eb-E0),mu,T))

        def Gamma_21(Ea,Eb,mu,T): 
            summe=0
            nu=0
            Ca=C[np.where(E2==Ea)[0][0]]
            while nu < NEcut:
                 summe=summe+Ca[np.where(E1==Eb)[0][0]][nu]
                 nu=nu+1
            return(Gamma(Ea,Eb,V)*(np.absolute(summe))**2*(1-fermi((Ea-Eb),mu,T)))

        def Gamma_10(Ea,mu,T): 
            return(Gamma(Ea,0,V)*(1-fermi(Ea-E0,mu,T)))


        Gamma_R=np.zeros((1+np.size(E1)+np.size(E2),1+np.size(E1)+np.size(E2)))   #defining the matrix that contains all rates of transmissions through the right barrier
        Gamma_L=np.zeros((1+np.size(E1)+np.size(E2),1+np.size(E1)+np.size(E2)))   #defining the matrix that contains all rates of transmissions through the left barrier

        Gamma_R=np.zeros((1+np.size(E1)+np.size(E2),1+np.size(E1)+np.size(E2)))
        Gamma_L=np.zeros((1+np.size(E1)+np.size(E2),1+np.size(E1)+np.size(E2)))

        i_=0
        for i in E1:
                j_=0
                for j in E2:
                        j_ =np.where(E2==j)[0][0]
                        Gamma_L[i_+1][j_+1+np.size(E1)]=Gamma_12(i,j,muL,T)
                        Gamma_L[j_+1+np.size(E1)][i_+1]=Gamma_21(j,i,muL,T)
                        Gamma_R[i_+1][j_+1+np.size(E1)]=Gamma_12(i,j,muR,T)
                        Gamma_R[j_+1+np.size(E1)][i_+1]=Gamma_21(j,i,muR,T)
                        j_=j_+1
                Gamma_L[0][i_+1]=Gamma_10(i,muL,T)
                Gamma_R[0][i_+1]=Gamma_10(i,muR,T)
                Gamma_L[i_+1][0]=Gamma_01(i,muL,T)
                Gamma_R[i_+1][0]=Gamma_01(i,muR,T)
                i_=1+i_
        print(Gamma_L,Gamma_R)
        return(Gamma_L,Gamma_R)
        
