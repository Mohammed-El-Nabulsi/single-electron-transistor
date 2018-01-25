#the following units are being used throughout the HaPPPy-project: [energy]=meV, [length]=nm, [time]=ps, [temperature]=K

import math
import numpy as np
from HaPPPy.OneBody import OneBodySolver
from HaPPPy.TwoBody import TwoBodySolver
from HaPPPy.Transmission import TransmissionCalculator
from HaPPPy.MasterEquation import MasterEquationSolver

__docformat__ = 'reStructuredText'

class RateCalculator:
    
    """
    
    This class calculates the transition rates between the possible one-body
    and two-body states of an SET.

    It gives out the transition rates of events where
    electrons tunnel through a barrier (here through the Left barrier) onto the dot:
    
    :math:`\\Gamma^L_{\\alpha_N \\rightarrow \\beta_{N+1}} =
    \Gamma | \langle \\beta_{N+1} | \sum_{\\nu}`
    :math:`c_{\\nu}^{\\dagger} | \\alpha_N \\rangle |^2 n_F ( E_{\\beta_{N+1}} - E_{\\alpha_N} -\\mu_L )`

    where
    :math:`\\beta_{N+1}` is the final one-body or two-body state and :math:`\\alpha_N` is the initial vacuum state or one-body state of the dot.
    :math:`\\Gamma` includes the transmission coefficient of given tunneling event and the density of states function (DOS) on the dot
    (see the description below). :math:`n_F` is the fermifunction.
    
    :math:`c_{\\nu}^{\\dagger}` is the creation operator. Applied to the :math:`| \\alpha_N \\rangle` - state, it creates an
    artificial :math:`N + 1`-body state that is being projected onto the :math:`| \\beta_N \\rangle` - state.

    The rates of the transitions where the number of electrons on the dot decrease are similiarly calculated:

    :math:`\\Gamma^L_{\\alpha_N \\rightarrow \\beta_{N-1}} =
    \Gamma | \langle \\beta_{N-1} | \sum_{\\nu}`
    :math:`c_{\\nu} | \\alpha_N \\rangle |^2 ( 1 - n_F ( E_{\\beta_{N-1}} - E_{\\alpha_N} -\\mu_L ))`.  
    :math:`\\beta_{N-1}` and :math:`\\alpha_{N}` are, again, the final and initial states. Here we use the annihilation operator :math:`c_{\\nu}`.

    The Rates through the right barrier are calculated the same way (one simply exchanges :math:`\\mu_L` with :math:`\\mu_R`). 

    """
    
    def __init__(self):
        """ The constructor.
        """

        print("Hello from the RateCalculator, what can I do for you?")
     
    def doCalculation(self, E1, E2, muL, muR, T, pot, C, TCalc, Density, E0, L):

        """This function calculates the transition rates.

        :param E1: Array of eigenvalues of the single-particle-states on the dot
                   (calculated in the One Body Module).
        :type E1: np.array
        :param E2: Array of eigenvalues of the two-particle-states on the dot
                   (calculated in the Two Body Module).
        :type E2: np.array
        :param muL: Chemical potential of the source (main input).
        :type muL: float
        :param muR: Chemical potential of the drain (main input).
        :type muR: float
        :param T: Temperature (main input).
        :type T: float
        :param pot: Potential of the barriers.
        :type pot: matrix
        :param V: tunnel-barrier-potential in form of an array for input into :code:`TransmissionCalculator`
                  (main input)
        :type V: np.array
        :param C: array of coefficients-matrices (calculated in the Two Body Module), where each matrix
                  represents a two-particle-state in the single-particle-state-product-basis,
                  that must be sorted according to the order of the E2 array (E2[i] must be the energy
                  of the state represented by C[i]).
        :type C: np.array
        :param TCalc: class in the Transmission Module that contains the function
                      :code:`calculate_transmission(E,V)` which calculates the transmission coefficient
                      for a given energy difference E and tunnel-barrier-potential V.
        :type TCalc: class
        :param Density: function that gives us the density of states on source and drain (main input).
        :type Density: function
        :param E0: Offset Value/ Vacuum energy.
        :type E0: float
        :param L: Length of tunnel barriers.
        :type L: float
        
        :return: Returns (:math:`\\Gamma_L` , :math:`\\Gamma_R`):
                 :math:`\\Gamma_L` and :math:`\\Gamma_R` are two nxn-matrices
                 (:math:`n`=size(E1)+size(E2)+1) that contain the rates through the left barrier on
                 the "source-side"(Gamma_L) and through the right barrier on the "drain-side"(Gamma_R).
                 The both output matrices have following format:
                 Gamma=[[0,Gamma_01,zeros],
                       [Gamma_10,zeros,Gamma_12]
                       [zeros, Gamma_21,zeros]] 
        :rtype: list with two matrices
        
       
        """
        NEcut = len(E1) #we determine the number of single-particle states that we use
        VG=np.diag(pot)
        E= int(0.5*np.size(VG))
        V = VG[0:E] #since the potential of both barriers is symmetric and we only tunnel through one barrier, we only use one half of the potential
        dx= L/(np.size(pot))

        #Following prints are for debugging purposes:
        #print("---------------------------------------------------------------------")
        #print("---------------------------------------------------------------------")
        #print("Hier beginnt die Ausgabe von Rates:")
        #print("---------------------------------------------------------------------")
        #print("V:", V)
        #print("E1:", E1)
        #print("E2:", E2)
        #print("C:", C)

        kB=0.08629 #Boltzmann constant in meV/K
        
        
        def fermi(E,mu,T):
            """This fermi-function tells us with which likelyhood a state with given energy E is occupied on the lead
            E(float): energy difference between the initial and the final state that the tunneling electron has to carry
            mu(float): chemical potential of either drain(muR) or source(muL)
            T(float): temperature
            """
            f=1/(math.exp((E-mu)/(kB*T) )+1)
            return(f)
          

	#This function is called by the Gamma_ij-equations and includes the transmission-coefficient for each tunnelling-event
        #and the density of state function of the source and drain. 
        def Gamma(Ea,Eb,V):
             """:math:`\\Gamma` includes the transmission coefficient and DOS: :math:`\Gamma = | t |^2 * \\delta`

             Ea(float): energy of initial state
             Eb(float): energy of final state
             V(np.array): barrier potential
             """
             print(Ea)
             print(V)
             return (np.absolute(TCalc.calculate_transmission(Ea,V,dx)*Density.calculate_DensityofStates(np.absolute(Ea-Eb))))
                     
        #These next four functions are used to calculate the transition rates.Each function for a different ind of transition:
        #We distinguish between transitions, in which the number of electrons on the dot changes from one to two(Gamma_12) and reverse(Gamma_21).
        #And between transitions in which the number of electrons on the dot change from zero to one(Gamma_01) and reverse(Gamma_10).

        def Gamma_12(Ea,Eb,mu,T):
            """Calculates the rate of a transition from a one body state to a two body state

            Ea(float): energy of initial state
            Eb(float): energy of final state
            mu(float): chemical potential of either drain(muR) or source(muL)
            T(float): temperature
            """
            summe=0
            j=0
            Cb=C[np.where(E2==Eb)[0][0]]
            while j< NEcut:
                summe=Cb[np.where(E1==Ea)[0][0]][j]+summe
                j=j+1
            return(Gamma(Ea,Eb,V)*(np.absolute(summe))**2*fermi((Eb-Ea),mu,T))


        def Gamma_01(Eb,mu,T):
            """Calculates the rate of a transition from the vacuum state to a one body state

            Eb(float): energy of final state
            mu(float): chemical potential of either drain(muR) or source(muL)
            T(float): temperature
            """
            return(Gamma(0,Eb,V)*fermi((Eb-E0),mu,T))

        def Gamma_21(Ea,Eb,mu,T):
            """Calculates the rate of a transition from a two body state to a one body state

            Ea(float): energy of initial state
            Eb(float): energy of final state
            mu(float): chemical potential of either drain(muR) or source(muL)
            T(float): temperature
            """
            summe=0
            nu=0
            Ca=C[np.where(E2==Ea)[0][0]]
            while nu < NEcut:
                 summe=summe+Ca[np.where(E1==Eb)[0][0]][nu]
                 nu=nu+1
            return(Gamma(Ea,Eb,V)*(np.absolute(summe))**2*(1-fermi((Ea-Eb),mu,T)))

        def Gamma_10(Ea,mu,T):
            """Calculates the rate of a transition from a one body state to the vacuum

            Ea(float): energy of initial state
            Eb(float): energy of final state
            mu(float): chemical potential of either drain(muR) or source(muL)
            T(float): temperature
            """
            return(Gamma(Ea,0,V)*(1-fermi(Ea-E0,mu,T)))

        #creating the output matrices that later contain all the transition rates through either
        #the left or the right barrier
        Gamma_R=np.zeros((1+np.size(E1)+np.size(E2),1+np.size(E1)+np.size(E2)))
        Gamma_L=np.zeros((1+np.size(E1)+np.size(E2),1+np.size(E1)+np.size(E2)))

        #using a loop to fill the output matrices with transition rates.
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

        #print("Gamma_L und Gamma_R:")
        #print(Gamma_L,Gamma_R)
        #print("-----------------------------------------------------------------------")
        #print("---------------------------------------------------------------------")
        return(Gamma_L,Gamma_R)
        
