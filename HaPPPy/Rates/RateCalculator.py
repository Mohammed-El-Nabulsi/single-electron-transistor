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

    This class returns all transition rates belonging to events in which only one electron tunnels through a barrier.
    We assume there are only up to 2 electrons on the dot. Therefore we only consider one-body states,
    two-body states and the vacuum state.

    The rate of events in which an electron tunnels through the left barrier onto the dot and changes the dot-state from :math:`| \\alpha_N \\rangle` to :math:`| \\beta_{N+1} \\rangle` looks like this:
    
    :math:`\\Gamma^L_{\\alpha_N \\rightarrow \\beta_{N+1}} =
    \Gamma | \langle \\beta_{N+1} | \sum_{\\nu}`
    :math:`c_{\\nu}^{\\dagger} | \\alpha_N \\rangle |^2 n_F ( E_{\\beta_{N+1}} - E_{\\alpha_N} -\\mu_L )` ,

    where
    :math:`\\beta_{N+1}` is the final one-body or two-body dot-state and :math:`\\alpha_N` is the initial vacuum state or one-body state of the dot.
    :math:`\\Gamma` includes the transmission coefficient of given tunneling event and the density of states function (DOS)
    on the dot
    (see the description below). :math:`n_F` is the fermifunction.
    
    :math:`c_{\\nu}^{\\dagger}` is the creation operator. Applied to the :math:`| \\alpha_N \\rangle` - state, it creates an
    artificial :math:`N + 1`-body state *(if N=0: it would be a physical one-body state,
    if N=1: it would be a naive product of the* :math:`| \\alpha_N \\rangle` *-state and another one-body state
    with quantum-number* :math:`\\nu` *)* ,
    that is being projected onto the  final :math:`| \\beta_{N+1} \\rangle` eigenstate.
    :math:`E_{\\alpha_N}` and :math:`E_{\\beta_{N+1}}`
    are the energies of the initial and final states.

    The rates of the transitions in which the number of electrons on the dot decrease are similiarly calculated:

    :math:`\\Gamma^L_{\\alpha_N \\rightarrow \\beta_{N-1}} =
    \Gamma | \langle \\beta_{N-1} | \sum_{\\nu}`
    :math:`c_{\\nu} | \\alpha_N \\rangle |^2 ( 1 - n_F ( E_{\\beta_{N-1}} - E_{\\alpha_N} -\\mu_L ))`.  
    :math:`\\beta_{N-1}` and :math:`\\alpha_{N}` are, again, the final and initial states. Here we use the annihilation operator :math:`c_{\\nu}`.

    The rates of events in which electrons tunnel through the **right barrier** are calculated the same way (one simply exchanges :math:`\\mu_L` with :math:`\\mu_R`). 

    """
    
    def __init__(self):
        """ The constructor.
        """

        print("Hello from the RateCalculator, what can I do for you?")
     
    def doCalculation(self, E1, E2, muL, muR, T, pot, C, TCalc, Density, E0, L):

        """This function calculates the transition rates.
        It contains a number of "sub-function" to break down the :math:`\\Gamma^L_{\\alpha \\rightarrow \\beta}` - equation.

        :param E1: Array of energy-eigenvalues of the single-particle states on the dot
                   (calculated by the One Body Module).
        :type E1: np.array
        :param E2: Array of eigenvalues of the two-particle-states on the dot
                   (calculated by the Two Body Module).
        :type E2: np.array
        :param muL: Chemical potential of the source (main input).
        :type muL: float
        :param muR: Chemical potential of the drain (main input).
        :type muR: float
        :param T: Temperature (main input).
        :type T: float
        :param pot: Potential of the barriers (main input).
        :type pot: matrix
        :param C: array of coefficients-matrices (calculated in the Two Body Module), where each matrix
                  represents a two-body state in the one-body state product basis,
                  that must agree to the order of the E2 array *(E2[i] must be the energy
                  of the state represented by C[i])*.
        :type C: np.array
        :param TCalc: class in the Transmission Module that contains the function
                      :code:`calculate_transmission(E,V)` which calculates the transmission coefficient
                      for a given energy difference E and tunnel-barrier-potential V.
        :type TCalc: class
        :param Density: Function that gives us the density of states on source and drain (main input).
        :type Density: function (with one parameter: energy (float))
        :param E0: Offset Value/ Vacuum energy (main input).
        :type E0: float
        :param L: length of the potential (main input).
        :type L: float
        
        :return: Returns (:math:`\\Gamma_L` , :math:`\\Gamma_R`):
        
                 :math:`\\Gamma_L` and :math:`\\Gamma_R` are two :math:`(n*n)`-matrices
                 :math:`(n = n1 + n2 + 1)`, *(where* :math:`n1` *and* :math:`n2` *are the lengths of the E1 -and E2-arrays)*.

                 They contain the rates through the left barrier on the source-side (:math:`\\Gamma_L`)
                 and through the right barrier on the drain-side (:math:`\\Gamma_R`).
                 Both output matrices have following format:
                 

                 .. math::
                     \\Gamma_{L or R}
                     =
                     \\begin{pmatrix}
                         0                       & \\Gamma_{10} (E1_1)          & \\dots  & \\Gamma_{10} (E1_{n1})
                             & 0                           & \\dots  & 0                                  \\\\
                         \\Gamma_{01} (E1_1)     & 0                            & \\dots  & 0
                             & \\Gamma_{12} (E1_1,E2_1)    & \\dots  & \\Gamma_{12} (E1_1,E2_{n2})        \\\\
                         \\vdots                 & \\vdots                      & \\ddots & \\vdots
                             & \\vdots                     & \\ddots & \\vdots                            \\\\
                         \\Gamma_{01} (E1_{n1})  & 0                            & \\dots  & 0
                             & \\Gamma_{12} (E1_{n1},E2_1) & \\dots  & \\Gamma_{12} ( E1_{n1} , E2_{n2} ) \\\\
                         0                       & \\Gamma_{21} (E2_1,E1_1)     & \\dots  & \\Gamma_{21} (E2_1,E1_{n1})
                             & 0                           & \\dots  & 0                                  \\\\
                         \\vdots                 & \\vdots                      & \\ddots & \\vdots
                             & \\vdots                     & \\ddots & \\vdots                            \\\\
                         0                       & \\Gamma_{21} (E2_{n2},E1_1) & \\dots  & \\Gamma_{21} (E2_{n2},E1_{n1})
                             & 0                           & \\dots  & 0                                  \\\\
                         
                     \\end{pmatrix}
                     ,


                
                 where :math:`\\Gamma_ab(E1_i)`  is the value for the transition rate from either an one-body state to the vacuum  :math:`(ab=10)`  or the other way around :math:`(ab=01)`.
                 :math:`E1_i` is the energy of the one-body state (which is the :math:`i`-th value in the E1-array). 

                 :math:`\\Gamma_ab(Ea_i,Eb_j)`,  analogously, is the value for a transition from either an one-body state to a two-body state :math:`(ab=12)`  or the other way around :math:`(ab=21)`.
                 :math:`Ea_i` is the energy of the initial state and :math:`Eb_i` is the energy of the final state (again, read out of either the E1 -or E2-array).

                 :math:`\\Gamma_L`  and :math:`\\Gamma_R` have the same format, that is shown above.
                 But, of course, the values for each transition differs, due to the different chemical potentials in source and drain.

                 **Example** *(we refer to the state with the energy-eigenvalue* :code:`E1[n-1]`  *or*  :code:`E2[n-1]`  *as the n-th one-body or two-body state)*:
                 
                 if the dot is in the 5th one-body state and one wants to have the rate of the transition in which an electron tunnels through the left barrier and puts the dot into the 2nd two-body state,
                 one must read out the matrix element:
                 
                 :math:`(\\Gamma_L)_{1+5,1+n1+2}`

                 

            
        :rtype: list with two numpy matrices: :math:`(\\Gamma_L` , :math:`\\Gamma_R)`
        
       
        """
        NEcut = len(E1) #we determine the number of single-particle states that we use
        VG=np.diag(pot)
        E= int(0.5*np.size(VG))
        V = VG[0:E] #since the potential of both barriers is symmetric and we only tunnel through one barrier. Therefore we only use one half of the potential.
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
            """This fermi-function tells us with which likelyhood a state with an E is occupied on the lead.
            E(float): energy difference between the initial and the final state that the tunneling electron has to carry.
            mu(float): chemical potential of either drain(muR) or source(muL).
            T(float): temperature.
            """
            if (E-mu)/T > 600:
                    f=0
				
            else:
                    f=1/(math.exp((E-mu)/(kB*T) )+1)
            return(f)
          

	#This function is called by the Gamma_ij-equations and includes the transmission-coefficient for each tunnelling-event
        #and the density of state function of the source and drain. 
        def Gamma(Ea,Eb,V):
             """:math:`\\Gamma` includes the transmission coefficient and DOS: :math:`\Gamma = | t |^2 * DOS`

             Ea(float): energy of initial state
             Eb(float): energy of final state
             V(np.array): barrier potential
             """
             #print(Ea)
             #print(V)
             return (np.absolute(TCalc.calculate_transmission(Ea,V,dx))**2*Density.calculate_DensityofStates(np.absolute(Ea-Eb)))
                     
        #These next four functions are used to calculate the transition rates.Each function for a different kind of transition:
        #We distinguish between transitions, in which the number of electrons on the dot changes from one to two(Gamma_12) and reverse(Gamma_21).
        #And between transitions in which the number of electrons on the dot change from zero to one(Gamma_01) and reverse(Gamma_10).

        def Gamma_12(Ea,Eb,mu,T):
            """Calculates the rate of a transition from a one body state to a two body state.

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
            """Calculates the transition rate from the vacuum state to a one-body state.

            Eb(float): energy of final state
            mu(float): chemical potential of either drain(muR) or source(muL)
            T(float): temperature
            """
            return(Gamma(E0,Eb,V)*fermi((Eb-E0),mu,T))

        def Gamma_21(Ea,Eb,mu,T):
            """Calculates the rate of a transition from a two body state to a one body state.

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
            """Calculates the rate of a transition from a one body state to the vacuum state.

            Ea(float): energy of initial state            
            mu(float): chemical potential of either drain(muR) or source(muL)
            T(float): temperature
            """
            return(Gamma(Ea,E0,V)*(1-fermi((Ea-E0),mu,T)))

        #creating the output matrices that later contain all the transition rates through either
        #the left or the right barrier
        Gamma_R=np.zeros((1+np.size(E1)+np.size(E2),1+np.size(E1)+np.size(E2)))
        Gamma_L=np.zeros((1+np.size(E1)+np.size(E2),1+np.size(E1)+np.size(E2)))

        #using a loop to fill the output matrices with transition rates.
        i_=0
        for i in E1:
                j_=0
                for j in E2:
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
        
