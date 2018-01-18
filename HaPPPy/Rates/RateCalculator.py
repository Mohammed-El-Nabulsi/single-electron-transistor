#alle folgende Einheiten werden benutzt: [Energie]=meV, [Laenge]=nm, [Zeit]=ps

import math
import numpy as np
from HaPPPy.OneBody import OneBodySolver
from HaPPPy.TwoBody import TwoBodySolver
from HaPPPy.Transmission import TransmissionCalculator
from HaPPPy.MasterEquation import MasterEquationSolver

#E1=Energie_Einteilchen # OBSolver.doCalculation()?
#V=Vinput[:half]
#E2 = Energie_Zweiteilchen 
#C= Koeffizienten_Zweiteilchen #TBSolver.doCalculation()
#muL, muR,T, Vinput aus main importieren 
#Transimissionsmatrix wird in def Gamma importiert
class RateCalculator:
    
    def __init__(self):
        """ The constructor.
        """

        print("Hello from the RateCalculator, what can I do for you?")
    
   

    #benoetigte Konstanten 

    kB=1 #Boltzmannkonstante in eV
    

    #definiere Fermifunktion
    def doCalculation(self, E1, E2, muL, muR, T, V, C):
        kB=1
        def fermi(E,mu,T):
            f=1/(math.exp((E-mu)/(kB*T) )+1)
            return(f)

        def D(A):
            if A == 0:
                return(1)
            else:
                return(0)
        def Gamma(Ea,Eb,V):
            return (np.absolute(t(np.absolute(Eb-Ea),V))**2*D(np.absolute(Ea-Eb)))

        def t(E, V):
            return(0.1)
        NEcut= np.size(E2)
          
        #Um Tunnelraten zu berechnen, die durch die linke Tunnerlbariere gehen, muss als Parameter mu das chemische Potential der linken Tunnelbariere eingesetzt werden. Umgekehrt symmetrisch fuer die rechte Tunnelbariere

        def Gamma_12(Ea,Eb,mu,T): 
            summe=0
            #Skalarproduktsumme (10.29)
            j=0
            Cb=C[np.where(E2==Eb)[0][0]]
            while j< NEcut:
                summe=Cb[np.where(E1==Ea)[0][0]][j]+summe
                j=j+1
            return(Gamma(Ea,Eb,V)*(np.absolute(summe))**2*fermi(np.absolute(Eb-Ea),mu,T))


        def Gamma_01(Eb,mu,T): 
            return(Gamma(0,Eb,V)*fermi(Eb,mu,T))

        def Gamma_21(Ea,Eb,mu,T): 
            summe=0
            nu=0
            Ca=C[np.where(E2==Ea)[0][0]]
            while nu < NEcut:
                 summe=summe+Ca[np.where(E1==Eb)[0][0]][nu]
                 nu=nu+1
            return(Gamma(Ea,Eb,V)*(np.absolute(summe))**2*(1-fermi(np.absolute(Eb-Ea),mu,T)))

        def Gamma_10(Ea,mu,T): #anfangszustand a einteilchen
            return(Gamma(Ea,0,V)*(1-fermi(Ea,mu,T)))



        Gamma_12L=[[0 for i in E2]for i in E1]
        Gamma_12R=[[0 for i in E2]for i in E1]
        Gamma_21L=[[0 for i in E2]for i in E1]
        Gamma_21R=[[0 for i in E2]for i in E1]
        for i in E1:
            for j in E2:
               Gamma_12L[np.where(E1==i)[0][0]][np.where(E2==j)[0][0]]=Gamma_12(i,j,muL,T)
               Gamma_12R[np.where(E1==i)[0][0]][np.where(E2==j)[0][0]]=Gamma_12(i,j,muR,T)
               Gamma_21L[np.where(E1==i)[0][0]][np.where(E2==j)[0][0]]=Gamma_21(j,i,muL,T)
               Gamma_21R[np.where(E1==i)[0][0]][np.where(E2==j)[0][0]]=Gamma_21(j,i,muR,T)

        Gamma_01L=[0 for i in E1]
        Gamma_01R=[0 for i in E1]
        Gamma_10L=[0 for i in E1]
        Gamma_10R=[0 for i in E1]
        for i in E1:
            Gamma_01L[np.where(E1==i)[0][0]]=Gamma_01(i,muL,T)
            Gamma_01R[np.where(E1==i)[0][0]]=Gamma_01(i,muR,T)
            Gamma_10L[np.where(E1==i)[0][0]]=Gamma_10(i,muL,T)
            Gamma_10R[np.where(E1==i)[0][0]]=Gamma_10(i,muR,T)
        print(Gamma_01L,Gamma_01R,Gamma_10L,Gamma_10R,Gamma_12L,Gamma_12R,Gamma_21L,Gamma_21R)
        print('reihenfolge 01L,01R,10L,10R,12L,12R,21L,21R')
        return(Gamma_01L,Gamma_01R,Gamma_10L,Gamma_10R,Gamma_12L,Gamma_12R,Gamma_21L,Gamma_21R)
        





































































































































































































































































