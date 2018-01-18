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
     
    def doCalculation(self, E1, E2, muL, muR, T, V, C, TCalc, Density):
      

        kB=0.08629 #Boltzmannkonstante in meV/K

        def fermi(E,mu,T):
            f=1/(math.exp((E-mu)/(kB*T) )+1)
            return(f)
          
	
        def Gamma(Ea,Eb,V):
             return (np.absolute(TCalc.calculate_transmission(np.absolute(Eb-Ea),V))**2*Density.calculate_DensityofStates(np.absolute(Ea-Eb)))

 
        NEcut= np.size(E2)
          
        #Um Tunnelraten zu berechnen, die durch die linke Tunnelbarriere gehen, muss als Parameter mu das chemische Potential der linken Tunnelbariere eingesetzt werden. Umgekehrt symmetrisch fuer die rechte Tunnelbariere

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



#        Gamma_12L=[[0 for i in E2]for i in E1]
#        Gamma_12R=[[0 for i in E2]for i in E1]
#        Gamma_21L=[[0 for i in E2]for i in E1]
#        Gamma_21R=[[0 for i in E2]for i in E1]
#        for i in E1:
#            for j in E2:
#               Gamma_12L[np.where(E1==i)[0][0]][np.where(E2==j)[0][0]]=Gamma_12(i,j,muL,T)
#               Gamma_12R[np.where(E1==i)[0][0]][np.where(E2==j)[0][0]]=Gamma_12(i,j,muR,T)
#               Gamma_21L[np.where(E1==i)[0][0]][np.where(E2==j)[0][0]]=Gamma_21(j,i,muL,T)
#               Gamma_21R[np.where(E1==i)[0][0]][np.where(E2==j)[0][0]]=Gamma_21(j,i,muR,T)
#
#        Gamma_01L=[0 for i in E1]
#        Gamma_01R=[0 for i in E1]
#        Gamma_10L=[0 for i in E1]
#        Gamma_10R=[0 for i in E1]
#        for i in E1:
#            Gamma_01L[np.where(E1==i)[0][0]]=Gamma_01(i,muL,T)
#            Gamma_01R[np.where(E1==i)[0][0]]=Gamma_01(i,muR,T)
#            Gamma_10L[np.where(E1==i)[0][0]]=Gamma_10(i,muL,T)
#            Gamma_10R[np.where(E1==i)[0][0]]=Gamma_10(i,muR,T)

#        Gamma_00=np.array([0])
#        Gamma_11=np.zeros((np.size(E1), np.size(E1)))
#        Gamma_22=np.zeros((np.size(E2), np.size(E2)))
#        Gamma_02=np.zeros((1,np.size(E2)))
#        Gamma_20=np.zeros((np.size(E2),1))
#        Gamma_10L = np.transpose(Gamma_10L)
#        Gamma_21L = np.transpose(Gamma_21L)
#        Gamma_L = np.array([[Gamma_00 , Gamma_01L, Gamma_02 ], [Gamma_10L, Gamma_11 , Gamma_12L], [Gamma_20 , Gamma_21L, Gamma_22 ]])
#        Gamma_10R = np.transpose(Gamma_10R)
#        Gamma_21R = np.transpose(Gamma_21R)
#        Gamma_R = np.array([[Gamma_00 , Gamma_01R, Gamma_02 ],[Gamma_10R, Gamma_11 , Gamma_12R],[Gamma_20 , Gamma_21R, Gamma_22 ]])

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
                Gamma_L[0][i_+1]=Gamma_01(i,muL,T)
                Gamma_R[0][i_+1]=Gamma_01(i,muR,T)
                Gamma_L[i_+1][0]=Gamma_10(i,muL,T)
                Gamma_R[i_+1][0]=Gamma_10(i,muR,T)
                i_=1+i_
        print(Gamma_L)
        print(Gamma_R)
        return(Gamma_L,Gamma_R)
       
















































































































































































































































