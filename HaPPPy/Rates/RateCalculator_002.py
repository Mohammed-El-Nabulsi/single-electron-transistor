#alle folgende Einheiten werden benutzt: [Energie]=meV, [Länge]=nm, [Zeit]=ps

import math
import numpy as np

#Input
NEcut=100          #Anzahl der energiezustände bis schnitt
muL=0
muR=0
Ebeta=1
Ealpha=0.5
summe=0
T=400

#benötigte Konstanten 

kB=1 #Boltzmannkonstante in eV



#provisorische potential
V=5


#provisorisches t gruppe3
def t(E,V):
    return(float(0.8))

#provisorische Energieeigenwerte einteilchen gruppe1
E1=np.arange(0,NEcut,1)
 
#provisorische Energieeigenwerte zweiteilchen gruppe2
E2=np.arange(0,NEcut,1)

#provisorische koeffizient zweiteilchen matrix in matrix äußere matrix besitzt mehrteilchenbasis und innere matrix produktbasis
C=[[(np.arange(0,NEcut,1)/500000) for q in range(NEcut)] for x in range(NEcut)]


#provisorisches D
def D(E):
    return(0.7)



#definiere Fermifunktion

def fermi(E,mu,T):
    f=1/(math.exp((E-mu)/(kB*T) )+1)
    return(f)

#Energie-Kronecker delta in Ratengleich für schnellere Schreibweise E1=EzielzustandDot
    #E2=Eleadszustand E3=EStartzustandDot


def Gamma(Ea,Eb,V):
    return (np.absolute(t(np.absolute(Eb-Ea),V))**2*D(np.absolute(Ea-Eb)))

  
#Um Tunnelraten zu berechnen, die durch die linke Tunnerlbariere gehen, muss als Parameter mu das chemische Potential der linken Tunnelbariere eingesetzt werden. Umgekehrt symmetrisch für die rechte Tunnelbariere

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

print(Gamma_12L)











































































































































































































































































