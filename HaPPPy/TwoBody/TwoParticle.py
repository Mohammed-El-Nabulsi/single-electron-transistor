import numpy as np 
from .MatrixElement import getMatrixElement, testFunction

## Set up  
n=6 #number of one-particle eigenmodes 
X=100 # number of grid points 
V = [0 for n in range(0,X)]

### Testfunctions for the eigenemodes 
#	createTwoParticleData(n,);

def test(x, n):
    return np.sin(x*n)

def getOneParticleStuff(n,V):
	W=np.zeros((n, X))
	e=np.zeros(n)
	for i in range(n):
		W[i,:] = testFunction(-1.0, 1.0, i+1, X)
		e[i] = float(i+1)**2
		#for j in range(X):
		#	W[i, j]=test(j, i)
	#e=[1, 2, 3, 4, 5, 6]
	return [e,W]


#write the Eigenenergies and Eigenstates to a file
def save(E,S):
	return None

#returns the i-th eigenenergy from the file created by createTwoParticleData. 
#Gives an error if such a file does not exist (i.e. the creation function was not called)
def getEigenenergy(i):
	return 0

#returns the i-th eigenvector from the file created by createTwoParticleData.
#Gives an error if such a file does not exist.
def getEigenvector(i):
	return zeros(X)

# Calculate and store the two-electron eigenfunctions in the potential V from the first n single-electron eigenfunctions. Returns None.
def createTwoParticleData(n,V):
	[SP_EE,SP_EV] = getOneParticleStuff(n,V)
	
	# Basis: 
	#       Singuletts: n - terms: |1>|1>, ...|n>|n>
	#                   1/2 n(n-1) terms: |1>|2>, ....|1>|n>, |2>|3>,..., |2>|n>, ..., |n-1>|n>  
	#       Tripletts: 1/2 n(n-1) terms: |1>|2>, ....|1>|n>, |2>|3>,..., |2>|n>, ..., |n-1>|n>
	I=np.zeros((int(1/2*n*(n-1)),2))
	k=0
	for i in range(n): 
		for j in range(i+1, n):
			I[k, 0]=i 
			I[k, 1]=j
			k=k+1

	# Setup the Matrix in 4 parts: 
	#Symmetric spatial wavefunctions: A (<m,m|V|n,n,>), B(<m,m|V|a,b>), C(<ab|V|cd>)
	# Antisymmetric spatial wavefunctions D (<a,b|V|cd>)

	# Matrix A
	A=np.zeros((n,n))
	for i in range (n):
		for j in range(n):
			A[i,j]=getMatrixElement(SP_EV[i, :], SP_EV[i, :], SP_EV[j, :], SP_EV[j, :])

	# Matrix B und Matrix B transponiert 
	L=int(1/2*n*(n-1))
	B=np.zeros((L, n))

	for i in range(n):
		for j in range(L):
			a=int(I[j, 0])
		   
			b=int(I[j, 1])
			B[j, i]=np.sqrt(2)*getMatrixElement(SP_EV[a, :], SP_EV[b,:],SP_EV[i, :], SP_EV[i, :] )

	# Matrix C und D: 
	C=np.zeros((L, L))
	D=np.zeros((L, L))

	for i in range(L): 
		for j in range(i+1):
			a=int(I[i, 0])
			b=int(I[i, 1])
			c=int(I[j, 0])
			d=int(I[j, 1])
			PartA=getMatrixElement(SP_EV[a, :], SP_EV[b, :],SP_EV[c, :], SP_EV[d, :])+getMatrixElement(SP_EV[b,:], SP_EV[a,:],SP_EV[d,:], SP_EV[c,:])
			PartB=getMatrixElement(SP_EV[a,:], SP_EV[b,:],SP_EV[d,:], SP_EV[c,:])+getMatrixElement(SP_EV[b, :], SP_EV[a,:],SP_EV[c,:], SP_EV[d,:])
			C[i,j]=1/2*(PartA+PartB)
			C[j,i]=C[i,j]
			D[i,j]=1/2*(PartA-PartB)
			D[j,i]=D[i,j]

	# Eigenenergien der Einteilchen kombiniert 

	Energies=np.zeros( n**2)
	for i in range(n): 
		Energies[ i]=2*SP_EE[i]
	for i in range(L): 
		a=int(I[i, 0])
		b=int(I[i,1])
		Energies[ n+i]=SP_EE[a]+SP_EE[b]
		Energies[n+L+i]=Energies[ n+i]

	MatrixEnergies=np.diag(Energies)

	# Fülle die Matrix MatrixAll

	MatrixAll=np.zeros((n**2, n**2))

	MatrixAll[0:n, 0:n]=A
	MatrixAll[n:n+L,0:n]=B 
	MatrixAll[0:n, n:n+L]=np.transpose(B)
	MatrixAll[n:n+L, n:n+L]=C
	MatrixAll[n+L:n**2, n+L:n**2]=D
	MatrixAll=MatrixAll+MatrixEnergies

	# Berechne die Eigenwerte und Eigenvektoren des Problems 
	[Eigenenergies, Eigenvectors]=np.linalg.eig(MatrixAll) # Eigenvalues and eigenvectors 

	#Transform the eigenvectors back to the product state basis

	Eigenvectors_Productbasis = np.zeros((n**2,n**2))
	
	Z = np.zeros((n**2,n**2))
	for i in range(n**2):
		if(i<n):
			Z[i,i]=1
		elif(i<n+L):
			Z[i,i]=1/np.sqrt(2)
			Z[i,i+L]=1/np.sqrt(2)
		else:
			Z[i,i]=(-1)/np.sqrt(2)
			Z[i,i-L]=1/np.sqrt(2)

	Eigenvectors_Productbasis = np.dot(Eigenvectors, Z)
	#Transform  Eigenvectors_Productbasis into n**2 Matrices with the coeeficients 
	Q=np.zeros((n,n, n**2))
	for i in range(n):
		for j in range(n):
			if i==j:
				Q[i,j,:]=Eigenvectors_Productbasis[i, :]
			elif i<j:
				Q[i,j,:]=Eigenvectors_Productbasis[int((2*n-i-1)*i/2)+j-1,:]
			else:
				Q[i,j,:]=Eigenvectors_Productbasis[int((2*n-i-1)*i/2)+j-1+L,:]				
	

	save(Eigenenergies, Eigenvectors_Productbasis)
	print(Eigenenergies)
	#print(Eigenvectors)#_Productbasis)
	return Q	

createTwoParticleData(3,None)
