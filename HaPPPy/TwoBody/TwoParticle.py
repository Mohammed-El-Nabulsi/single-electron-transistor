import numpy as np 
from .MatrixElement import getMatrixElement, testFunction
from .OneParticleLoader import SpectrumData

_SQRT2 = np.sqrt(2)

def createTwoParticleData(opData):
	""" Calculate and return the two-electron eigenvalues and eigenvectors from the single-electron eigenfunctions.
	
	Arguments:
	opData -- the SpectrumData object containing the single-electron spectrum data to use for calculation
	
	Returns [E,Q], where E is an array of (scalar) eigen values and Q is an array of matrices of the shape [i,j,n], where Q[i,j,n] is the coefficient of the nth eigenvector belonging to the |i,j> product basis function. These arrays do not have a special ordering, but Q[:,:,n] is the eigenvector belonging to E[n].
	"""
	SP_EE = opData.energies[:]
	SP_EV = opData.waves[:,:]
	n = opData.m
	X = opData.n
	dx = opData.dx

	# Basis: 
	#       Singuletts: n - terms: |1>|1>, ...|n>|n>
	#                   1/2 n(n-1) terms: |1>|2>, ....|1>|n>, |2>|3>,..., |2>|n>, ..., |n-1>|n>  
	#       Tripletts: 1/2 n(n-1) terms: |1>|2>, ....|1>|n>, |2>|3>,..., |2>|n>, ..., |n-1>|n>
	#
	# The index set I is the ordering of the mixed singlet and triplet terms: I[k] is a pair of indices [i,j] meaning the basis element corresponding to |i,j>.
	# In the case of singlets/triplets, this is the function 1/sqrt(2)(|i,j> +/- |j,i>)
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
			vi = SP_EV[:, i]
			vj = SP_EV[:, j]
			A[i, j] = getMatrixElement(dx, X, vi, vi, vj, vj)

	# Matrix B und Matrix B transponiert 
	L=int(1/2*n*(n-1))
	B=np.zeros((L, n))

	for i in range(n):
		for j in range(L):
			va = SP_EV[:, int(I[j, 0])]
			vb = SP_EV[:, int(I[j, 1])]
			vi = SP_EV[:, i]
			B[j, i] = _SQRT2 * getMatrixElement(dx, X, va, vb, vi, vi)

	# Matrix C und D: 
	C=np.zeros((L, L))
	D=np.zeros((L, L))

	for i in range(L):
		print(i, "of", L)
		for j in range(i+1):
			va = SP_EV[:, int(I[i, 0])]
			vb = SP_EV[:, int(I[i, 1])]
			vc = SP_EV[:, int(I[j, 0])]
			vd = SP_EV[:, int(I[j, 1])]
			PartA = getMatrixElement(dx, X, va, vb, vc, vd) + getMatrixElement(dx, X, vb, va, vd, vc)
			PartB = getMatrixElement(dx, X, va, vb, vd, vc) + getMatrixElement(dx, X, vb, va, vc, vd)
			C[i,j] = 1/2 * (PartA+PartB)
			C[j,i] = C[i,j]
			D[i,j] = 1/2 * (PartA-PartB)
			D[j,i] = D[i,j]

	# Eigenenergien der Einteilchen kombiniert 

	Energies=np.zeros(n**2)
	for i in range(n):
		Energies[ i]=2*SP_EE[i]
	for i in range(L): 
		a=int(I[i, 0])
		b=int(I[i,1])
		Energies[ n+i]=SP_EE[a]+SP_EE[b]
		Energies[n+L+i]=Energies[ n+i]

	MatrixEnergies=np.diag(Energies)

	# Fuelle die Matrix MatrixAll

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
	
	#Z is the Matrix of the basis change from the product basis to the singlet/triplet-basis
	Z = np.zeros((n**2,n**2))
	for i in range(n**2):
		if(i<n):
			Z[i,i]=1
		elif(i<n+L):
			Z[i,i]=1/_SQRT2
			Z[i,i+L]=1/_SQRT2
		else:
			Z[i,i]=(-1)/_SQRT2
			Z[i,i-L]=1/_SQRT2

	#Eigenvectors are the eigenvectors in the singlet/triplet-basis, so the basis transformation Z gives the eigenvectors in the product basis.
	Eigenvectors_Productbasis = np.dot(Eigenvectors, Z)
	
	#Transform  Eigenvectors_Productbasis into n**2 Matrices with the coefficients 
	Q=np.zeros((n,n, n**2))
	for i in range(n):
		for j in range(n):
			if i==j:
				Q[i,j,:]=Eigenvectors_Productbasis[i, :]
			elif i<j:
				Q[i,j,:]=Eigenvectors_Productbasis[int((2*n-i-1)*i/2)+j-1,:]
			else:
				Q[i,j,:]=Eigenvectors_Productbasis[int((2*n-i-1)*i/2)+j-1+L,:]				
	

	return [Eigenenergies, Q]
