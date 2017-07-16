import numpy as np
import m_linalg as linalg
import scipy.sparse.linalg as sp

## This module contains several time propagator functions ---------------------
## 1. Explicit method
## 2. Crank-Nicolson method
## Input:
## Output: 
## 170617 Dongjae Shin, EEWS, KAIST
## ----------------------------------------------------------------------------

## Time evolution of a function using explicit method (not used)------
def update_explicit(i, PHI, NSTEP, DT, DTH):
## ---------------------------------------------------------------------
	# Loop over time steps using explicit method
	POLD = 0.
	for IX in range(1,NSTEP):
		PNEW = PHI[IX]+DTH*(POLD+PHI[IX+1]-2*PHI[IX])
		POLD = PHI[IX]	# roll POLD value
		PHI[IX] = PNEW	# store new value
	return PHI

## Time evolution of wavefunction using Crank-Nicolson method ----------
def update_CN(i, PHI, V, NSTEP, DT, H, eps, solver=1):
## ---------------------------------------------------------------------
	## Loop over time steps
	## Initialization of imtermediate function(CHI), [A], {b}
	A 			= np.zeros((NSTEP-1, NSTEP-1), dtype=complex)
	b 			= np.zeros(NSTEP-1, dtype=complex)
	CHI 		= np.zeros(NSTEP-1, dtype=complex)
	CHI_0 		= 0. + 0j # boundary condition
	CHI_NSTEP 	= 0. + 0j

	## Construction of tridiagonal matrix [A] and vector {b}
	for j in range(0,NSTEP-1): # loop from 1 to NSTEP-1
		if j == 0:
			A[j,j] 		= -2. + 2*1j*H**2/DT - H**2*V[j+1]
			A[j,j+1] 	=  1.
			b[j] 		= 4*1j*H**2/DT * PHI[j+1] - CHI_0
		elif j == NSTEP-2:
			A[j,j-1] 	= 1.
			A[j,j] 		= -2. + 2*1j*H**2/DT - H**2*V[j+1]
			b[j] 		= 4*1j*H**2/DT * PHI[j+1] - CHI_NSTEP
		else:
			A[j,j-1] 	= 1.
			A[j,j] 		= -2. + 2*1j*H**2/DT - H**2*V[j+1]
			A[j,j+1] 	= 1.
			b[j] 		= 4*1j*H**2/DT * PHI[j+1]
	if   solver == 1:
		CHI = sp.spsolve(A,b) 
	elif solver == 2:
		CHI = linalg.gaussSeidel(A,b,eps) 
		# too slow!
	else:
		print "update_CN: inappropiate argument for solver!"
		exit()

	## Update the function using intermediate function, CHI
	for IX in range(1,NSTEP):
		PHI[IX] = CHI[IX-1] - PHI[IX]

    ## Normalization on each time step
	NORM 	 = 0.
	SQRTNORM = 0.
	for IX in range(0,NSTEP): # integration of probability (trapezoidal)
		NORM += np.absolute(PHI[IX])**2 + np.absolute(PHI[IX+1])**2
	NORM *= H/2.0
	SQRTNORM = np.sqrt(NORM)
	for IX in range(1,NSTEP): # normalize the function
		PHI[IX] /= SQRTNORM

	return PHI

