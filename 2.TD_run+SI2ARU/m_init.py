import numpy as np
import m_funcs as func

ARU_2_FS = 4.8377687e-2

def get_input():

	K_E0 	= float(input('K_E0:\n')) # kinetic energy
	X0 	= float(input('X0:\n')) # center of GWP
	SIGMA 	= float(input('SIGMA:\n'))
	V_xi 	= float(input('V_xi:\n')) # initial position of V0
	V_xf 	= float(input('V_xf:\n')) # final position of V0
	V0 	= float(input('V0:\n')) # potential height

	## Getting spatial parameters
	LEN 	= float(input('Spatial length:\n'))
	NSTEP 	= input('No. of grids:\n')
	H 		= LEN/NSTEP # = dX

	## Getting temporal parameters
	LEN_TIME 	= float(input('Total time:\n')) #Total time
	DT 		= float(input('Time interval:\n'))
	if DT == 0.:
		print 'error) time step should be positive real.'
		exit()
	NITER 	= int(LEN_TIME/DT) # the no. of time iteration

	# SI unit into ARU
	# Atomic Rydberg units (ARU) defined by:
	# h_bar = 2*m_e = e^2/2 = 1, a_o = 4*pi*eps_o = 1
	# spatial range is initially in ARU

	VELOCITY = np.sqrt(4 * K_E0)
	DTH 	 = DT/H**2

	## Getting other parameters
	FREQ	= int(input('FREQ:\n')) # output frequency
	eps	= float(input('eps:\n')) # tolerance for convergence

	PHI 	= np.ones(NSTEP+1, dtype=complex) # array for wavefunction 
	V 	= np.ones(NSTEP+1, dtype=complex) # array for potential
	X	= np.ones(NSTEP+1) # array for spatial domain
	for IX in range(0,NSTEP+1):
		X[IX] = IX * H

	return K_E0, X0, SIGMA, V_xi, V_xf, V0,\
	       LEN, NSTEP, H, LEN_TIME, DT, NITER, DTH, VELOCITY,\
	       FREQ, eps,\
	       PHI, V, X

def set_wavefunction(X, PHI, NSTEP, SIGMA, X0, H, VELOCITY):
	PHI[0]		= 0. + 0j
	PHI[NSTEP] 	= 0. + 0j
	NORM 		= 0.

	## Initial shape of wavefunction
	for IX in range(1,NSTEP): # loop from 1 to NSTEP-1
		PHI[IX] = func.GWP(X[IX], SIGMA, X0, VELOCITY)
		# Gussian wavepacket

	## Normalization of wavefunction
	for IX in range(0,NSTEP): # integration of probability (trapezoidal)
		NORM += np.absolute(PHI[IX])**2 + np.absolute(PHI[IX+1])**2
	NORM *= H/2.0
	SQRTNORM = np.sqrt(NORM)
	for IX in range(1,NSTEP): # normalize the function
		PHI[IX] /= SQRTNORM

	return PHI

def set_ext_potential(X, V, NSTEP, V_xi, V_xf, V0):
	## External potential
	for IX in range(0,NSTEP+1):
		V[IX] = func.POTENTIAL_SQW_2(X[IX], V_xi, V_xf, V0)
	
	return V
