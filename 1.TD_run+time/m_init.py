import numpy as np
import m_funcs as func

def get_input():
	## Getting spatial parameters
	NSTEP 	= input('No. of grids:\n')
	LEN 	= float(input('Spatial length:\n'))
	H 		= LEN/NSTEP

	## Getting temporal parameters
	DT 		= float(input('Time interval:\n'))
	if DT == 0.:
		print 'error) time step should be positive real.'
		exit()
	LEN_TIME = float(input('Total time:\n'))
	NITER 	 = int(LEN_TIME/DT)

	DTH 	= DT/H**2

	## Getting parameters for initial wavefunction
	K_factor= float(input('K0 factor:\n'))
	K0 		= K_factor*np.sqrt(1./DT)
	X0 		= float(input('X0:\n'))
	SIGMA 	= float(input('Sigma:\n'))

	## Getting parameters for external potential
	W 		= float(input('W:\n'))
	V0 		= float(input('V0:\n'))
	X0_V 	= float(input('X0_V:\n'))

	## Getting other parameters
	FREQ	= int(input('FREQ:\n')) # output frequency

	PHI 	= np.ones(NSTEP+1, dtype=complex) # array for wavefunction 
	V 	= np.ones(NSTEP+1, dtype=complex)


	return PHI, V,\
	       NSTEP, LEN, H, DT, LEN_TIME, NITER,\
	       K0, X0, SIGMA, W, V0, X0_V,\
	       FREQ

def set_wavefunction(PHI, NSTEP, K0, X0, SIGMA, H, DT):
	PHI[0]		= 0. + 0j
	PHI[NSTEP] 	= 0. + 0j
	NORM 		= 0.

	## Initial shape of wavefunction
	for IX in range(1,NSTEP): # loop from 1 to NSTEP-1
		PHI[IX] = func.GAUSS_packet(IX*H, K0, X0, SIGMA, DT)
		# Gussian wavepacket

	## Normalization of wavefunction
	for IX in range(0,NSTEP): # integration of probability (trapezoidal)
		NORM += np.absolute(PHI[IX])**2 + np.absolute(PHI[IX+1])**2
	NORM *= H/2.0
	SQRTNORM = np.sqrt(NORM)
	for IX in range(1,NSTEP): # normalize the function
		PHI[IX] /= SQRTNORM

	return PHI

def set_ext_potential(V, NSTEP, W, V0, X0_V, H):
	## External potential
	for IX in range(0,NSTEP+1):
		V[IX] = func.POTENTIAL_SQW(IX*H, X0_V, W, V0)
	
	return V
