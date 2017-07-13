#!/usr/bin/python
import numpy as np
import m_init as init
import m_propagator as propa

## Time-propagation of Wavefunction following Schrodinger eqn-----------
## using Crank-Nicolson method with Gauss-Seidel method
## 170614 Modifying the normalization and propagation part for Schrodin-
## ger eqn. (considering that wavefunction is complex)
## 170617 Modularization of the code
## (Dongjae Shin, EEWS, KAIST)

## Getting input parameters--------------------------------------------
PHI, V,\
NSTEP, LEN, H, DT, LEN_TIME, NITER,\
K0, X0, SIGMA, W, V0, X0_V,\
FREQ = init.get_input()

## Initialization of wavefunction and external potentiali--------------
PHI  = init.set_wavefunction(PHI, NSTEP, K0, X0, SIGMA, H, DT)
V    = init.set_ext_potential(V, NSTEP, W, V0, X0_V, H)

## Print initial wavefunction(prob. density) and potential (stdout)----
print '### Initial wavefunction:\n'
print 'x    Re[psi]   Imag[psi]:\n'
for k in range(0, NSTEP+1):
	print '{0:10.5f}{1:10.5f}{2:10.5f}\n'.\
	format(k*H, PHI[k].real, PHI[k].imag)

print 'x    probability density:\n'
for k in range(0, NSTEP+1):
	print '{0:10.5f}{1:10.5f}\n'.\
	format(k*H, np.absolute(PHI[k])**2)

print '### External potential:\n'
print 'x    V(x):\n'
for k in range(0, NSTEP+1):
	print '{0:10.5f}{1:10.5f}\n'.\
	format(k*H, V[k].real)

## Print initial wavefunction (output file)----------------------------
f = open("./output/{0:d}.WF".\
	format(0), 'w')
for k in range(0,NSTEP+1):
	f.write("{0:10.5f}{1:10.5f}\n".\
	format(k*H, np.absolute(PHI[k])**2))
f.close()

## Print external portential ------------------------------------------
f = open("./output/V.PT", 'w')
for k in range(0,NSTEP+1):
	f.write("{0:10.5f} {1:10.5f}\n".\
	format(k*H, V[k].real))
f.close()


## Print wavefunctions w/r/t time ------------------------------------
print '###Time-evolution start!'

for i in range(1, NITER+1):
	PHI = propa.update_CN(i, PHI, V, NSTEP, DT, H)
	if i%10 == 0:
		print 'iteration = {0:10d}, time = {1:10.5f} \n'.\
		format(i, i*DT)
		f = open("./output/{0:d}.WF".\
			format(i), 'w')
		for k in range(0,NSTEP+1):
			f.write("{0:10.5f}{1:10.5f}\n".\
			format(k*H, np.absolute(PHI[k])**2))
		f.close()

print '###TIME-evolution finished!'
