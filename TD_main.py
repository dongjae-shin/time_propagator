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
K_E0, X0, SIGMA, V_xi, V_xf, V0,\
LEN, NSTEP, H, LEN_TIME, DT, NITER, DTH, K0,\
FREQ, eps,\
PHI, V, X = init.get_input()

## Initialization of wavefunction and external potentiali--------------
PHI  = init.set_wavefunction(X, PHI, NSTEP, SIGMA, X0, H, K0)
V    = init.set_ext_potential(X, V, NSTEP, V_xi, V_xf, V0)

## Print initial wavefunction(prob. density) and potential (stdout)----
print '### Initial wavefunction:\n'
print 'x [Bohr]   Re[psi]   Imag[psi]:\n'
for k in range(0, NSTEP+1):
	print '{0:12.7f}{1:12.7f}{2:12.7f}\n'.\
	format(X[k], PHI[k].real, PHI[k].imag)

print 'x [Bohr]   probability density:\n'
for k in range(0, NSTEP+1):
	print '{0:12.7f}{1:12.7f}\n'.\
	format(X[k], np.absolute(PHI[k])**2)

print '### External potential:\n'
print 'x [Bohr]   V(x) [Ryd]:\n'
for k in range(0, NSTEP+1):
	print '{0:12.7f}{1:12.7f}\n'.\
	format(X[k], V[k].real)

## Print initial wavefunction (output file)----------------------------
f = open("./output/{0:d}.WF".\
	format(0), 'w')
for k in range(0,NSTEP+1):
	f.write("{0:12.7f}{1:12.7f}\n".\
	format(X[k], np.absolute(PHI[k])**2))
f.close()

## Print external portential ------------------------------------------
f = open("./output/V.PT", 'w')
for k in range(0,NSTEP+1):
	f.write("{0:12.7f} {1:12.7f}\n".\
	format(X[k], V[k].real))
f.close()


## Print iwavefunctions w/r/t time ------------------------------------
print '###Time-evolution start!'

for i in range(1, NITER+1):
	PHI = propa.update_CN(i, PHI, V, NSTEP, DT, H, eps)
	if i%FREQ == 0:
		print 'iteration = {0:10d}, time = {1:10.5f} fs\n'.\
		format(i, i*DT*init.ARU_2_FS)
		f = open("./output/{0:d}.WF".\
			format(i), 'w')
		for k in range(0,NSTEP+1):
			f.write("{0:12.7f}{1:12.7f}\n".\
			format(X[k], np.absolute(PHI[k])**2))
		f.close()

print '###TIME-evolution finished!'
