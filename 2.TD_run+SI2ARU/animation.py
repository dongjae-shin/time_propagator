#!/usr/bin/python
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import m_init as init

## Getting input parameters--------------------------------------------
K_E0, X0, SIGMA, V_xi, V_xf, V0,\
LEN, NSTEP, H, LEN_TIME, DT, NITER, DTH, VELOCITY,\
FREQ, eps,\
PHI, V, X = init.get_input()
# PHI and V are not used

## Fill the file names with zeros and sort
width = len(str(NITER)) + 4

os.chdir('./output')
filenames = glob.glob('*.WF')
for name in filenames:
	os.system('mv {0} {1}'.format(name, name.zfill(width)))
for i in range(0,len(filenames)):
	filenames[i] = filenames[i].zfill(width)
filenames.sort()

print filenames

X *= init.ANG_2_BOHR

## Set up the initial figure
fig 	= plt.figure()
ax 	= plt.axes(xlim=(X[0],X[NSTEP]), ylim=(0,0.1))
line, 	= ax.plot([], [], lw=2)
plt.ylabel("|psi(x)|^2")
plt.xlabel("x [Bohr]")

## Initialization of wavefunction (for animation) ----------------------
def init():
	x, V = np.loadtxt("V.PT", unpack=True)
	line.set_data(x, V)
	return line,
## Time evolution of wavefunction using output wavefunctions -----------
def update(i):
## ---------------------------------------------------------------------
	x, psi = np.loadtxt("{0}".\
		format(filenames[i]), unpack=True)
	line.set_data(x, psi)
	return line,

length 	= int(len(filenames)) # because of FREQ, length < NITER
ani 	= animation.FuncAnimation(fig, update, init_func=init,\
     	 frames=length, interval=1, blit=True, repeat=True)
#ani.save('Crank-Nicolson.gif', dpi=80, writer='imagemagick')

plt.show()