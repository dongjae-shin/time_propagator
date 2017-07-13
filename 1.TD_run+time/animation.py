#!/usr/bin/python
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import m_init as init

## Getting input parameters--------------------------------------------
PHI, V,\
NSTEP, LEN, H, DT, LEN_TIME, NITER,\
K0, X0, SIGMA, W, V0, X0_V,\
FREQ = init.get_input()
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

## Set up the initial figure
## 170712 time_text was added.
fig 	= plt.figure()
ax 	= plt.axes(xlim=(10,40), ylim=(0,1))
line, 	= ax.plot([], [], lw=2)
time_text = ax.text(0.05, 0.95,'',horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)
plt.ylabel("|psi(x)|^2")
plt.xlabel("x")

## Initialization of wavefunction (for animation) ----------------------
def init():
	x, V = np.loadtxt("V.PT", unpack=True)
	line.set_data(x, V)
	time_text.set_text('')
	return line, time_text,
## Time evolution of wavefunction using output wavefunctions -----------
def update(i):
## ---------------------------------------------------------------------
	x, psi = np.loadtxt("{0}".\
		format(filenames[i]), unpack=True)
	line.set_data(x, psi)
	time_text.set_text('time = {0:10.5f}'.format(i*DT))
	return line, time_text,

length 	= int(len(filenames)) # because of FREQ, length < NITER
ani 	= animation.FuncAnimation(fig, update, init_func=init,\
     	 frames=length, interval=10, blit=True, repeat=True)
#ani.save('Crank-Nicolson.gif', dpi=80, writer='imagemagick')

plt.show()
