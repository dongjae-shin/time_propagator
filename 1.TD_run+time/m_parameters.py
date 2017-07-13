import numpy as np

## This module contains global variables --------------------------------------
##
## 170617 Dongjae Shin, EEWS, KAIST
## ----------------------------------------------------------------------------

## Global variables (default values) ------------------------------------------
NSTEP 	= 100. # no. of grid
LEN 	= 50. # length of grid
H 	= LEN/NSTEP

DT 		= 0.0008 # time interval
LEN_TIME 	= 10. # total simulation time
NITER 		= int(LEN_TIME/DT)
T 		= 0.0 # current time

DTH 		= DT/H**2

## Parameters -----------------------------------------------------------------
NORM 	= 0.0 # the integration value of probability density

K0 	= 0.0 # related to momentum of wavepacket (may be transferred to...)
X0 	= 0.0 # the center of wavepacket (may be transferred to...)
SIGMA 	= 0.0 # width of wavepacket (may be transferred to...)

W 	= 0.0 # the width of external potential  (may be transferred to...)
V0 	= 0.0 # height/depth of external potential (may be transferred to...)
X0_V 	= 0.0 # center of external potential (may be transferred to...)


