import numpy as np
# ----------------------------------------------------------------------
# 170710 Newer form of Gaussian wavepacket was added.
# ----------------------------------------------------------------------
# Defining analytical solution------------------------------------------
def GAUSS(X, T):
	value = np.exp(-20.*(X-0.5)**2/(1.+80*T))/np.sqrt(1+80*T)
	return value
def EXACT(X, T):
	value = GAUSS(X, T)-GAUSS(X-1., T)-GAUSS(X+1., T)
	return value
def GWP(X, SIGMA, X0, K0):
	value = (2*np.pi*SIGMA**2)**(-0.25)*np.exp(1j*K0*(X-X0))\
		*np.exp(-(X-X0)**2/(2*SIGMA**2))
	return value
# Description of oscillations of the initial wavepacket-----------------
def GAUSS_packet(X,K0,X0,SIGMA,DT):
# ----------------------------------------------------------------------
	if DT >= 1/K0**2: # -sqrt(1/DT) < K0 < sqrt(1/DT)
		print 'error) DT < 1 / K0**2 should be satisfied.'
		print 'error) K0 should be b/w -{0} and {1}'\
			  .format(np.sqrt(1./DT), np.sqrt(1./DT))
		exit()

	value = np.exp(1j*K0*X) * np.exp(-(X-X0)**2/(2*SIGMA**2))
	return value

# Potential function----------------------------------------------------
def POTENTIAL_SQW(X,X0,W,V0):
# ----------------------------------------------------------------------
	if X < X0-W/2 or X >= X0+W/2:
		value = 0. 
	elif X >= X0-W/2 and X < X0+W/2:
		value = V0
	return value
def POTENTIAL_SQW_2(X, V_xi, V_xf, V0):
	if X < V_xi or X >= V_xf:
		value = 0.
	elif X >= V_xi and X < V_xf:
		value = V0
	return value


