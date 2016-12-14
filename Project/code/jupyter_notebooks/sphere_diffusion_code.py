################################################################################################
# Commented code can be found in the sphere_diffusion.ipynb Jupyter notebook. This has been
# copied and pasted so I can import it elsewhere.
################################################################################################
import math
from math import exp
from math import cosh
from math import sinh
from math import sqrt
import numpy as np
from numpy.polynomial.legendre import leggauss
import matplotlib.pyplot as plt
################################################################################################
# REFLECTED SPHERE W/ UNIFORM SOURCE
def det_constants(R1, R2, xs_a1, xs_a2, D1, D2, s):
	L1 = sqrt(D1/xs_a1)
	L2 = sqrt(D2/xs_a2)
	# extrapolation distance
	a = 2.0/3.0/xs_a2
	A = exp(2*(R2+a)/L2)
	B = 2*sinh(R1/L1)/R1
	C = (R1*cosh(R1/L1)-L1*sinh(R1/L1))/(L1*R1**2)
	D = (A*exp(-R1/L2)+exp(R1/L2))/L2/R1 - (exp(R1/L2)-A*exp(-R1/L2))/R1**2
	const3 = (2*D1*s/xs_a1)/((2*D1/R1)*(exp(R1/L2)-A*exp(-R1/L2))-D*B/C)
	const1 = ((const3/R1)*(exp(R1/L2)-A*exp(-R1/L2))-(s/xs_a1))/B
	return(const1, const3)
#----------------------------------------------------------------------------------------------#
def det_reflected_flux(r, R1, R2, xs_tot1, xs_s1, xs_tot2, xs_s2, s, mu):
	# mu is assumed to be the same in both the core and reflector
	# Determine properties for core
	xs_a1 = xs_tot1-xs_s1
	xs_tr1 = xs_a1 - (1-mu)*xs_s1
	D1 = 1.0/3.0/xs_tr1
	L1 = sqrt(D1/xs_a1)
	# Determine properties for reflector
	xs_a2 = xs_tot2-xs_s2
	xs_tr2 = xs_a2 - (1-mu)*xs_s2
	D2 = 1.0/3.0/xs_tr2
	L2 = sqrt(D2/xs_a2)
	# calculate equation constants
	C1, C3 = det_constants(R1, R2, xs_a1, xs_a2, D1, D2, s)
	# extrapolation distance
	a = 2.0/3.0/xs_a2
	if(r<R1):
		phi = (2*C1/r)*sinh(r/L1)+s/xs_a1
	#elif(r>=R1 and r<R2):
	else:
		phi = C3*(exp(r/L2)-exp(2*(R2+a)/L2)*exp(-r/L2))/r
	return(phi)
################################################################################################
# BARE SPHERE W/ UNIFORM SOURCE
def det_bare_flux(r, R, source, xs_tot, xs_s, mu):
	xs_a = xs_tot-xs_s
	xs_tr = xs_a - (1-mu)*xs_s
	D = 1.0/3.0/xs_tr
	L = sqrt(D/xs_a)
	a = 2.0/3.0/xs_a
	c1 = -source*(R+a)/xs_a/2/sinh((R+a)/L)
	phi = (2*c1/r)*sinh(r/L)+source/xs_a
	return(phi)
################################################################################################
# BARE SPHERE W/ POINT SOURCE AT ORIGIN
def diffusion_point_source(r, R, total_xsec, scattering_xsec, source, mu):
	#abs_xsec = total_xsec -scattering_xsec
	#mu = 2.0/3.0/A
	D = 1.0/3.0/(total_xsec-(1-mu)*scattering_xsec)
	L = sqrt(D/total_xsec)
	flux = source*sinh((R-r)/L)/4/math.pi/D/sinh(R/L)
	return(flux)
################################################################################################





