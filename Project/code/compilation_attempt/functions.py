import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import legendre
import matplotlib.pyplot as plt
import input
#####################################################################################################
def det_shell_surface_area(r, h):
	# r is the spherical shell center point
	# h is the mesh spacing, which allows us to determine the i+1/2
	# and i-1/2 shell surface areas.
	Aplus = 4*np.pi*(r+(h/2))**2
	Aminus = 4*np.pi*(r-(h/2))**2
	return(Aplus, Aminus)
#####################################################################################################
def det_shell_volume(r, h):
	return(np.pi*((r+(h/2))**3-(r-(h/2))**3)/0.75)
#####################################################################################################
def det_alpha(anm, mu, w):
	anp = anm - mu*w
	return(anp)
#####################################################################################################
def det_starting_flux(r, h, sigma_t, incoming_s_flux, Qcell):
	return((2*incoming_s_flux+h*Qcell)/(2+sigma_t*h))
#####################################################################################################
def det_cell_flux(r, h, mu, w, sigma_t, ap_angle, am_angle, incoming_s_flux, incoming_a_flux, src, vals=False):
	# sigma_t is the total cross section
	# ap and am_angle are the angular differencing coefficients alpha_n+1/2 and alpha_n-1/2a
	Ap, Am = det_shell_surface_area(r, h)
	V = det_shell_volume(r, h)
	if(mu<0):
		denom = 2*abs(mu)*Am+(4.0/w)*(Ap-Am)*ap_angle+V*sigma_t
	else:
		denom = 2*abs(mu)*Ap+(4.0/w)*(Ap-Am)*ap_angle+V*sigma_t
	a = abs(mu)*(Ap+Am)*incoming_s_flux
	b = (2.0/w)*(Ap-Am)*(ap_angle+am_angle)*incoming_a_flux
	c = V*src
	flux = (a+b+c)/denom
	return(flux)
#####################################################################################################
def spatial_difference(prev_cell_flux, prev_half_flux):
	return(2*prev_cell_flux-prev_half_flux)
#####################################################################################################
def angular_difference(flux_ni, flux_nmi):
	flux_npi = 2*flux_ni - flux_nmi
	return(flux_npi)
#####################################################################################################
def det_core_volume_mesh(R, nCells):
	# R is the radius of the spherical core
	# determine cell-center x values
	cell_volume = np.pi*R**3/0.75/nCells
	# specify radial positions of the spherical shells
	# that bound the equal volume cells
	shells = np.zeros(int(nCells)+1)
	# cell_center contains the radial positions of the
	# centers of the spherical shells
	cell_center = np.zeros(int(nCells))
	shells[1]=(cell_volume*0.75/np.pi)**0.333333333333333333333
	cell_center[0]=shells[1]/2
	for i in range(2, nCells+1):
		shells[i]=((0.75*cell_volume/np.pi)+shells[i-1]**3)**0.333333333333333333333
		cell_center[i-1]=shells[i-1]+(shells[i]-shells[i-1])/2
	return(cell_center, shells)
#####################################################################################################
def det_refl_volume_mesh(core_radius, refl_thickness, nCells=1):
	refl_volume = np.pi*((core_radius+refl_thickness)**3-core_radius**3)/0.75
	if(refl_volume!=0):
		cell_volume = refl_volume/nCells
		# specify radial positions of the spherical shells
		# that bound the equal volume cells
		shells = np.zeros(int(nCells)+1)
		# cell_center contains the radial positions of the
		# centers of the spherical shells in the reflector
		cell_centers = np.zeros(int(nCells))
		shells[0]=core_radius
		for i in range(1, nCells+1):
			shells[i]=((0.75*cell_volume/np.pi)+shells[i-1]**3)**0.333333333333333333333
			cell_centers[i-1]=shells[i-1]+(shells[i]-shells[i-1])/2
			# print(4*np.pi*(shells[i]**3-shells[i-1]**3)/3)
		return(cell_centers, shells)
	else:
		return(np.array([]), np.array([]))
#####################################################################################################
def det_core_uniform_mesh(R, nCells):
	# R is the radius of the spherical core
	# determine cell-center x values
	# the spacing between the cells
	h = float(R)/nCells
	shells = np.zeros(int(nCells)+1)
	# centers contains the radial positions of the
	# centers of the spherical shells
	centers = np.zeros(int(nCells))
	for i in range(0, nCells):
		centers[i]=(i*h)+h/2
		shells[i+1]=(i+1)*h
	return(centers, shells)
#####################################################################################################
def det_refl_uniform_mesh(core_radius, refl_thickness, nCells=1):
	refl_volume = np.pi*((core_radius+refl_thickness)**3-core_radius**3)/0.75
	if(refl_volume!=0):
		h = refl_thickness/float(nCells)
		# specify radial positions of the spherical shells
		# that bound the equal volume cells
		shells = np.zeros(int(nCells)+1)
			# cell_center contains the radial positions of the
		# centers of the spherical shells in the reflector
		cell_centers = np.zeros(int(nCells))
		shells[0]=core_radius
		for i in range(1, nCells+1):
			shells[i]=shells[i-1]+h
			cell_centers[i-1]=shells[i]-h/2
		return(cell_centers, shells)
	else:
		return(np.array([]), np.array([]))
#####################################################################################################
def integrate_quadrature(angularflux, muVector, wVector):
	x = {}
	for i in range(0, len(muVector)):
		x[muVector[i]] = angularflux[muVector[i]]*wVector[i]
	return(np.sum(x.values(), axis=0)/2)
#####################################################################################################
def update_scattering_source(order, angularflux, sigma_s, muVector, wVector):
#---------------------------------------------------------------------------------------------------#
	#Function to generate the flux moments
	def det_legendre_moments(L, angflux, muVec, wVec):
		phi_l = {}
		for l in range(0, L+1):
			Pl = legendre(l)
			x = {}
			for i in range(0, len(muVec)):
				x[muVec[i]] = angflux[muVec[i]]*wVec[i]*Pl(muVec[i])
			phi_l[l]=np.sum(x.values(), axis=0)/2
		return(phi_l)
#---------------------------------------------------------------------------------------------------#
	phi_l = det_legendre_moments(order, angularflux, muVector, wVector)
	s_source = {}
	for n in range(0, len(muVector)):
		qn={}
		for l in range(0, order+1):
			Pl = legendre(l)
			qn[l] = (2*l+1)*sigma_s*Pl(muVector[n])*phi_l[l]
		s_source[muVector[n]] = np.sum(qn.values(), axis=0)/2
	return(s_source)
#####################################################################################################
def init_uniform_source(Q, nCells, cell_centers, cell_bound, muVec):
	ext = np.ones(nCells)
	for i in range(0, nCells):
		inner, outer = cell_bound[i], cell_bound[i+1]
		h = outer - inner
		ext[i]=Q*4*np.pi*(outer**3-inner**3)/det_shell_volume(cell_centers[i],h)/3
	angular_ext = {}
	for i in range(0, len(muVec)):
		angular_ext[muVec[i]]=ext
	return(angular_ext)
#####################################################################################################
def init_point_source(Q, nCells, cell_centers, cell_bound, muVec, pt_src_loc):
	# pt_src_dir specifies the direction that the point source sends particles
	ext = np.zeros(nCells)
	ext[pt_src_loc]=Q
	angular_ext = {}
	for i in range(0, len(muVec)):
		angular_ext[muVec[i]]=ext
	return(angular_ext)
#####################################################################################################
def init_wdd(cell_centers, cell_shells, mu1, Q, xs_tot, hf):
	half_angle_flux = np.zeros(len(cell_centers))
	for j in range(0, len(cell_centers)):
		r = cell_centers[::-1][j]
		h = cell_shells[::-1][j]-cell_shells[::-1][j+1]
		qc = Q[mu1][len(cell_centers)-j-1]#/len(muvec)
		xs = xs_tot[len(cell_centers)-j-1]
		haf = det_starting_flux(r, h, xs, hf, qc)
		half_angle_flux[len(cell_centers)-j-1] = haf
		#print('incoming hf = '+str(hf))
		#print('haf = '+str(haf))
		hf = spatial_difference(haf, hf)
		#print ('outgoing hf = '+str(hf))
	#incoming_center_hf = hf
	return(half_angle_flux, hf)
#####################################################################################################
def wdd(cell_center_flux, half_angle_flux, hf, cell_centers, cell_shells, muvec, wvec, Q, xs_tot, incoming_hf_at_center):
#---------------------------------------------------------------------------------------------------#
	# START IN THE NEGATIVE DIRECTION, AT THE OUTSIDE OF THE SPHERE MOVING IN TOWARD THE CENTER
	# Iterate in angle
	am_angle = 0
	#temp = {}
	for i in range(0, len(muvec)/2):
		mu = muvec[i]
		w = wvec[i]
		ap_angle = det_alpha(am_angle, mu, w)
		#print('mu = '+str(mu))
#---------------------------------------------------------------------------------------------------#
		# Iterate in space
		for j in range(0, len(cell_centers)):
			r = cell_centers[::-1][j]
			h = cell_shells[::-1][j]-cell_shells[::-1][j+1]
			haf = half_angle_flux[len(cell_centers)-j-1]
			qc = Q[mu][len(cell_centers)-j-1]
			xs = xs_tot[len(cell_centers)-j-1]
			#test[mu][len(cell_centers)-j-1]=qc/((2*mu/r)+(2*(ap_angle-am_angle)/w/r)+xs)
			#print('r = '+str(r)+'\t incoming hf = '+str(hf))
			cf = det_cell_flux(r, h, mu, w, xs, ap_angle, am_angle, hf, haf, qc, vals=False)
			hf = spatial_difference(cf, hf)
			# NEGATIVE FLUX FIXUP
			# if the outgoing flux from a cell is less than 0, recalculate the flux in that
			# cell by setting the incoming flux to 0. Only need to do this one time in 1D.
			if(hf < 0.0):
				cf = det_cell_flux(r, h, mu, w, xs, ap_angle, am_angle, 0.0, haf, qc)
				hf = spatial_difference(cf, hf)
			cell_center_flux[mu][len(cell_centers)-j-1] = cf
			#print('r = '+str(r)+'\t cf = '+str(cf))
			#print('r = '+str(r)+'\t outgoing hf = '+str(hf))
			half_angle_flux = angular_difference(cell_center_flux[mu], half_angle_flux)
			#temp[mu]=half_angle_flux
			am_angle = ap_angle
		#print('\n')
#---------------------------------------------------------------------------------------------------#
	# NOW STEP IN THE POSITIVE DIRECTION, FROM THE CENTER OF THE SPHERE OUT TO THE EDGE
	# Iterate in angle
	for i in range(len(muvec)/2, len(muvec)):
		mu = muvec[i]
		w = wvec[i]
		hf = incoming_hf_at_center
		#print('mu = '+str(mu))
		ap_angle = det_alpha(am_angle, mu, w)
#---------------------------------------------------------------------------------------------------#
		# Iterate in space
		for j in range(0, len(cell_centers)):
			r = cell_centers[j]
			h = cell_shells[j+1]-cell_shells[j]
			haf = half_angle_flux[j]
			qc = Q[mu][j]
			xs = xs_tot[j]
			#print('r = '+str(r)+'\t incoming hf = '+str(hf))
			#test[mu][j]=qc/((2*mu/r)+(2*(ap_angle-am_angle)/w/r)+xs)
			cf = det_cell_flux(r, h, mu, w, xs, ap_angle, am_angle, hf, haf, qc, vals=False)
			hf = spatial_difference(cf, hf)
			# NEGATIVE FLUX FIXUP
			# if the outgoing flux from a cell is less than 0, recalculate the flux in that
			# cell by setting the incoming flux to 0. Only need to do this one time in 1D.
			if(hf < 0.0):
				#print('negative flux at r = '+str(r))
				cf = det_cell_flux(r, h, mu, w, xs, ap_angle, am_angle, 0.0, haf, qc)
				hf = spatial_difference(cf, 0.0)
			cell_center_flux[mu][j] = cf
			#print('r = '+str(r)+'\t cf = '+str(cf))
			#print('r = '+str(r)+'\t outgoing hf = '+str(hf))
			half_angle_flux = angular_difference(cell_center_flux[mu], half_angle_flux)
			am_angle = ap_angle
			#temp[mu]=half_angle_flux
			#print('\n')
	return(cell_center_flux)#, temp)
#####################################################################################################
def iterate_flux(cell_centers, cell_shells, muvec, wvec, hf, external_source, L, sigma_totv, sigma_sv, tol=1e-4, iteration_cap = 100):
#---------------------------------------------------------------------------------------------------#
	# Function to assess error; because the scattering source varies over angle, I take error
	# to be the maximum difference between old and new scattering sources.
	def check_error(new_s_source, old_s_source):
		error = {}
		for mu in new_s_source:
			error[mu]=abs(np.linalg.norm(new_s_source[mu] - old_s_source[mu]))
		return(max(error.values()))
#---------------------------------------------------------------------------------------------------#
	# Initialize angular flux dictionary
	angularflux = {}
	for i in range(0, len(muvec)):
		angularflux[muvec[i]]=np.zeros(len(cell_centers))
#---------------------------------------------------------------------------------------------------#
	# Initialize scattering source
	scattering_source = update_scattering_source(L, angularflux, sigma_sv, muvec, wvec)
#---------------------------------------------------------------------------------------------------#
	# Initialize source dictionary
	source = {}
	for mu in external_source:
		source[mu]=external_source[mu]+scattering_source[mu]
#---------------------------------------------------------------------------------------------------#
	# Iterate to converge the scattering source
	error = 1.0
	iteration_counter = 0
	while(error > tol and iteration_counter < iteration_cap):
		iteration_counter += 1
		haf, icf = init_wdd(cell_centers, cell_shells, muvec[0], source, sigma_totv, hf)
		angularflux = wdd(angularflux, haf, hf, cell_centers, cell_shells, muvec, wvec, source, sigma_totv, icf)
		new_scattering_source = update_scattering_source(L, angularflux, sigma_sv, muvec, wvec)
		error = check_error(new_scattering_source, scattering_source)
		scattering_source = new_scattering_source
		for mu in source:
			source[mu]=external_source[mu]+scattering_source[mu]
	if(error<tol):
		print('converged in '+str(iteration_counter)+' iterations')
		scalarflux = integrate_quadrature(angularflux, muvec, wvec)
		return(scalarflux, angularflux)
	else:
		print('not converged in '+str(iteration_cap)+' iterations')
	print(angularflux)
#####################################################################################################
def det_pt_src_loc(location, cc):
	a = [abs(r-location) for r in cc]
	idx = cell_centers.index(min(location))
	return(idx)
#####################################################################################################
def initialize_inputs():
	reflector=input.reflector
	sourcetype=input.sourcetype
	if(reflector==False):
		try:
			core_radius = input.core_radius
			nCells = input.nCells
			degree = input.degree
			muvec, wvec = leggauss(degree)
			hf = input.hf
			Q = input.Q
			L = input.L
			sigma_totv = np.ones(nCells)*input.sigma_tot
			sigma_sv = np.ones(nCells)*input.sigma_s
			cell_centers, cell_shells = det_core_uniform_mesh(core_radius, nCells)
			if(sourcetype=='uniform'):
				external_source = init_uniform_source(Q, nCells, cell_centers, cell_shells, muvec)
			elif(sourcetype=='point source'):
				pt_src_loc = input.pt_src_loc
				pt_src_idx = det_pt_src_loc(pt_src_loc, cell_centers)
				external_source = init_point_source(Q, nCells, cell_centers, cell_shells, muvec)
			else:
				print('incorrect source specification, please see comments in input file')
			return(cell_centers, cell_shells, muvec, wvec, hf, external_source, L, sigma_totv, sigma_sv)
		except:
			print('input error, please ensure that all inputs are correct')
	elif(reflector==True):
		try:
			core_radius = input.core_radius
			refl_thickness = input.refl_thickness
			nCells_core = input.nCells_core
			nCells_refl = input.nCells_refl
			degree = input.degree
			muvec, wvec = leggauss(degree)
			sigma_tot_core = input.sigma_tot_core
			sigma_s_core = input.sigma_s_core
			sigma_tot_refl = input.sigma_tot_refl
			sigma_s_refl = input.sigma_s_refl
			Q = input.Q
			hf = input.hf
			L = input.L
			sigma_totv = np.concatenate((np.ones(nCells_core)*sigma_tot_core, np.zeros(nCells_refl)*sigma_tot_refl))
			sigma_sv = np.concatenate((np.ones(nCells_core)*sigma_s_core, np.ones(nCells_refl)*sigma_s_refl))
			core_cell_centers, core_cell_shells = det_core_uniform_mesh(core_radius, nCells_core)
			refl_cell_centers, refl_cell_shells = det_refl_uniform_mesh(core_radius, refl_thickness, nCells_refl)
			cell_centers = np.concatenate((core_cell_centers, refl_cell_centers))
			if(sourcetype=='uniform'):
				external_source = init_uniform_source(Q, nCells_core, core_cell_centers, core_cell_shells, muvec)
			elif(sourcetype=='point source'):
				pt_src_loc = input.pt_src_loc
				pt_src_idx = det_pt_src_loc(pt_src_loc, core_cell_centers)
				external_source = init_point_source(Q, nCells, core_cell_centers, core_cell_shells, muvec)
			else:
				print('incorrect source specification, please see comments in input file')
			for mu in external_source:
				external_source[mu] = np.concatenate((external_source[mu], np.zeros(nCells_refl)))
			return(cell_centers, cell_shells, muvec, wvec, hf, external_source, L, sigma_totv, sigma_sv)
		except:
			print('input error, please ensure that all inputs are correct')
	else:
		print('make sure you either want or do not want a reflector, you cannot have it both ways')
#####################################################################################################




