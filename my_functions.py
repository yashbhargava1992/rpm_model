import numpy as np
import scipy as sc
import scipy.optimize as opt

pi 	= np.pi
G 	= 6.67e-11		# in SI units
M_sun	= 1.989e30		# In SI units
c	= 2.99e8		# In SI units

def nu_phi (r,*pars):
	"""
	Takes the mass spin and the radius of interest and returns the orbital frequency at that radius in the Kerr geometry. M and a are packed as *pars

	INPUT:
	------------------
	M		: Mass of the BH in solar mass
	a		: Spin of the BH in dimensionless units
	r		: The radius in gravitation radii. GM/c^2

	OUTPUT:
	------------------
	orb_freq	: Orbital frequency as given by formula in Motta et al 2014 GRO 1655-40 paper



	"""
	M,a = pars
	#M = float(M)
	#a = float(a)
	#r = float(r)
	#if a>=0	:orb_freq = (1/(2*pi))*(M/r**3)**0.5 * (1/(1+a*(M/r)**1.5))
	#else 	:orb_freq =-(1/(2*pi))*(M/r**3)**0.5 * (1/(1+a*(M/r)**1.5))
	
	if a>=0	: orb_freq= (1/(2*pi))* (c**6 /((G*M_sun*M)**2*r**3))**0.5 * 1/(1+a*(1.0/r)**1.5)
	else	: orb_freq=-(1/(2*pi))* (c**6 /((G*M_sun*M)**2*r**3))**0.5 * 1/(1+a*(1.0/r)**1.5)
	
	return orb_freq

def nu_per (r,*pars):
	"""
	Takes the mass spin and the radius of interest and returns the periastron frequency at that radius in the Kerr geometry. M and a are packed as *pars

	INPUT:
	------------------
	M		: Mass of the BH in solar mass
	a		: Spin of the BH in dimensionless units
	r		: The radius in gravitation radii. GM/c^2

	OUTPUT:
	------------------
	per_freq	: Periastron requency as given by formula in Motta et al 2014 GRO 1655-40 paper


	"""
	M,a = pars
	orb_freq = nu_phi(r,M,a)
	#per_freq = orb_freq*(1-(1- 6*(M/r) - 3*a**2*(M/r)**2 + 8*a*(M/r)**1.5)**0.5)
	per_freq = orb_freq*(1-(1- 6*(1.0/r) - 3*a**2*(1.0/r)**2 + 8*a*(1.0/r)**1.5)**0.5)
	return per_freq

def nu_nod (r,*pars):
	"""
	Takes the mass spin and the radius of interest and returns the nodal frequency at that radius in the Kerr geometry. M and a are packed as *pars

	INPUT:
	------------------
	M		: Mass of the BH in solar mass
	a		: Spin of the BH in dimensionless units
	r		: The radius in gravitation radii. GM/c^2

	OUTPUT:
	------------------
	nod_freq	: Nodal requency as given by formula in Motta et al 2014 GRO 1655-40 paper


	"""
	M,a = pars
	orb_freq = nu_phi(r,M,a)
	#nod_freq = orb_freq* (1 - (1+3*a**2*(M/r)**2-4*a*(M/r)**1.5)**0.5)
	nod_freq = orb_freq* (1 - (1+3*a**2*(1.0/r)**2-4*a*(1.0/r)**1.5)**0.5)
	return nod_freq





def newton_solver (r0,*pars):
	"""
	The solver will accept the variable to be computed and pars.
	The assumed form to be solved is y - f(x,theta) = 0 where x is the variable to be solved
	y is the given value of the frequency and f is the assumed form of the frequency. Theta is the set of additional parameters required. 
	The pars should be in following order.
	pars[0]		: The function
	pars[1]		: The 'y' value. (in our case would be the observed frequency)
	pars[2:]	: The other parameters the function might depend upon



	"""
	#print r0,pars
	func = pars[0]
	root_func = lambda r0,*par: par[0]-func(r0,*par[1:])
	r = opt.newton(root_func,r0, args=pars[1:],tol=1e-7)
	return r


def radius_compute_and_compare(nu1,nu2,nu3,mass_array,spin_array,ms_index,r_guess=10):
	"""

	This will accept single triplet of nu1<nu2<nu3 and use the RPM model to get a radius estimates. This function is to be called in parallel to do this for multiple mass and spin estimates. 

	INPUT
	nu1		: Triplet in the form par,min_err,plus_err
	nu2		: Triplet in the form par,min_err,plus_err
	nu3		: Triplet in the form par,min_err,plus_err
	mass_array	: Array of mass points from which mass will be sourced
	spin_array	: Array of spin points from which spin will be sourced	
	ms_index	: Index tuple of mass and spin

	OUTPUT
	r_orb
	r_per
	r_nod
	flag_sel
	"""
	m,s = ms_index
	mass = mass_array[m]
	spin = spin_array[s]
	try:
		r_orb_arr = newton_solver(r_guess,nu_phi,nu3[0],mass,spin)
		r_orb_max = newton_solver(r_guess,nu_phi,nu3[0]-nu3[1],mass,spin)
		r_orb_min = newton_solver(r_guess,nu_phi,nu3[0]+nu3[1],mass,spin)
	except RuntimeError:
		#print "Runtime error"
		r_orb_arr = np.nan
		#continue
		r_orb_max = np.nan
		r_orb_min = np.nan
	try:
		r_per_arr = newton_solver(r_guess,nu_per,nu2[0],mass,spin)
		r_per_max = newton_solver(r_guess,nu_per,nu2[0]-nu2[1],mass,spin)
		r_per_min = newton_solver(r_guess,nu_per,nu2[0]+nu2[1],mass,spin)
	except RuntimeError:
		r_per_arr = np.nan
		r_per_max = np.nan
		r_per_min = np.nan
	try:
		r_nod_arr = newton_solver(r_guess,nu_nod,nu1[0],mass,spin)
		r_nod_max = newton_solver(r_guess,nu_nod,nu1[0]-nu1[1],mass,spin)
		r_nod_min = newton_solver(r_guess,nu_nod,nu1[0]+nu1[1],mass,spin)
	except RuntimeError:
		r_nod_arr = np.nan
		r_nod_max = np.nan
		r_nod_min = np.nan
	# Computing the range of the radius. Since frequency is decreasing monotonic func of radius, a low freq will give higher radius. 
			
	# if the largest of the mins is greater than the smallest of the maxs then there is no common interval. So that mass spin pair is not chosen. 
	if np.max([r_orb_min,r_per_min,r_nod_min]) <= np.min([r_orb_max,r_per_max,r_nod_max]): 
		flag_sel = 1
	else: flag_sel = 0

#	return flag_sel
	return r_orb_arr,r_per_arr,r_nod_arr,flag_sel


def chi_sq(y,y_mod,del_y):
	return (y-y_mod)**2/del_y**2
