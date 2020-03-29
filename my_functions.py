import numpy as np
import scipy as sc
import scipy.optimize as opt

pi = np.pi

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
	if a>=0	:orb_freq = (1/(2*pi))*(M/r**3)**0.5 * (1/(1+a*(M/r)**1.5))
	else 	:orb_freq =-(1/(2*pi))*(M/r**3)**0.5 * (1/(1+a*(M/r)**1.5))

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
	per_freq = orb_freq*(1-(1- 6*(M/r) - 3*a**2*(M/r)**2 + 8*a*(M/r)**1.5)**0.5)
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
	nod_freq = orb_freq* (1 - (1+3*a**2*(M/r)**2-4*a*(M/r)**1.5)**0.5)
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
	r = opt.newton(root_func,r0, args=pars[1:])
	return r

