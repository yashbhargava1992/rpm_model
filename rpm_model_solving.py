import numpy as np
import matplotlib.pyplot as plt
import my_functions as mf


freqs = np.loadtxt("../freqs.txt",skiprows=1,unpack=True)

#print freqs[3]

for i in range(np.shape(freqs)[0]):
	plt.plot(freqs[2],freqs[i],'.')
plt.xscale('log')
plt.yscale('log')
#plt.show()
plt.clf()


nu0 = freqs[0]
nu1 = freqs[2]	# nu0,1,2,3 are sorted in freq
nu2 = freqs[1]
nu3 = freqs[3]


# Use one frequency and a given Mass and spin. Try to determine r. Use the error on the frequency to get a range of r. 
#For same mass and spin use the other freq and formula to compute r. Find where three 3 lies within the conf interval of each other

### Use these ranges once the code is finalised

mass_guess = np.linspace(3.0,20.0,200)		# In solar masses
spin_guess = np.linspace(0.,0.998,100) 		#For now only positive spins are considered. The dimensionless spin paramter J/Mc2

#mass = 7
#spin = 0.5 

# The looping over the mass and spin shall be done using iterables. 

##### Make a plot of all the freqs for a migven mass and spin to see which freq could correspond to which eq.

r_guess = 10

for mass in mass_guess:
	for spin in spin_guess:
		try:
			r_orb = mf.newton_solver(r_guess,mf.nu_phi,nu2[0],mass,spin)
		except RuntimeError:
			#print "Runtime error"
			r_orb = np.nan
			#continue
		try:
			r_per = mf.newton_solver(r_guess,mf.nu_per,nu1[0],mass,spin)
		except RuntimeError:
			r_per = np.nan

		try:
			r_nod = mf.newton_solver(r_guess,mf.nu_nod,nu0[0],mass,spin)
		except RuntimeError:
			r_nod = np.nan

		#print mass,spin,r_guess,r_orb,r_per,r_nod
		if not(np.isnan(r_nod)) : print r_guess,r_orb,r_per,r_nod
