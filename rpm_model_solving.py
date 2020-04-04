import numpy as np
import matplotlib.pyplot as plt
import my_functions as mf


freqs = np.loadtxt("../pha_from_pds/freqs_l3_sel_20200404.txt",unpack=True)

#print freqs[3]

for i in range(np.shape(freqs)[0]/2):
	plt.plot(freqs[0],freqs[i*2],'.')
plt.xscale('log')
plt.yscale('log')
plt.show()
plt.clf()


nu0 = freqs[0]
nu1 = freqs[2]	# nu0,1,2,3 are sorted in freq
nu2 = freqs[1]
nu3 = freqs[3]


# Use one frequency and a given Mass and spin. Try to determine r. Use the error on the frequency to get a range of r. 
#For same mass and spin use the other freq and formula to compute r. Find where three 3 lies within the conf interval of each other

### Use these ranges once the code is finalised

mass_guess = np.linspace(7.0,10.0,100)		# In solar masses
spin_guess = np.linspace(0.0,0.998,500) 		#For now only positive spins are considered. The dimensionless spin paramter J/Mc2

#mass = 7
#spin = 0.5 

# The looping over the mass and spin shall be done using iterables. 

##### Make a plot of all the freqs for a migven mass and spin to see which freq could correspond to which eq.

r_guess = 10
sig = 1-1e-2

# use different guesses to get i
for i in range(len(nu0)):

	r_orb_arr = np.zeros((len(mass_guess),len(spin_guess)))
	r_per_arr = np.zeros((len(mass_guess),len(spin_guess)))
	r_nod_arr = np.zeros((len(mass_guess),len(spin_guess)))
#print np.shape(r_orb_arr)

	for m,mass in enumerate(mass_guess):
		for s,spin in enumerate(spin_guess):
			try:
				r_orb_arr[m,s] = mf.newton_solver(r_guess,mf.nu_phi,nu3[i],mass,spin)
			except RuntimeError:
				#print "Runtime error"
				r_orb_arr[m,s] = np.nan
				#continue
			try:
				r_per_arr[m,s] = mf.newton_solver(r_guess,mf.nu_per,nu2[i],mass,spin)
			except RuntimeError:
				r_per_arr[m,s] = np.nan

			try:
				r_nod_arr[m,s] = mf.newton_solver(r_guess,mf.nu_nod,nu0[i],mass,spin)
			except RuntimeError:
				r_nod_arr[m,s] = np.nan

		#print mass,spin,r_guess,r_orb,r_per,r_nod
	#	if not(np.isnan(r_nod)) : print r_guess,r_orb,r_per,r_nod

	ratio_rad_orb_nod = r_orb_arr/r_nod_arr
	ratio_rad_per_nod = r_per_arr/r_nod_arr
	
	ind_orb_nod = np.where((ratio_rad_orb_nod< 1.0/sig) & (ratio_rad_orb_nod>sig))
	ind_per_nod = np.where((ratio_rad_per_nod< 1.0/sig) & (ratio_rad_per_nod>sig))
	ind_all	    = np.where((ratio_rad_per_nod< 1.0/sig) & (ratio_rad_per_nod>sig) & (ratio_rad_orb_nod< 1.0/sig) & (ratio_rad_orb_nod>sig) )

	#print np.shape(ind_orb_nod), ind_orb_nod[0], np.shape(ind_per_nod), ind_per_nod
	#print mass_guess[ind_orb_nod[0]],spin_guess[ind_orb_nod[1]]
	#print mass_guess[ind_per_nod[0]],spin_guess[ind_per_nod[1]]

	#plt.plot( mass_guess[ind_orb_nod[0]],spin_guess[ind_orb_nod[1]],'oC0')
	
	#plt.plot( mass_guess[ind_per_nod[0]],spin_guess[ind_per_nod[1]],'dC1')
	print i, np.shape(ind_all),
	if np.shape(ind_all)[1] == 0:
		print "Too tight constraint. No fit found. "
	else :
		print np.shape(ind_all)[1], "pairs of mass and spin satisfy, plotting position of the first of them"
		plt.plot(mass_guess[ind_all[0]][0],spin_guess[ind_all[1]][0],'.')
		plt.annotate(i,(mass_guess[ind_all[0]][0],spin_guess[ind_all[1]][0]))
	plt.plot( mass_guess[ind_all[0]],spin_guess[ind_all[1]], '.')
	plt.xlabel("Mass")
	plt.ylabel("Spin")
	#plt.show()
	#plt.clf()
plt.xscale('log')
plt.yscale('log')
plt.show()
