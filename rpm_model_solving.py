import numpy as np
import matplotlib.pyplot as plt
import my_functions as mf


freqs = np.loadtxt("../pha_from_pds/freqs_l3_sel_20200410.txt",unpack=True)
pds_avg_info = np.loadtxt("../segs_of_pds_avg.txt",unpack=True)
#print pds_avg_info
#print freqs[3]

#for i in range(np.shape(freqs)[0]/2):
#	plt.plot(freqs[0],freqs[i*2],'.')
#plt.xscale('log')
#plt.yscale('log')
#plt.show()
#plt.clf()

segs	= freqs[0] 
nu0 	= freqs[1]
nu0_e	= freqs[2]
nu1 	= freqs[5]	# nu0,1,2,3 are sorted in freq
nu1_e	= freqs[6]
nu2 	= freqs[3]
nu2_e	= freqs[4]
nu3 	= freqs[7]
nu3_e	= freqs[8]


# Use one frequency and a given Mass and spin. Try to determine r. Use the error on the frequency to get a range of r. 
#For same mass and spin use the other freq and formula to compute r. Find where three 3 lies within the conf interval of each other

### Use these ranges once the code is finalised

mass_guess = np.linspace(2.0,20.0,600)		# In solar masses
spin_guess = np.linspace(0.0,0.998,500) 		#For now only positive spins are considered. The dimensionless spin paramter J/Mc2

#mass = 7
#spin = 0.5 

# The looping over the mass and spin shall be done using iterables. 

##### Make a plot of all the freqs for a given mass and spin to see which freq could correspond to which eq.

r_guess = 10
sig = 1-1e-2

mass_spin_all = np.zeros((len(mass_guess),len(spin_guess)))


# use different guesses to get i
#for i in range(1):
for i in range(len(nu0)):

	r_orb_arr = np.zeros((len(mass_guess),len(spin_guess)))
	r_per_arr = np.zeros((len(mass_guess),len(spin_guess)))
	r_nod_arr = np.zeros((len(mass_guess),len(spin_guess)))
#print np.shape(r_orb_arr)
	flag_sel = np.zeros((len(mass_guess),len(spin_guess)))
	mass_spin_test = np.zeros((len(mass_guess),len(spin_guess)))
	for m,mass in enumerate(mass_guess):
		for s,spin in enumerate(spin_guess):
			try:
				r_orb_arr[m,s] = mf.newton_solver(r_guess,mf.nu_phi,nu3[i],mass,spin)
				r_orb_max = mf.newton_solver(r_guess,mf.nu_phi,nu3[i]-nu3_e[i],mass,spin)
				r_orb_min = mf.newton_solver(r_guess,mf.nu_phi,nu3[i]+nu3_e[i],mass,spin)
			except RuntimeError:
				#print "Runtime error"
				r_orb_arr[m,s] = np.nan
				#continue
				r_orb_max = np.nan
				r_orb_min = np.nan
			try:
				r_per_arr[m,s] = mf.newton_solver(r_guess,mf.nu_per,nu2[i],mass,spin)
				r_per_max = mf.newton_solver(r_guess,mf.nu_per,nu2[i]-nu2_e[i],mass,spin)
				r_per_min = mf.newton_solver(r_guess,mf.nu_per,nu2[i]+nu2_e[i],mass,spin)
			except RuntimeError:
				r_per_arr[m,s] = np.nan
				r_per_max = np.nan
				r_per_min = np.nan

			try:
				r_nod_arr[m,s] = mf.newton_solver(r_guess,mf.nu_nod,nu1[i],mass,spin)
				r_nod_max = mf.newton_solver(r_guess,mf.nu_nod,nu1[i]-nu1_e[i],mass,spin)
				r_nod_min = mf.newton_solver(r_guess,mf.nu_nod,nu1[i]+nu1_e[i],mass,spin)
			except RuntimeError:
				r_nod_arr[m,s] = np.nan
				r_nod_max = np.nan
				r_nod_min = np.nan
			# Computing the range of the radius. Since frequency is decreasing monotonic func of radius, a low freq will give higher radius. 
			
		# if the largest of the mins is greater than the smallest of the maxs then there is no common interval. So that mass spin pair is not chosen. 
			if np.max([r_orb_min,r_per_min,r_nod_min]) <= np.min([r_orb_max,r_per_max,r_nod_max]): 
				flag_sel[m,s] = 1

			

		#print mass,spin,r_guess,r_orb,r_per,r_nod
	#	if not(np.isnan(r_nod)) : print r_guess,r_orb,r_per,r_nod
	print i, np.shape(flag_sel), np.sum(flag_sel)
	plot_flag = np.sum(flag_sel)>0
	print plot_flag
	flag_sel = flag_sel.astype(bool)
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
	#print i, np.shape(ind_all),
	#if np.shape(ind_all)[1] == 0:
	#	print "Too tight constraint. No fit found. "
	#else :
	#	print np.shape(ind_all)[1], "pairs of mass and spin satisfy, plotting position of the first of them"
	#	plt.plot(mass_guess[ind_all[0]][0],spin_guess[ind_all[1]][0],'.')
	#	plt.annotate(i,(mass_guess[ind_all[0]][0],spin_guess[ind_all[1]][0]))
	#plt.plot( mass_guess[ind_all[0]],spin_guess[ind_all[1]], '.')
	mass_spin_meshgrid = np.meshgrid(mass_guess,spin_guess,indexing='ij')
	
	mass_spin_test[flag_sel] = 1  
	mass_spin_all += mass_spin_test			# Storing all mass, spin pairs recorded in this iteration
	#mass_spin_test[~flag_sel] = np.nan
	np.savetxt("mass_spin_sampling_20200415_nu1_nod/{:02}_mass_spin_flag.dat".format(i),mass_spin_test)
#	hist_2d = np.histogram2d(,mass_guess,bins=10,weights=mass_spin_test)
#	print hist_2d
	#plt.plot(mass_spin_meshgrid[0][flag_sel],mass_spin_meshgrid[1][flag_sel],'.',alpha=0.2)
	#plt.contourf(mass_spin_meshgrid[0][flag_sel],mass_spin_meshgrid[1][flag_sel])
	if plot_flag: 
		#plt.contour(spin_guess,mass_guess,mass_spin_test)
		plt.hist2d(mass_spin_meshgrid[0][flag_sel],mass_spin_meshgrid[1][flag_sel], bins=20)
		plt.colorbar()
		plt.xlabel("Mass")
		plt.ylabel("Spin")
		obsid = 1200120000+segs[i]
		number_pds = pds_avg_info[1][np.where(pds_avg_info[0]==obsid)]
		plt.title("{0}: {1}".format(int(segs[i]),int(number_pds)))
		plt.ylim(0,0.998)
		plt.xlim(np.min(mass_guess),np.max(mass_guess))
#		plt.savefig("{}.png".format(i))
		#plt.show()

	plt.clf()

plt.imshow(mass_spin_all.T, origin='low',extent=[np.min(mass_guess),np.max(mass_guess),np.min(spin_guess),np.max(spin_guess)])
plt.colorbar()
plt.xlabel("Mass")
plt.ylabel("Spin")
#plt.xscale('log')
#plt.yscale('log')
#plt.savefig("mass_spin_distribution_all_pds_20200414.pdf")
#plt.show()
plt.clf()
