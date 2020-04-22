import numpy as np
import matplotlib.pyplot as plt
import my_functions as mf

from joblib import Parallel, delayed
import multiprocessing as mpg

import itertools as it
import os
import time

freqs = np.loadtxt("../pha_from_pds/freqs_l4_sel_20200419.txt",unpack=True)
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

mass_guess = np.linspace(2.0,20.0,60)		# In solar masses
spin_guess = np.linspace(0.0,0.998,50) 		#For now only positive spins are considered. The dimensionless spin paramter J/Mc2

#mass = 7
#spin = 0.5 

# The looping over the mass and spin shall be done using iterables. 

##### Make a plot of all the freqs for a given mass and spin to see which freq could correspond to which eq.

r_guess = 10
sig = 1-1e-2

mass_spin_all = np.zeros((len(mass_guess),len(spin_guess)))
num_cores = mpg.cpu_count()


out_dir="mass_spin_sampling_20200422_para_v5_nu0_nod_nu2_per_nu3_orb/"
os.system("mkdir -p {}".format(out_dir))
# use different guesses to get i
#for i in range(1):
for i in range(66,len(nu0)):

	r_orb_arr = list(np.zeros((len(mass_guess),len(spin_guess))))
	r_per_arr = list(np.zeros((len(mass_guess),len(spin_guess))))
	r_nod_arr = list(np.zeros((len(mass_guess),len(spin_guess))))
#print np.shape(r_orb_arr)
	flag_sel = list(np.zeros((len(mass_guess),len(spin_guess))))
	mass_spin_test = np.zeros((len(mass_guess),len(spin_guess)))
		#print mass,spin,r_guess,r_orb,r_per,r_nod
	#	if not(np.isnan(r_nod)) : print r_guess,r_orb,r_per,r_nod
	nu0_doub = [nu0[i],nu0_e[i]]
	nu1_doub = [nu1[i],nu1_e[i]]
	nu2_doub = [nu2[i],nu2_e[i]]
	nu3_doub = [nu3[i],nu3_e[i]]
	
	begin = time.time()

	outs = Parallel(n_jobs=4)(delayed(mf.radius_compute_and_compare)(nu0_doub,nu2_doub,nu3_doub,mass_guess,spin_guess,ms) for ms in it.product(range(len(mass_guess)),range(len(spin_guess))) )
	outs = np.array(outs)
	flag_sel=outs[:,-1]
	print time.time()-begin, "s passed for an observation"
	flag_sel = np.reshape(flag_sel,(len(mass_guess),len(spin_guess)))
	print i, np.shape(flag_sel), np.sum(flag_sel),
	plot_flag = np.sum(flag_sel)>0
	print plot_flag
	flag_sel = np.array(flag_sel)
	flag_sel = flag_sel.astype(bool)
	#ratio_rad_orb_nod = r_orb_arr/r_nod_arr
	#ratio_rad_per_nod = r_per_arr/r_nod_arr
	
	mass_spin_meshgrid = np.meshgrid(mass_guess,spin_guess,indexing='ij')
	
	mass_spin_test[flag_sel] = 1  
	mass_spin_all += mass_spin_test			# Storing all mass, spin pairs recorded in this iteration
	#mass_spin_test[~flag_sel] = np.nan
	np.savetxt(out_dir+"{:02}_mass_spin_flag.dat".format(i),mass_spin_test)
#	hist_2d = np.histogram2d(,mass_guess,bins=10,weights=mass_spin_test)
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
		print segs[i], number_pds
		plt.title("{0}: {1}".format(int(segs[i]),int(number_pds)))
		plt.ylim(0,0.998)
		plt.xlim(np.min(mass_guess),np.max(mass_guess))
		plt.savefig(out_dir+"{}.png".format(i))
		#plt.show()

	plt.clf()

plt.imshow(mass_spin_all.T, origin='low',extent=[np.min(mass_guess),np.max(mass_guess),np.min(spin_guess),np.max(spin_guess)],aspect='auto')
plt.colorbar()
plt.xlabel("Mass")
plt.ylabel("Spin")
#plt.xscale('log')
#plt.yscale('log')
plt.savefig(out_dir+"mass_spin_distribution_all_pds_20200417.pdf")
#plt.show()
plt.clf()
