## In this branch, the nodal freq will be used to compute the radius for a given pair of mass and spin. 
## Once we obtain the radius then the radius will be used to compute the other frequencies for the mass and spin.
## Compute a chi^2 for the these two frequencies and choose the mass and spin pair which satisfies the criteria of minimum chi^2



import numpy as np
import matplotlib.pyplot as plt
import my_functions as mf
import time
import pandas as pd 
parallelising =True
parallelising =False

if parallelising: 
	from joblib import Parallel, delayed
	import multiprocessing as mpg

import itertools as it
import os
import time

#freqs = np.loadtxt("../pha_from_pds/freqs_l4_sel_20200419.txt",unpack=True)
#~ freqs = np.loadtxt("../pha_from_pds/freqs_l3_sel_20200513_onlycorr.txt",unpack=True)
#freqs = np.loadtxt("../pha_from_pds/freqs_l3_sel_20200609_width_as_error_v2.txt",unpack=True)

pds_avg_info = np.loadtxt("../segs_of_pds_avg.txt",unpack=True)


#~ segs	= freqs[0] 
#~ nu0 	= freqs[1]
#~ nu0_e	= freqs[2]
#~ nu1 	= freqs[5]	# nu0,1,2,3 are sorted in freq
#~ nu1_e	= freqs[6]
#~ nu2 	= freqs[3]
#~ nu2_e	= freqs[4]
#~ nu3 	= freqs[7]
#~ nu3_e	= freqs[8]


freqs = pd.read_csv("../pha_from_pds/freqs_l3_sel_20200513_onlycorr.txt",header=None, delim_whitespace=True, names = ['segs','nu0','nu0_e','nu2','nu2_e','nu1','nu1_e','nu3','nu3_e'])


mass_guess = np.linspace(7.0,10.0,100)		# In solar masses
spin_guess = np.linspace(0.4,0.998,600) 		#For now only positive spins are considered. The dimensionless spin paramter J/Mc2

mass_guess = np.atleast_1d([7.76])

r_guess = 10
sig = 1-1e-2

mass_spin_all = np.zeros((len(mass_guess),len(spin_guess)))
if parallelising:num_cores = mpg.cpu_count()

out_dir = 'global_fit_single_mass_7p76/'
os.system("mkdir -p {}".format(out_dir))
os.system("cp rpm_model_solving.py my_function*py {}".format(out_dir))

chi_2_total = np.zeros((len(mass_guess),len(spin_guess)))

nu1_sel = nu1
nu1_e_sel = nu1_e
print len(nu1_sel)

if not parallelising:
	begin = time.time()
	for ms in it.product(range(len(mass_guess)),range(len(spin_guess))):
			m,s = ms
			mass,spin = mass_guess[m],spin_guess[s]
			r_nod_arr = np.zeros(len(nu1_sel))
			for i in range(0,len(nu1_sel)):
				try:
					r_nod_arr[i] = mf.newton_solver(r_guess,mf.nu_nod,nu1_sel[i],mass,spin)
				except RuntimeError:
					r_nod_arr[i] = np.nan
			nu_per_model = mf.nu_per(r_nod_arr,mass,spin)
			nu_orb_model = mf.nu_phi(r_nod_arr,mass,spin)
			chi_2_per = np.sum(mf.chi_sq(nu2,nu_per_model,nu2_e))
			chi_2_orb = np.sum(mf.chi_sq(nu3,nu_orb_model,nu3_e))
			chi_2_total[m,s] = chi_2_per+chi_2_orb

else:		
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
	#print i, np.shape(flag_sel), np.sum(flag_sel),
	plot_flag = np.sum(flag_sel)>0
	#print plot_flag
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
		#print segs[i], number_pds
		plt.title("{0}: {1}".format(int(segs[i]),int(number_pds)))
		plt.ylim(0,0.998)
		plt.xlim(np.min(mass_guess),np.max(mass_guess))
		plt.savefig(out_dir+"{}.png".format(i))
		#plt.show()

	plt.clf()

np.savetxt(out_dir+"chi_2.dat", chi_2_total)
plt.imshow(np.log10(chi_2_total.T), origin='low',extent=[np.min(mass_guess),np.max(mass_guess),np.min(spin_guess),np.max(spin_guess)],aspect='auto')
#pchi_2_perchi_2_perlt.imshow(mass_spin_all.T, origin='low',extent=[np.min(mass_guess),np.max(mass_guess),np.min(spin_guess),np.max(spin_guess)],aspect='auto')

plt.colorbar()
plt.xlabel("Mass")
plt.ylabel("Spin")
#plt.xscale('log')
#plt.yscale('log')
#plt.savefig(out_dir+"mass_spin_distribution_all_pds_20200417.pdf")
plt.show()
plt.clf()
