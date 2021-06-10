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

import emcee as mc
import corner as cr
import multiprocess as mpg


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


#~ mass_guess = np.linspace(7.0,10.0,100)		# In solar masses
#~ spin_guess = np.linspace(0.4,0.998,600) 		#For now only positive spins are considered. The dimensionless spin paramter J/Mc2

#~ mass_guess = np.atleast_1d([7.76])

sr_guess = 10
#~ sig = 1-1e-2

#~ mass_spin_all = np.zeros((len(mass_guess),len(spin_guess)))
#~ if parallelising:num_cores = mpg.cpu_count()

out_dir = 'bayesian_test2_parallelise/'
os.system("mkdir -p {}".format(out_dir))
os.system("cp rpm_model_solving.py my_function*py {}".format(out_dir))

#~ chi_2_total = np.zeros((len(mass_guess),len(spin_guess)))

#~ nu1_sel = nu1
#~ nu1_e_sel = nu1_e
#~ print len(nu1_sel)

ndim, nwalkers = 2,100
iters = 1e3

guess = [8,0]


pos = [guess+ 1e-4 * np.random.randn(ndim) for i in range (nwalkers)]

start_time = time.time()

args_for_prior = [7.76,9.27, 0.4, 0.998]

with mpg.Pool() as pool:
	sampler = mc.EnsembleSampler(nwalkers,ndim, mf.lnprob, args = (freqs, args_for_prior, 'uniform', mf.total_chisq), threads=4, pool=pool)
	sampler.run_mcmc(pos, iters)

print ("EMCEE done. Time required:\t", time.time()-start_time)

samples = sampler.chain[:,100:,:].reshape(-1,ndim)

np.save(out_dir+'samples_text.npy', samples)

#~ np.savetxt(out_dir+"chi_2.dat", chi_2_total)
#~ plt.imshow(np.log10(chi_2_total.T), origin='low',extent=[np.min(mass_guess),np.max(mass_guess),np.min(spin_guess),np.max(spin_guess)],aspect='auto')
#pchi_2_perchi_2_perlt.imshow(mass_spin_all.T, origin='low',extent=[np.min(mass_guess),np.max(mass_guess),np.min(spin_guess),np.max(spin_guess)],aspect='auto')

#~ plt.colorbar()
#~ plt.xlabel("Mass")
#~ plt.ylabel("Spin")
#plt.xscale('log')
#plt.yscale('log')
#plt.savefig(out_dir+"mass_spin_distribution_all_pds_20200417.pdf")
#~ plt.show()
#~ plt.clf()
