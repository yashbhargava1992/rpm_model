import numpy as np
import matplotlib.pyplot as plt
import my_functions as mf
import time

parallelising =True
parallelising =False

if parallelising: 
	from joblib import Parallel, delayed
	import multiprocessing as mpg

import itertools as it
import os
import time

freqs = np.loadtxt("../pha_from_pds/freqs_l3_sel_20200410.txt",unpack=True)		# l3 has to be nodal
pds_avg_info = np.loadtxt("../segs_of_pds_avg.txt",unpack=True)

#print pds_avg_info
#print freqs[3]
G 	= 6.67e-11		# in SI units
M_sun	= 1.989e30		# In SI units
c	= 2.99e8		# In SI units


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


### USE ONLY nu 123
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
if parallelising:num_cores = mpg.cpu_count()

#out_dir = "mass_spin_sampling_20200415_nu1_nod/"

#out_dir="mass_spin_sampling_20200422_para_v5_nu0_nod_nu2_per_nu3_orb/"
out_dir = 'analytical_test1/'
os.system("mkdir -p {}".format(out_dir))
# use different guesses to get i
#for i in range(1):


## Using the terminology of Ingram and Motta 2014
## nu1: nod; nu2: per; nu3: orb
Gamma = (1-nu2/nu3)**2
Delta = (1-nu1/nu3)**2
#print Gamma,Delta
beta = c**3/(2*np.pi*G*M_sun)

r = 2/3.*((6-Delta-5*Gamma+2*np.sqrt(2*(Delta-Gamma)*(3-Delta-2*Gamma)))/(Delta+Gamma-2)**2)

a = r**1.5/4.*(Delta+Gamma-2+6.0/r)

weird_spin =  np.where(a>0.998)
print segs[weird_spin]
M = beta/(nu3*(r**1.5+a))

plt.plot (segs,r)
plt.plot (segs[weird_spin],r[weird_spin], 'or')
#plt.show()
plt.clf()


plt.plot (segs,a)
plt.plot (segs[weird_spin],a[weird_spin], 'or')
plt.show()
plt.clf()


plt.plot (segs,M)
plt.plot (segs[weird_spin],M[weird_spin], 'or')
#plt.show()
plt.clf()


A = (Gamma+Delta-2)
brack = A+6.0/r
prec = r**1.5/4
plt.plot(segs,brack)
plt.plot(segs,prec)
plt.show()


