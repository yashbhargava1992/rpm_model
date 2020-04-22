import numpy as np
import matplotlib.pyplot as plt
import glob


list_of_mass_spin_files = np.sort(glob.glob("mass_spin_sampling_20200419_para_v4_nu0_nod_nu2_per_nu3_orb/*_mass_spin_flag.dat"))

data_all = np.loadtxt(list_of_mass_spin_files[0])
#data_all/=np.sum(data_all)
for i,filename in enumerate(list_of_mass_spin_files):
	if i==0: continue
	data  =np.loadtxt(filename)
	data_all+=data#/np.sum(data)

data_all/=np.shape(list_of_mass_spin_files)[0]

print np.shape(data_all),np.where(data_all>=0.7)
plt.imshow(data_all,origin='lower',extent=[0,0.998,2,20],aspect='auto')
plt.xlabel("Spin")
plt.ylabel("Mass")
plt.colorbar()
#plt.savefig("mass_spin_distribution_overall_20200417.pdf")
plt.show()
		

