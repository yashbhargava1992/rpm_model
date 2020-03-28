import numpy as np
import matplotlib.pyplot as plt

freqs = np.loadtxt("freqs.txt",skiprows=1,unpack=True)

#print freqs[3]

for i in range(np.shape(freqs)[0]):
	plt.plot(freqs[2],freqs[i],'.')
plt.xscale('log')
plt.yscale('log')
plt.show()

nu0 = freqs[0]
nu1 = freqs[2]	# nu0,1,2,3 are sorted in freq
nu2 = freqs[1]
nu3 = freqs[3]


# Use one frequency and a given Mass and spin. Try to determine r. Use the error on the frequency to get a range of r. 
#For same mass and spin use the other freq and formula to compute r. Find where three 3 lies within the conf interval of each other

### Use these ranges once the code is finalised

#mass_guess = np.arange(5,20,100)		# In solar masses
#spin_guess = np.arange(-0.998,0.998,100) 	# The dimensionless spin paramter J/Mc2



