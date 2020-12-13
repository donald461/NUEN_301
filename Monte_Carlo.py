import numpy as np

# slab description: width and total sigma (macroscopic cross section)
width = 20
sigt_t = 0.2
# numerical parameters
n_bins = 12
bin_width = width/n_bins

# number of histories to follow
neutron_histories = int(input("Input numebr of neutron histories to run: "))
# array to record track length left by neutrons
flux = np.zeros(n_bins)
# create value to count how many neutrons leak out
n_leak = 0

for i in range(neutron_histories):
    absorbed = 0

    # assume random number given
    x = np.random.uniform(0,1,1)
    # distance to be traveld
    distance = -1/sigt_t * np.log(x)
    
    if distance > width:
        # if the neutron gets out of material
        flux[:] += bin_width
        n_leak += 1
    else:
        n_bins_traversed = int(np.floor(distance/bin_width))
        flux[0:n_bins_traversed] += bin_width
        # remainder of distance
        distance_remainder = distance - bin_width * n_bins_traversed
        flux[n_bins_traversed] += distance_remainder
        
# compute fraction of leaked neutons
frac = n_leak/neutron_histories
print('Fraction of neutrons that leked out of slab: ')        
print(frac)
# transform track length tally ito a flux
flux /= bin_width 
# normalize statistics this is now I(x)/I(0)
# Recall that I(x)/I(0) should be equal to: exp(-sigt*x)
flux /= neutron_histories

# plots
import matplotlib.pyplot as plt
x = np.linspace(0,width,n_bins+1)
xx = x[0:-1] + bin_width/2
plt.plot(xx, flux)
plt.grid()










