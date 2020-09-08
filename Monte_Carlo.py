import numpy as np

# slab description: width and total sigma (macroscopic cross section)
width = 0.5
sigt_t = 1
sigt_a = .4
# numerical parameters
n_bins = 10
bin_width = width/n_bins

# number of histories to follow
neutron_histories = 1000
# array to record track length left by neutrons
flux = np.zeros(n_bins)

for i in range(neutron_histories):
    absorbed = 0

    # assume random number given
    x = np.random.uniform(0,1,1)
    a = np.random.uniform(0,1,1)
    d = np.random.uniform(0,1,1)
    # distance to be traveld
    distance = -1/sigt_t * np.log(x)
    while (absorbed != 1):
        if a <= sigt_a:
            # neutron is absorbed
            if distance > width:
                flux[:] += bin_width
                break
            
            else:
                bins_fully_traversed = int(np.floor(distance/bin_width))
                flux[0:bins_fully_traversed] += bin_width
                
                #  track_remainder = distance % bin_width
                track_remainder = distance - bin_width * bins_fully_traversed
                flux[bins_fully_traversed] += track_remainder
                break
        else:
            # neutron scatters
            if distance > width:
                flux[:] += bin_width
                break
            else:
                bins_fully_traversed = int(np.floor(distance/bin_width))
                flux[0:bins_fully_traversed] += bin_width
                
                #  track_remainder = distance % bin_width
                track_remainder = distance - bin_width * bins_fully_traversed
                flux[bins_fully_traversed] += track_remainder
                if d < 0.5:
                    # neutrons going to the left
                    n = np.random.uniform(0,1,1)
                    distance_new = -1/sigt_t * np.log(x)
                    if distance - distance_new < 0:
                        flux[0:distance_new] += bin_width
                        break
                    else:
                        #neutron interacts again
                        bins_fully_traversed_new = bins_fully_traversed - int(np.floor(distance_new/bin_width))
                        flux[bins_fully_traversed_new:bins_fully_traversed] += bin_width
                            
                        #  track_remainder = distance % bin_width
                        track_remainder = distance_new - bin_width * bins_fully_traversed_new
                        flux[bins_fully_traversed_new] += track_remainder  
                        distance = distance_new
                        bins_fully_traversed_new = bins_fully_traversed
                else:
                    # neutrons going to the right
                    n = np.random.uniform(0,1,1)
                    distance_new = -1/sigt_t * np.log(x)
                    if distance - distance_new < 0:
                        flux[0:distance_new] += bin_width
                        break
                    else:
                        bins_fully_traversed_new = bins_fully_traversed + int(np.floor(distance_new/bin_width))
                        flux[bins_fully_traversed:bins_fully_traversed_new] += bin_width
                            
                        #  track_remainder = distance % bin_width
                        track_remainder = distance_new - bin_width * bins_fully_traversed_new
                        flux[bins_fully_traversed_new] += track_remainder  
                        distance = distance_new
                        bins_fully_traversed_new = bins_fully_traversed
                
        
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
