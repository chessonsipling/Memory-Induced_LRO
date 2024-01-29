import numpy as np
import random
import matplotlib.pyplot as plt

from param_unwrapper import *

#Plots the evolution of num_of_spins random spins in the lattice
#This function allows for plotting spins in which time evolution is complete (with necessary [NUM].txt file)
def spin_plotting(general_params, spin, num_of_spins, plot_duration, instance_num):

    sz, gamma, delta, zeta, xini, T, dt, transient, param_string = param_unwrapper(general_params)

    #Establishes time array over which we'll plot spin data
    time_data = np.array([i*dt for i in range(int(plot_duration/dt))])

    #Plots the spin evolution profile of n randomly selected [i, j] lattice sites
    i_list = random.sample(range(len(spin[0])), num_of_spins)
    j_list = random.sample(range(len(spin[0][0])), num_of_spins)
    for k in range(num_of_spins):
        try:
            plt.plot(time_data, np.array([spin[l][i_list[k]][j_list[k]] for l in range(int(plot_duration/dt))]))
        except:
            print('Plot duration out of range!')

    #Plots vertical line indicating the chosen transient cutoff
    try:
        plt.axvline(x = time_data[int(transient/dt)], color = 'red', linestyle = 'dashed')
    except:
        print('Transient out of range!')

    plt.xlabel('Time')
    plt.ylabel('Spin')
    plt.title('Spin Evolution')
    plt.savefig(param_string + '/' + str(instance_num) + '.jpg')
    plt.clf()

    print('Instance ' + str(instance_num) + ' plotted')
