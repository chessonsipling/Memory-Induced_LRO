import numpy as np
import random
import matplotlib.pyplot as plt

from param_unwrapper import *

#Plots the evolution of num_spins random spins in the lattice, and 1 random spin s_i with its 4 memories x_{i j} for a given instance_num
def data_plotting(general_params, spin, memory_up, memory_right, num_spins, plot_duration, final_time, instance_num):

    sz, gamma, delta, zeta, xini, T, dt, transient, param_string = param_unwrapper(general_params)


    params = {'font.size':16,
        'figure.figsize':[6.4, 4.8],
        'text.usetex':True,
        'font.family':'serif',
        'font.serif':['Computer Modern Serif']}
    plt.rcParams.update(params)
    plt.subplots_adjust(bottom=0.15, left=0.15)

    #Establishes time array over which we'll plot data
    time_data = np.array([i*dt for i in range(int(plot_duration/dt))])

    #Plots the spin evolution profile at n randomly selected [i, j] lattice sites
    i_list = random.sample(range(len(spin[0])), num_spins)
    j_list = random.sample(range(len(spin[0][0])), num_spins)
    for k in range(num_spins):
        try:
            plt.plot(time_data, np.array([spin[l][i_list[k]][j_list[k]] for l in range(int(plot_duration/dt))]))
        except:
            print('Plot duration out of range!')
            plt.plot(time_data, np.array([spin[l][i_list[k]][j_list[k]] for l in range(int(final_time/dt))]))

    #Plots spin values about which about spins will oscillate for \gamma >> 1 (quasiperiodic orbits)
    ###plt.axhline(y=-np.sqrt(2*delta - 1), color='r', linestyle='dashed')
    ###plt.axhline(y=np.sqrt(2*delta - 1), color='r', linestyle='dashed')

    plt.xlabel('Time (arb. units)')
    plt.ylabel('Spin')
    plt.savefig(param_string + '/' + str(instance_num) + 'Spins.jpg', dpi=600)
    plt.clf()


    params = {'font.size':26,
          'figure.figsize':[8.0, 5.0],
          'text.usetex':True,
          'font.family':'serif',
          'font.serif':['Computer Modern Serif']}
    plt.rcParams.update(params)
    plt.subplots_adjust(bottom=0.25, left=0.15)

    #Choses a single random lattice site (which, for simplicity, is not at the edges)
    i = random.randint(1, sz[0]-2)
    j = random.randint(1, sz[1]-2)

    #Plots 1 s_i spin and its 4 x_{i j} memories (for the j nearest-neighbors of i)
    try:
        plt.plot(time_data, np.array([memory_up[k][i][j] for k in range(int(plot_duration/dt))]), alpha=0.6, color='red', label=r'$x_{i j}$') #x_{i j} memory, for j above i
        plt.plot(time_data, np.array([memory_up[k][i+1][j] for k in range(int(plot_duration/dt))]), alpha=0.6, color='red') #x_{i j} memory, for j below i
        plt.plot(time_data, np.array([memory_right[k][i][j] for k in range(int(plot_duration/dt))]), alpha=0.6, color='red') #x_{i j} memory, for j to the right of i
        plt.plot(time_data, np.array([memory_right[k][i][j-1] for k in range(int(plot_duration/dt))]), alpha=0.6, color='red') #x_{i j} memory, for j to the left of i
        plt.plot(time_data, np.array([spin[k][i][j] for k in range(int(plot_duration/dt))]), color='blue', label=r'$s_i$') #s_i spin
    except:
        print('Plot duration out of range!')
        plt.plot(time_data, np.array([memory_up[k][i][j] for k in range(int(final_time/dt))]), alpha=0.6, color='red', label=r'$x_{i j}$') #x_{i j} memory, for j above i
        plt.plot(time_data, np.array([memory_up[k][i+1][j] for k in range(int(final_time/dt))]), alpha=0.6, color='red') #x_{i j} memory, for j below i
        plt.plot(time_data, np.array([memory_right[k][i][j] for k in range(int(final_time/dt))]), alpha=0.6, color='red') #x_{i j} memory, for j to the right of i
        plt.plot(time_data, np.array([memory_right[k][i][j-1] for k in range(int(final_time/dt))]), alpha=0.6, color='red') #x_{i j} memory, for j to the left of i
        plt.plot(time_data, np.array([spin[k][i][j] for k in range(int(final_time/dt))]), color='blue', label=r'$s_i$') #s_i spin

    plt.xlabel('Time (arb. units)')
    plt.legend()
    plt.savefig(param_string + '/' + str(instance_num) + 'SpinsAndMemories.jpg', dpi=600)
    plt.clf()

    print('Instance ' + str(instance_num) + ' plotted')
