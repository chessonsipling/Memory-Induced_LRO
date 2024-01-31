import os
import numpy as np
import random
import matplotlib.pyplot as plt

from avalanche_extraction import *
from spin_plotting import *
from param_unwrapper import *

#Initializes and time evolves spin lattice with provided parameters according to nearest-neighbor interactions, along with the presence of an external "memory" field
#Interactions are of strength 1 and randomly generated throughout the lattice as either ferro- or antiferromagnetic
#Memory field (bound between 0 and 1) increases in time as long as spins (pairwise) are in an energetically unfavorable configuration (i.e. s_1 = -1, s_2 = -1, while J_12 = -1 (antiferromagnetic))
def spin_dynamics(general_params, check_initialization, spin_tracker, instance_num):

    sz, gamma, delta, zeta, xini, T, dt, transient, param_string = param_unwrapper(general_params)

    #Defines d, the dimension of the system
    d = len(sz)

    if d == 2:

        #Tries to make the necessary directory, if it does not exist
        try:
            os.makedirs(param_string)
        except:
            print('Directory already exists!')

        #Randomly initializes the spins to either -1 or 1 in an array
        spin = np.array([[random.choice([-1, 1]) for j in range(sz[1])] for i in range(sz[0])])

        #Initializes the spin and memory coupling arrays (J_{ij} and x_{ij})
        #Sets J^{up}_{ij} and J^{right}_{ij} to be random arrays of -1s, and 1s (each INTERACTION is either ferromagnetic or antiferromagnetic)
        #J^{down}_{ij} and J^{left}_{ij} are determined by shifting the J^{up} and J^{right} arrays up and to the right, respectively (to preserve symmetry of interactions)
        #Sets x_{ij} to a constant array of xinis (symmetry constraint is automatically enforced since all entries are xini)
        #This approach makes both J_{ij} and x_{ij} ONLY non-zero for nearest-neighbor interactions
        #NOTE: the pair (i,j) corresponds to the 2-index lattice site
        spin_coupling_up = np.array([[random.choice([-1, 1]) for j in range(sz[1])] for i in range(sz[0])])
        spin_coupling_down = np.roll(spin_coupling_up, -1, axis=0)
        spin_coupling_right = np.array([[random.choice([-1, 1]) for j in range(sz[1])] for i in range(sz[0])])
        spin_coupling_left = np.roll(spin_coupling_right, 1, axis=1)

        #Saves these initialization as binary cmaps (to view initial conditions, verifying randomness, and to verify symmetry of interactions. "Up/down" and "right/left" cmap pairs should be identical up to a shift)
        if check_initialization:
            plt.imsave(param_string + '/' + str(instance_num) + 'initial_spins.png', spin, cmap='binary')
            plt.imsave(param_string + '/' + str(instance_num) + 'initial_spincouplingup.png', spin_coupling_up, cmap='binary')
            plt.imsave(param_string + '/' + str(instance_num) + 'initial_spincouplingdown.png', spin_coupling_down, cmap='binary')
            plt.imsave(param_string + '/' + str(instance_num) + 'initial_spincouplingright.png', spin_coupling_right, cmap='binary')
            plt.imsave(param_string + '/' + str(instance_num) + 'initial_spincouplingleft.png', spin_coupling_left, cmap='binary')

        memory = xini*np.ones((sz[0], sz[1]))
        memory_up = memory
        memory_down = memory
        memory_right = memory
        memory_left = memory #these arrays are identical now, but will take on different values as the memory variables evolve

        #Creates the file that this instance's spin dynamics will be written to and catalogues the initial spins
        spin_file = open(param_string + '/' + str(instance_num) + '.txt', 'a')
        for i in range(sz[0]):
            for j in range(sz[1]):
                spin_file.write(str(spin[i][j]) + ',') #writes the initialized spin values to this file
        spin_file.write('\n')


        t = 0 #initializes t, corresponding to the physical time that has passed
        while np.round(t, 4) < T:

            if spin_tracker:
                #Initializes a feature which keeps track of the number of spins within some region of 0 (those which are "close" to a flip)
                #This is intended to stop the dynamics immediately if a pseudo-steady-state is reached (all spins stop flipping)
                interesting_spins = 0

                #Saves most recent spin and memory values to a file within the corresponding (parameter-based) directory
                for i in range(sz[0]):
                    for j in range(sz[1]):
                        #Adds 1 to the "interesting spins" counter if a spin is within [-0.5, 0.5]
                        if (spin[i][j] <= 0.5) and (spin[i][j] >= -0.5):
                            interesting_spins += 1

                #Stops the dynamics if no spins are "close" to flipping
                #This indicates that a pseudo-steady-state has been found (which we do not care about)
                if (interesting_spins == 0) and (t > 0):
                    t += dt
                    break

            print('Simulating instance ' + str(instance_num) + ', Time: ' + str(np.round(t, 4)) + ' out of ' + str(T))

            #Creates arrays which, at each point, indicate the value of the spin in some direction adjacent to that point
            #For example, the (2, 3) entry in "spin_right" would correspond to the current spin value at (2, 4) in the spin lattice array (since (2, 4) is directly to the right of (2, 3))
            #Since the dynamics are governed by nearest-neighbor interactions, this is a relatively efficient way to evolve the spins
            spin_up = np.roll(spin, 1, axis=0)
            spin_up[0,:] = 0
            spin_down = np.roll(spin, -1, axis=0)
            spin_down[-1,:] = 0
            spin_right = np.roll(spin, -1, axis=1)
            spin_right[:,-1] = 0
            spin_left = np.roll(spin, 1, axis=1)
            spin_left[:,0] = 0

            #Explicitly updates both the spin and memory variables based on their derivatives
            spin_deriv = (np.multiply(spin_coupling_up,spin_up) + np.multiply(spin_coupling_down,spin_down) + np.multiply(spin_coupling_right,spin_right) + np.multiply(spin_coupling_left,spin_left) - zeta*np.multiply(memory_up + memory_down + memory_right + memory_left,spin))
            memory_up_deriv = gamma*((1/2)*(np.multiply(np.multiply(spin_coupling_up,spin),spin_up) + 1) - delta)
            memory_down_deriv = gamma*((1/2)*(np.multiply(np.multiply(spin_coupling_down,spin),spin_down) + 1) - delta)
            memory_right_deriv = gamma*((1/2)*(np.multiply(np.multiply(spin_coupling_right,spin),spin_right) + 1) - delta)
            memory_left_deriv = gamma*((1/2)*(np.multiply(np.multiply(spin_coupling_left,spin),spin_left) + 1) - delta)

            spin = np.clip(spin + spin_deriv*dt, -1, 1) #updates spins, bounding them between [-1, 1]
            memory_up = np.clip(memory_up + memory_up_deriv*dt, 0, 1) #updates memories, bounding them between [0, 1]
            memory_down = np.clip(memory_down + memory_down_deriv*dt, 0, 1)
            memory_right = np.clip(memory_right + memory_right_deriv*dt, 0, 1)
            memory_left = np.clip(memory_left + memory_left_deriv*dt, 0, 1)


            #Saves most recent spin values to a file within the corresponding (parameter-based) directory and to an accessible array
            for i in range(sz[0]):
                for j in range(sz[1]):
                    spin_file.write(str(spin[i][j]) + ',')
            spin_file.write('\n')

            t += dt

        #Writes the timestep at which the dynamics stopped to the corresponding spin file (may be t < T if dynamics stopped early due to "interesting spins" feature)
        spin_file.write(str(t))
        spin_file.close()
