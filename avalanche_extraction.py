import numpy as np
import matplotlib.pyplot as plt

from param_unwrapper import *
from avalanche_writing import *

#Extracts distribution of avalanches' sizes, given the trajectory of all lattice spins
#Does so by tracking each spin flip (the change of a continuously relaxed spin's sign), and updating an existing avalanche if that spin flip occured adjacent to recently flipped spin in that avalanche
def avalanche_extraction(general_params, spin, time_window, movie, instance_num, T_min=None, T_max=None):

    sz, gamma, delta, zeta, xini, T, dt, transient, param_string = param_unwrapper(general_params)
    
    #Updates T_min and T_max to transient and T, respectively, if no arguments are passed
    if T_min == None:
        T_min = transient
    if T_max == None:
        T_max = T
    
    avalanche_sizes = [] #distribution of spatial avalanche sizes

    #Initializes an array which keeps track of the avalanche history at previous timesteps
    #This will allow us to verify when spin flips are induced by their nearest neighbors (a branching process)
    #Each element characterizes a single avalanche, which has a running size, a list of sites which flipped recently ("active" spins), a list of number of time steps since each active spin flipped, and a list of sites which flipped long ago ("passive" spins)
    recent_avalanches = []

    #Iterates over each timestep of the spin configuration's evolution
    for t in range(np.int(T_min/dt), np.int(T_max/dt)):

        #Generates "movie" file at each timestep to easily visualize spin flips
        if movie:
            movie_array = np.zeros((sz[0], sz[1]))
            plt.imsave(param_string + '/' + str(instance_num) + 'movie_spin' + str(t) + '.png', spin[t], cmap='RdYlBu')

        #Establishes array of spin flips at each time step, with 1 indicating a flip and 0 indicating no flip
        spin_flips = np.where(np.multiply(spin[t+1], spin[t]) < 0 , 1, 0)

        #Only investigates a given timestep if a spin flip occurs
        if any(element == 1 for element in np.ndarray.flatten(spin_flips)) or any(element == -1 for element in np.ndarray.flatten(spin_flips)):

            #Establishes indices at which all flips occur at this timestep
            flipping_array = np.where(spin_flips)

            #Iterates over all flipped spins
            for k in range(len(flipping_array[0])):
                i = flipping_array[0][k] #extracts row index
                j = flipping_array[1][k] #extracts column index
                ###print('Spin flip at [' + str(i) + ', ' + str(j) + ']!')

                #Provides list of nearest-neighbors given a spin's location in the lattice
                if (i==0 and j==0): #top left corner
                    neighbor_list = [[i, j+1], [i+1, j], [i, sz[1]-1], [sz[0]-1, j]]
                elif (i==0 and j==sz[1]-1): #top right corner
                    neighbor_list = [[i, 0], [i+1, j], [i, j-1], [sz[0]-1, j]]
                elif (i==sz[0]-1 and j==0): #bottom left corner
                    neighbor_list = [[i, j+1], [0, j], [i, sz[1]-1], [i-1, j]]
                elif (i==sz[0]-1 and j==sz[1]-1): #bottom right corner
                    neighbor_list = [[i, 0], [0, j], [i, j-1], [i-1, j]]
                elif (i==0): #top edge
                    neighbor_list = [[i, j+1], [i+1, j], [i, j-1], [sz[0]-1, j]]
                elif (i==sz[0]-1): #bottom edge
                    neighbor_list = [[i, j+1], [0, j], [i, j-1], [i-1, j]]
                elif (j==0): #left edge
                    neighbor_list = [[i, j+1], [i+1, j], [i, sz[1]-1], [i-1, j]]
                elif (j==sz[1]-1): #right edge
                    neighbor_list = [[i, 0], [i+1, j], [i, j-1], [i-1, j]]
                else:
                    neighbor_list = [[i, j+1], [i+1, j], [i, j-1], [i-1, j]]

                ###print('Neighbors: ' + str(neighbor_list))

                #Indicates that a spin has flipped in the movie visualization
                if movie:
                    movie_array[i][j] = 1

                no_neighbors = True
                existing_avalanches = []
                for avalanche in recent_avalanches:
                    #Updates existing avalanches if the [i, j] flip has just occurred adjacent to one of their "active" (recently flipped) spins 
                    if any(site in avalanche[1] for site in neighbor_list):
                        #Verifies that the flipping spin is not already an active or passive spin in this avalanche
                        if ([i, j] not in avalanche[1]) and ([i, j] not in avalanche[3]):
                            ###print('Neighbor found: ' + str(avalanche))
                            avalanche[0] += 1 #adds one to that avalanche's size
                            avalanche[1].append([i, j]) #Adds the newly flipped site to its corresponding list in the avalanche it is a part of
                            avalanche[2].append(0) #Specifies that 0 time has passed since this spin flipped
                            no_neighbors = False
                            existing_avalanches.append(avalanche)

                #Creates new, combined avalanche which combines the size and recent spins flips from two constituent avalanches
                #Process allows for the combination of up to 4 avalanches at the same time step
                if len(existing_avalanches) > 1:
                    combined_avalanche = [1, [[i, j]], [0], max([existing_avalanches[i][3] for i in range(len(existing_avalanches))]), []] #initializes the new, combined avalanche
                    for avalanche in existing_avalanches:
                        #Adds the flipping sites and duration since each flip of all constituent avalanches to the combined avalanche (active and passive)
                        for k in range(len(avalanche[1])):
                            if (avalanche[1][k] != [i, j]) and (avalanche[1][k] not in combined_avalanche[1]): #doesn't add duplicates
                                combined_avalanche[0] += 1
                                combined_avalanche[1].append(avalanche[1][k])
                                combined_avalanche[2].append(avalanche[2][k])
                    for avalanche in existing_avalanches:
                        for k in range(len(avalanche[3])):
                            if (avalanche[3][k] not in combined_avalanche[3]) and (avalanche[3][k] not in combined_avalanche[1]): #doesn't add duplicates, taking care to not add a "passive" spin if it's active in a previously merged avalanche
                                combined_avalanche[0] += 1
                                combined_avalanche[3].append(avalanche[3][k])

                        recent_avalanches.remove(avalanche) #removes all of the avalanches which have now combined from the recent_avalanches list
                    recent_avalanches.append(combined_avalanche) #replaces them with the new, combined avalanche

                if no_neighbors:
                    ###print('Neighbor not found!')
                    recent_avalanches.append([1, [[i, j]], [0], []]) #creates a new avalanches of initial size 1, at site [i, j], in which the [i, j] site flip occurred 0 timesteps ago, and which has no "passive" spins


        ###print('Recent avalanches before updating: ' + str(recent_avalanches))

        #Saves an black and white image of the lattice at every timestep, with black corresponding to a spin flip at time t
        if movie:
            plt.imsave(param_string + '/' + str(instance_num) + 'movie_flip' + str(t) + '.png', movie_array, cmap='binary')

        #Updating procedure for all avalanches at the end of each timestep
        for avalanche in recent_avalanches:
            indices_to_remove = []
            for i in range(len(avalanche[2])):
                #Updates the list of indices corresponding to spins in each avalanche which haven't had any neighbor flips in a sufficiently long time, such that these spins are now considered "inactive"
                if avalanche[2][i] >= int(time_window/dt):
                    indices_to_remove.append(i)
                else:
                    avalanche[2][i] += 1
            #Removes each inactive spin, along with its duration indicator, from the list of active spins, adding them to the list of passive spins
            #(replaces each active spin with "-1", and then removes each "-1")
            for index in indices_to_remove:
                avalanche[3].append(avalanche[1][index])
                avalanche[1][index] = -1
                avalanche[2][index] = -1
            for i in range(len(indices_to_remove)):
                avalanche[1].remove(-1)
                avalanche[2].remove(-1)

        #If no spins are still "active" (its been sufficiently long since any neighbor of theirs flipped) in a given avalanche, extract its final size/duration
        for avalanche in recent_avalanches[:]:
            if avalanche[1] == [] and avalanche[2] == []:
                ###print('Avalanche ends, updating avalanche sizes!')
                avalanche_sizes.append(avalanche[0]) #adds final size of a given avalanche to the meta-counter
                recent_avalanches.remove(avalanche) #removes the avalanche from the list of recent avalanches

        ###print('Recent avalanches at end of timestep (T = ' + str(float((t+1)*dt)) + '): ' + str(recent_avalanches))
        ###print('Current avalanche sizes: ' + str(avalanche_sizes))


    #Write remaining avalanches at end of simulation to avalanche_sizes
    for avalanche in recent_avalanches:
        avalanche_sizes.append(avalanche[0])

    #Write avalanche sizes to a file for easy extraction/plotting later
    avalanche_writing(param_string, 'time_window' + str(time_window) + '_' + str(T_min) + '_to_' + str(T_max), avalanche_sizes)
 
    print('Instance ' + str(instance_num) + ' avalanches extracted')

    return avalanche_sizes
