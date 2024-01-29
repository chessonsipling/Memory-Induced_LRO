import numpy as np

from param_unwrapper import *

#Extracts spins from existing spin evolution file ([NUM].txt)
#First index extracts timestep, next index/indices give(s) lattice position
def spin_extraction(general_params, instance_num):

    sz, gamma, delta, zeta, xini, T, dt, transient, param_string = param_unwrapper(general_params)

    s_file = open(param_string + '/' + str(instance_num) + '.txt', 'r')
    s = s_file.readlines()
    final_time = float(s.pop(-1)) #"pops" time that evolution stopped from the last line in the spin file

    #Reformats spin data into a list of lists
    spin_intermediate = [s[index].split(',')[:-1] for index in range(len(s))]
    spin_intermediate = [[float(spin_intermediate[i][j]) for j in range(len(spin_intermediate[i]))] for i in range(len(spin_intermediate))]

    #If the array is 2-dimensional, reshapes the spins so that they are formatted in rows and columns (will make these arrays 3-dimensional, including time)
    if len(sz) == 2:
        spin = [np.reshape(spin_intermediate[i], (sz[0], sz[1])) for i in range(len(spin_intermediate))]
        spin = np.array(spin)
    elif len(sz) == 1:
    	spin = np.array(spin_intermediate)

    print('Instance ' + str(instance_num) + ' spins extracted')

    return spin, final_time
