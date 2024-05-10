import numpy as np

from param_unwrapper import *

#Extracts data from existing DOF evolution file (instance = str(instance_num) + 'Spin', str(instance_num) + 'MemUp', or str(instance_num) + 'MemRight')
#First index extracts timestep, next 2 indices give 2D lattice position
def data_extraction(general_params, instance):

    sz, gamma, delta, zeta, xini, T, dt, transient, param_string = param_unwrapper(general_params)

    data_file = open(param_string + '/' + instance + '.txt', 'r')
    raw_data = data_file.readlines()
    final_time = float(raw_data.pop(-1)) #"pops" time that evolution stopped from the last line in the spin file

    #Reformats spin-memory data into a list of lists
    data_intermediate = [raw_data[index].split(',')[:-1] for index in range(len(raw_data))]
    data_intermediate = [[float(data_intermediate[i][j]) for j in range(len(data_intermediate[i]))] for i in range(len(data_intermediate))]

    #Reshapes the flatted data so that it's formatted in rows and columns (will make these arrays 2+1 dimensional)
    data = [np.reshape(data_intermediate[i], (sz[0], sz[1])) for i in range(len(data_intermediate))]
    data = np.array(data)

    print('Instance ' + instance + ' data extracted')

    return data, final_time
