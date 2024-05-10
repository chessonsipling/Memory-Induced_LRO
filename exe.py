import time

from spin_dynamics import *
from data_extraction import *
from data_plotting import *
from avalanche_extraction import *
from avalanche_plotting import *

##########################PARAMETERS THAT CAN BE CHANGED ON A SIMULATION-BY-SIMULATION BASIS#########################
#PRIMARY SIMULATION PARAMETERS
size_list = [[16, 16], [32, 32], [64, 64]]
gamma = 0.25
delta = 0.75
g = 2.0
xini = 0.5

#SECONDARY SIMULATION PARAMETERS
#Miscellanous parameters
subdirectory = 'DATA_RAW'
num_instances = 100
#Spin dynamics parameters
check_initialization = False
spin_tracker = False
#Plotting parameters
num_spins = 10
#Avalanche extraction/plotting parameters
movie = False
time_ranges_list = [[0, 4], [4, 16], [16, 64], [64, 200]]
fit = 'sf'
#Finite size scaling parameters
fsscaling = True
#####################################################################################################################

#ADDITIONAL \gamma-DEPENDENT SIMULATION PARAMETERS
#gamma: T, dt, transient, plot_duration, time_window
parameter_dict = {#1e-4: [40000, 0.25, 7500, 40000, 125],
                  #1e-3: [4000, 0.25, 750, 4000, 50],
                  1e-2: [1000, 0.25, 100, 1000, 16],
                  1e-1: [200, 0.25, 15, 200, 5],
                  0.15: [200, 0.25, 15, 200, 5],
                  #0.2: [200, 0.25, 15, 200, 5],
                  0.25: [200, 0.25, 15, 200, 5],
                  #0.3: [200, 0.25, 15, 200, 5],
                  #0.35: [200, 0.25, 15, 200, 5],
                  0.4: [200, 0.25, 15, 200, 5],
                  #0.45: [200, 0.25, 15, 200, 5],
                  #0.5: [150, 0.25, 10, 150, 3],
                  1: [100, 0.15, 8, 100, 2],
                  #2: [100, 0.10, 5, 50, 2],
                  1e1: [100, 0.05, 3, 20, 2]}
                  #1e2: [100, 0.03, 3, 10, 2],
                  #1e3: [100, 0.01, 3, 6, 2]} #some gamma-dependent parameter estimates based on previous simulations (SLIGHT CHANGES MAY YIELD BETTER RESULTS)

variable_parameters = parameter_dict[gamma]
T = variable_parameters[0]
dt = variable_parameters[1]
transient = variable_parameters[2]
plot_duration = variable_parameters[3]
time_window = variable_parameters[4]

sf_param_dict = {1e-2: [2.02, 1.00],
                 1e-1: [1.77, 1.00],
                 0.15: [2.12, 1.00],
                 0.25: [1.97, 1.00],
                 0.4: [2.02, 1.00],
                 1: [1.98, 0.0],
                 1e1: [3.80, 0.0]} #some gamma-dependent parameter estimates based on previous simulations (SLIGHT CHANGES MAY YIELD BETTER RESULTS)

sf_parameters = sf_param_dict[gamma]
alpha = sf_parameters[0]
beta = sf_parameters[1]

#Iterates over lattice sizes
for sz in size_list:

    #Structure with general simulation parameters which will be passed to most functions
    general_params = [sz, gamma, delta, g, xini, T, dt, transient, subdirectory]

    #Iterates over instances
    for i in range(1, num_instances + 1):

        spin_dynamics(general_params, check_initialization, spin_tracker, i) #simulates dynamics
        spin, final_time = data_extraction(general_params, str(i)) #extracts spin evolutions
        if i <= 5:
            memory_up, final_time = data_extraction(general_params, str(i) + 'MemUp') #extracts memory evolutions (invoking J_{i j} i<->j symmetry)
            memory_right, final_time = data_extraction(general_params, str(i) + 'MemRight')
            data_plotting(general_params, spin, memory_up, memory_right, num_spins, plot_duration, final_time, i) #plots spin and memory dynamics for a few instances
        avalanche_extraction(general_params, spin, time_window, movie, i) #extracts avalanche distributions from T_min=transient to T_max=T from spin evolutions

        ###for time_range in time_ranges_list:
            ###avalanche_extraction(general_params, spin, time_window, movie, i, time_range[0], time_range[1]) #extracts avalanche distributions over different time ranges

    ###attractive_lro_plot(general_params, time_ranges_list, time_window, fit) #plots avalanche distributions over different time ranges (to show attractiveness of LRO phase)

    avalanche_plotting(general_params, time_window, fit) #plots avalanche distribution for a single size


common_general_params = general_params[1:]
finite_size_params = [fsscaling, alpha, beta]

plot_all_avalanches(size_list, common_general_params, time_window, fit, finite_size_params) #plots avalanche distributions over different sizes (to show persistence of LRO phase in thermodynamic limit)
