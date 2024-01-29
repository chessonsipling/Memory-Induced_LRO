import time

from spin_dynamics import *
from spin_extraction import *
from spin_plotting import *
from avalanche_extraction import *
from avalanche_plotting import *

#Spin dynamics/plotting and avalanche extraction (per instance) and avalanche plotting for (for entire ensemble)
#####################################################################################################################

#Initialize parameters
general_params = [[32, 32], 0.1, 0.8, 3, 0.5, 512, 0.25, 20, 'DATA, RAW']
num_instances = 10
#Spin dynamics parameters
check_initialization = False
spin_tracker = False
#Spin plotting parameters
num_of_spins = 10
plot_duration = 128
#Avalanche extraction/plotting params
time_window = 30*0.25
flip_axis = 0.0
movie = False
fit = False

for d in [0.65, 0.8, 0.95]:
    general_params[2] = d
    for z in [1.5, 3, 5]:
        general_params[3] = z
        for i in range(1, num_instances + 1):

            spin_dynamics(general_params, check_initialization, spin_tracker, i)
            spin, final_time = spin_extraction(general_params, i)
            if i <= 5:
                spin_plotting(general_params, spin, num_of_spins, plot_duration, i)
            for tw in [12*0.25, 15*0.25, 20*0.25]:
                time_window = tw
                avalanche_extraction(general_params, spin, time_window, flip_axis, movie, i)

        avalanche_plotting(general_params, 'spatial', 12*0.25, flip_axis, fit)
        avalanche_plotting(general_params, 'temporal', 12*0.25, flip_axis, fit)
        avalanche_plotting(general_params, 'spatial', 15*0.25, flip_axis, fit)
        avalanche_plotting(general_params, 'temporal', 15*0.25, flip_axis, fit)
        avalanche_plotting(general_params, 'spatial', 20*0.25, flip_axis, fit)
        avalanche_plotting(general_params, 'temporal', 20*0.25, flip_axis, fit)
#####################################################################################################################
