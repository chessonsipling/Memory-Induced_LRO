import time

from spin_dynamics import *
from spin_extraction import *
from spin_plotting import *
from avalanche_extraction import *
from avalanche_plotting import *

#Initialize parameters
gamma = 0.25
###delta =  0.8
###zeta = 3
xini = 0.5
T = 256
dt = 0.25
transient = 20
num_instances = 5
#Spin dynamics parameters
check_initialization = False
spin_tracker = False
#Spin plotting parameters
num_of_spins = 10
plot_duration = 256
#Avalanche extraction/plotting params
###time_window = 15*dt
flip_axis = 0.0
movie = False
fit = None

#Spin dynamics/plotting and avalanche extraction (per instance) and avalanche plotting for (for entire ensemble)
#####################################################################################################################
for delta in [0.6]: #[0.6, 0.75, 0.9]:
    for zeta in [2]: #[2, 3, 4]:
        for sz in [[32, 32], [64, 64], [96, 96]]:

            general_params = [sz, gamma, delta, zeta, xini, T, dt, transient, 'DATA_RAW']
                
            for i in range(1, num_instances + 1):
                spin_dynamics(general_params, check_initialization, spin_tracker, i)
                spin, final_time = spin_extraction(general_params, i)
                spin_plotting(general_params, spin, num_of_spins, plot_duration, i)
                avalanche_extraction(general_params, spin, 5*dt, flip_axis, movie, i)
                avalanche_extraction(general_params, spin, 15*dt, flip_axis, movie, i)
                avalanche_extraction(general_params, spin, 60*dt, flip_axis, movie, i)

            avalanche_plotting(general_params, 'spatial', 5*dt, flip_axis, fit)
            avalanche_plotting(general_params, 'temporal', 5*dt, flip_axis, fit)
            avalanche_plotting(general_params, 'spatial', 15*dt, flip_axis, fit)
            avalanche_plotting(general_params, 'temporal', 15*dt, flip_axis, fit)
            avalanche_plotting(general_params, 'spatial', 60*dt, flip_axis, fit)
            avalanche_plotting(general_params, 'temporal', 60*dt, flip_axis, fit)
#####################################################################################################################

#Comparing avalanche size distributions across different sized lattices
#####################################################################################################################
        common_general_params = [gamma, var1, var2, xini, T, dt, transient, 'DATA_RAW']

        finite_size = [True, 1.5, 0.0]

        plot_all_avalanches([[32, 32], [64, 64], [96, 96]], common_general_params, ['spatial', 5*dt, 0.0, False], finite_size)
        plot_all_avalanches([[32, 32], [64, 64], [96, 96]], common_general_params, ['temporal', 5*dt, 0.0, False], finite_size)
        plot_all_avalanches([[32, 32], [64, 64], [96, 96]], common_general_params, ['spatial', 15*dt, 0.0, False], finite_size)
        plot_all_avalanches([[32, 32], [64, 64], [96, 96]], common_general_params, ['temporal', 15*dt, 0.0, False], finite_size)
        plot_all_avalanches([[32, 32], [64, 64], [96, 96]], common_general_params, ['spatial', 60*dt, 0.0, False], finite_size)
        plot_all_avalanches([[32, 32], [64, 64], [96, 96]], common_general_params, ['temporal', 60*dt, 0.0, False], finite_size)
#####################################################################################################################
