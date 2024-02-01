import time

from spin_dynamics import *
from spin_extraction import *
from spin_plotting import *
from avalanche_extraction import *
from avalanche_plotting import *

#Initialize parameters
gamma = 0.1
delta = 0.75
zeta = 2
xini = 0.5
T = 4096
dt = 0.25
transient = 20
###num_instances = 10
#Spin dynamics parameters
check_initialization = False
spin_tracker = False
#Spin plotting parameters
num_of_spins = 10
plot_duration = 512
#Avalanche extraction/plotting params
###time_window = 15*dt
flip_axis = 0.0
movie = False
fit = None
#Finite size scaling parameters
###fsscaling = True
###alpha = 1.8
###beta = 0.0

#Spin dynamics/plotting and avalanche extraction (per instance) and avalanche plotting for (for entire ensemble)
######################################################################################################################################
for sz in [[16, 16], [128, 128]]:

    general_params = [sz, gamma, delta, zeta, xini, T, dt, transient, 'DATA_RAW_CUTOFF']
                
    for i in range(1, 10 + 1):
        spin_dynamics(general_params, check_initialization, spin_tracker, i)
        spin, final_time = spin_extraction(general_params, i)
        if i <= 5:
            spin_plotting(general_params, spin, num_of_spins, plot_duration, i)
        avalanche_extraction(general_params, spin, 15*dt, flip_axis, movie, i)
######################################################################################################################################

#Comparing avalanche size distributions across different sized lattices
######################################################################################################################################
'''for sz in [[16, 16], [32, 32], [64, 64], [96, 96], [128, 128]]:
    avalanche_plotting(general_params, 'spatial', 15*dt, flip_axis, fit)
    avalanche_plotting(general_params, 'temporal', 15*dt, flip_axis, fit)

common_general_params = [gamma, delta, zeta, xini, T, dt, transient, 'DATA_RAW_CUTOFF']

###finite_size = [fsscaling, alpha, beta]

plot_all_avalanches([[16, 16], [32, 32], [64, 64], [96, 96], [128, 128]], common_general_params, ['spatial', 15*dt, flip_axis, 'sf'], [True, 1.69, 0.0])
plot_all_avalanches([[16, 16], [32, 32], [64, 64], [96, 96], [128, 128]], common_general_params, ['temporal', 15*dt, flip_axis, 'sf'], [True, 1.83, 0.0])'''
######################################################################################################################################
