import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os

from get_string import get_sz_string
from distribution import *
from param_unwrapper import *

#Reads avalanche distribution data from a .txt file and returns it as useful x and y data, to be plotted
#Also, plots scale-free/gaussian fits, if desired
def single_avalanche_plot(param_string, characterizing_text, fit):

    #Opens the corresponding avalanche distribution file and extracts this data to a list
    avalanche_file = open(param_string + '/avalanche_distrb_' + characterizing_text + '.txt', 'r')
    avalanche_sizes = avalanche_file.readlines()
    avalanche_sizes = [int(element) for element in avalanche_sizes]

    #Stops if there are no avalanches to plot
    if avalanche_sizes == []:
        return

    #Extracts and prints an array with all represented avalanche sizes for this set of instances
    sizes = list(range(1, max(avalanche_sizes) + 1))

    #Extracts and prints the number of avalanches of each size, again for this set of instances
    counts = np.zeros(max(avalanche_sizes))
    for i in range(len(avalanche_sizes)):
        counts[avalanche_sizes[i] - 1] += 1 #offsets, since indexing starts at 0, but the first element in "counts" should correspond to the number of avalanches of size 1
    ###print('Max avalanche size: ' + str(max(avalanche_sizes)))
    ###print('Total number of avalanches: ' + str(int(sum(counts))))

    #Choose a logarithmic binning approach which gives reasonable results (see manuscript for a discussion of this particular binning approach)
    custom_bins = [i for i in range(0, 100, 2)] + [i for i in range(100, 1000, 20)] + [i for i in range(1000, 10000, 200)] + [i for i in range(10000, 100000, 2000)]
    custom_bins = [my_bin for my_bin in custom_bins if my_bin <= max(avalanche_sizes)] #truncates empty bins
    
    hist, bin_edges = np.histogram(avalanche_sizes, bins=custom_bins, density=True)
    bin_centers = [(bin_edges[i] + bin_edges[i+1])/2 for i in range(len(bin_edges) - 1)]

    #Fits data to scale free/gaussian distributions
    if fit == 'sf':
        try:
            #Performs parameter fit on avalanche distribution
            parameters_scale_free, covariance_scale_free = curve_fit(scale_free, bin_centers[3:], hist[3:], p0=[-2, 1], bounds=([-np.inf, 0], [0, np.inf]))
            stddev_error_scale_free = np.sqrt(np.diag(covariance_scale_free))
            print('Scale-Free Params: ' + str(parameters_scale_free))
            print('Standard Errors: ' + str(stddev_error_scale_free))
            #Generates data arrays from the fitted parameters
            x_fit = np.arange(1, max(avalanche_sizes) + 1, 1)
            scale_free_fit = []
            for i in range(len(x_fit)):
                scale_free_fit.append(scale_free(x_fit[i], parameters_scale_free[0], parameters_scale_free[1]))
            plt.plot(x_fit, scale_free_fit, color='red', linestyle='dashed')#, label='Scale Free Fit: a = %1.2f +/- %1.2f' % (parameters_scale_free[0], stddev_error_scale_free[0]))
            plt.text(30, 1e-2, '~ $s^{%1.2f}$' % parameters_scale_free[0], color='red')
        except:
            print('Fitting failed, try different fitting parameter guesses or adjust range over which fitted is attempted')

    elif fit == 'g':
        try:
            parameters_gaussian, covariance_gaussian = curve_fit(gaussian, bin_centers, hist, p0=[0, 2, 5], bounds=([-np.inf, 0, 0], [np.inf, np.inf, np.inf]))
            stddev_error_gaussian = np.sqrt(np.diag(covariance_gaussian))
            print('Gaussian Params: ' + str(parameters_gaussian))
            print('Standard Errors: ' + str(stddev_error_gaussian))
            x_fit = np.arange(1, max(avalanche_sizes) + 1, 1)
            gaussian_fit = []
            for i in range(len(x_fit)):
                gaussian_fit.append(gaussian(x_fit[i], parameters_gaussian[0], parameters_gaussian[1], parameters_gaussian[2]))
            plt.plot(x_fit, gaussian_fit, color='red', linestyle='dashed')#, label='Gaussian Fit: Mu =  %1.2f +/- %1.2f, Sigma = %1.2f +/- %1.2f' % (parameters_gaussian[0], stddev_error_gaussian[0], parameters_gaussian[1], stddev_error_gaussian[1]))
        except:
            print('Fitting failed, try different fitting parameter guesses or adjust range over which fitted is attempted')

    return bin_centers, hist


#Plots avalanche size distributions in a particular ensemble, along with fits (optional)
#Requires parameters describing how avalanche extraction ocurred (time_window)
def avalanche_plotting(general_params, time_window, fit, T_min=None, T_max=None):

    #Some formatting specfications (can be adjusted/added to)
    params = {'font.size':22,
          'legend.loc': 'upper right',
          'text.usetex':True,
          'font.family':'serif',
          'font.serif':['Computer Modern Serif']}
    plt.rcParams.update(params)
    plt.subplots_adjust(bottom=0.25, left=0.25)

    sz, gamma, delta, zeta, xini, T, dt, transient, param_string = param_unwrapper(general_params)

    #Updates T_min and T_max to transient and T, respectively, if no arguments are passed
    if T_min == None:
        T_min = transient
    if T_max == None:
        T_max = T

    characterizing_text = 'time_window' + str(time_window) + '_' + str(T_min) + '_to_' + str(T_max)

    #Extracts the x and y data from .txt files, performing fit if desired
    bin_centers, hist = single_avalanche_plot(param_string, characterizing_text, fit)

    plt.scatter(bin_centers, hist, label='Data')

    plt.xlabel('Avalanche Size (s)')
    plt.ylabel('P(s)')
    plt.legend(fontsize=14)
    if fit != 'g': #does not set axes to log-log iff the fit is to a gaussian
        plt.xscale('log')
        plt.yscale('log')
    plt.savefig(param_string + '/avalanche_distrb_' + characterizing_text + '.jpg')
    plt.clf()


#Plots multiple avalanche distributions of the same parameters on top of one another, but with avalanche extraction happening over different times during the dynamics
def attractive_lro_plot(general_params, time_ranges_list, time_window, fit):

    #Some formatting specfications (can be adjusted/added to)
    params = {'font.size':16,
          'figure.figsize':[6.4, 4.8],
          'legend.loc': 'upper right',
          'text.usetex':True,
          'font.family':'serif',
          'font.serif':['Computer Modern Serif']}
    plt.rcParams.update(params)
    plt.subplots_adjust(bottom=0.25, left=0.25)

    sz, gamma, delta, zeta, xini, T, dt, transient, param_string = param_unwrapper(general_params)
    for i in range(len(time_ranges_list)):
        characterizing_text = 'time_window' + str(time_window) + '_' + str(time_ranges_list[i][0]) + '_to_' + str(time_ranges_list[i][1])

        #Extracts the x and y data from .txt files, performing fit if desired
        if i == len(time_ranges_list) - 1: #only fits last distribution in time_ranges_list (which should be the "most" scale-free)
            bin_centers, hist = single_avalanche_plot(param_string, characterizing_text, fit)
        else:
            bin_centers, hist = single_avalanche_plot(param_string, characterizing_text, False)

        plt.scatter(bin_centers, hist, label=r'$T \in [%d, %d)$' % (time_ranges_list[i][0], time_ranges_list[i][1]))

    plt.xlabel('Avalanche Size (s)')
    plt.ylabel('P(s)')
    plt.legend(fontsize=14)
    ###plt.title(r'$\gamma = $' + str(general_params[1]) + ', ' + r'$N = %d^2$' % sz[0], loc='left', fontsize=16)
    if fit != 'g': #does not set axes to log-log iff the fit is to a gaussian
        plt.xscale('log')
        plt.yscale('log')
    plt.savefig(param_string + '/avalanche_distrb_time_window' + str(time_window) + '_attractiveLRO' '.jpg', dpi=600)
    plt.clf()


#Plots a list of avalanches distributions simultaneously, along with fits (both optional)
#Also can perform a scale-invariance analysis, provided the appropriate scaling exponents (typically "alpha" and "beta") are provided
#Currently, only allows for different sizes to be compared
def plot_all_avalanches(size_list, common_general_params, time_window, fit, finite_size_params, T_min=None, T_max=None):

    #Some formatting specfications (can be adjusted/added to)
    params = {'font.size':26,
          'legend.loc': 'lower left',
          'text.usetex':True,
          'font.family':'serif',
          'font.serif':['Computer Modern Serif']}
    plt.rcParams.update(params)
    plt.subplots_adjust(bottom=0.25, left=0.25)

    #Initializes list of arrays, providing comprehensive data on all avalanche distributions
    #Will be necessary to have this data to perform scale-invariance analysis
    list_of_arraysizes = []
    list_of_bin_centers = []
    list_of_hist = []

    #Updates T_min and T_max to transient and T, respectively, if no arguments are passed
    if T_min == None:
        T_min = common_general_params[6] #transient
    if T_max == None:
        T_max = common_general_params[4] #T

    #Organizes parameters which are common to all ensembles (excluding size) into useful strings
    common_general_string = 'g' + str(common_general_params[0]) + '_d' + str(common_general_params[1]) + '_z' + str(common_general_params[2]) + '_x' + str(common_general_params[3]) + '_T' + str(common_general_params[4]) + '_dt' + str(common_general_params[5])
    characterizing_text = 'time_window' + str(time_window) + '_' + str(T_min) + '_to_' + str(T_max)

    #Tries to make the necessary directory, if it does not exist
    try:
        os.makedirs(common_general_params[7] + '/' + common_general_string)
    except:
        print('Directory already exists!')

    for sz in size_list:

        list_of_arraysizes.append(sz)

        #Extracts the x and y data from .txt files, performing fit if desired
        if fit == 'sf':
            if sz == size_list[-1]:
                bin_centers, hist = single_avalanche_plot(get_param_string(sz, common_general_params[0], common_general_params[1], common_general_params[2], common_general_params[3], common_general_params[4], common_general_params[5], common_general_params[7]), characterizing_text, fit)
            else:
                bin_centers, hist = single_avalanche_plot(get_param_string(sz, common_general_params[0], common_general_params[1], common_general_params[2], common_general_params[3], common_general_params[4], common_general_params[5], common_general_params[7]), characterizing_text, False)
        else:
            bin_centers, hist = single_avalanche_plot(get_param_string(sz, common_general_params[0], common_general_params[1], common_general_params[2], common_general_params[3], common_general_params[4], common_general_params[5], common_general_params[7]), characterizing_text, fit)

        plt.scatter(bin_centers, hist, label=get_sz_string(sz)) #label=r'$%d^2$' % sz[0]

        list_of_bin_centers.append(bin_centers)
        list_of_hist.append(hist)


    plt.xlabel('Avalanche Size (s)')
    plt.ylabel('P(s)')
    plt.legend(fontsize=18)
    plt.title(r'$\gamma = $' + str(common_general_params[0]), loc='left', fontsize=26)
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(common_general_params[7] + '/' + common_general_string + '/' + characterizing_text + '_all_avalanches.jpg',  dpi=600)
    plt.clf()

    #Perform scale invariance analysis
    if finite_size_params[0] == True:
        plt.rcParams.update({'font.size':34})
        for i in range(len(size_list)):
            finite_size_x = [element/((list_of_arraysizes[i][0]*list_of_arraysizes[i][1])**finite_size_params[2]) for element in list_of_bin_centers[i]]
            finite_size_y = [a*b for a,b in zip([element**finite_size_params[1] for element in list_of_bin_centers[i]],list_of_hist[i])]
            plt.scatter(finite_size_x[:], finite_size_y[:], label=get_sz_string(size_list[i]))

        plt.xlabel(r'~$s/L^{%0.2f}$' % float(2*finite_size_params[2]))
        #plt.xlabel('s')
        plt.ylabel(r'~$s^{-%0.2f}$P(s)' % finite_size_params[1])
        plt.xscale('log')
        plt.yscale('log')
        plt.savefig(common_general_params[7] + '/' + common_general_string + '/' + characterizing_text + '_finite_size.jpg',  dpi=600)
        plt.clf()
