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
        counts[avalanche_sizes[i] - 1] += 1 #offsets, since indexing starts at 0, but the first element in "counts" should correspond to the number of avalanches of size 0
    ###print('Max avalanche size: ' + str(max(avalanche_sizes)))
    ###print('Total number of avalanches: ' + str(int(sum(counts))))
    percentages = counts/sum(counts) #normalizes counts array


    #Fits data to scale free/gaussian distributions
    if fit == 'sf':
        try:
            #Performs parameter fit on avalanche distribution
            parameters_scale_free, covariance_scale_free = curve_fit(scale_free, sizes[10:75], percentages[10:75], p0=[-1.7, 1], bounds=([-np.inf, 0], [0, np.inf]))
            stddev_error_scale_free = np.sqrt(np.diag(covariance_scale_free))
            #Generates data arrays from the fitted parameters
            x_fit = np.arange(1, max(avalanche_sizes) + 1, 1)
            scale_free_fit = []
            for i in range(len(x_fit)):
                scale_free_fit.append(scale_free(x_fit[i], parameters_scale_free[0], parameters_scale_free[1]))
            plt.plot(x_fit, scale_free_fit, color='red', linestyle='dashed')#, label='Scale Free Fit: a = %1.2f +/- %1.2f' % (parameters_scale_free[0], stddev_error_scale_free[0])) #plots this fitted data alongside the raw data of (sizes, counts)
            plt.text(45, 1e-3, '~ $s^{%1.2f}$' % parameters_scale_free[0], color='red')
        except:
            print('Fitting failed, try different fitting parameter guesses')

    elif fit == 'g':
        try:
            parameters_gaussian, covariance_gaussian = curve_fit(gaussian, sizes[1:-10], percentages[1:-10], p0=[sz[0]/10, sz[0]/20, 15000], bounds=([0, 0, 0], [np.inf, np.inf, np.inf]))
            stddev_error_gaussian = np.sqrt(np.diag(covariance_gaussian))
            x_fit = np.arange(1, max(avalanche_sizes) + 1, 1)
            gaussian_fit = []
            for i in range:
                gaussian_fit.append(gaussian(x_fit[i], parameters_gaussian[0], parameters_gaussian[1], parameters_gaussian[2]))
            plt.plot(x_fit, gaussian_fit, color='red', linestyle='dashed')#, label='Gaussian Fit: Mu =  %1.2f +/- %1.2f, Sigma = %1.2f +/- %1.2f' % (parameters_gaussian[0], stddev_error_gaussian[0], parameters_gaussian[1], stddev_error_gaussian[1]))
        except:
            print('Fitting failed, try different fitting parameter guesses')


    return sizes, percentages


#Plots avalanche size distributions (either spatial or temporal) in a particular ensemble, along with fits (optional)
#Requires parameters describing how avalanche extraction ocurred (time_window, flip_axis)
def avalanche_plotting(general_params, distribution_type, time_window, flip_axis, fit):

    sz, gamma, delta, zeta, xini, T, dt, transient, param_string = param_unwrapper(general_params)
    characterizing_text = 'time_window' + str(time_window) + 'flip_axis' + str(flip_axis) + '_' + distribution_type

    #Extracts the x and y data from .txt files, performing fit if desired
    sizes, percentages = single_avalanche_plot(param_string, characterizing_text, fit)

    plt.scatter(sizes, percentages, label='Data')
    plt.ylim(1e-5, 1)
    plt.xlabel('Avalanche Size')
    plt.ylabel('Percentage of Total Counts')
    plt.legend()
    plt.title('Avalanche Size Distribution (Number of Spin Flips)\n' + characterizing_text)
    if fit != 'g': #does not set axes to log-log iff the fit is to a gaussian
        plt.xscale('log')
        plt.yscale('log')
    plt.savefig(param_string + '/avalanche_distrb_' + characterizing_text + '.jpg')
    plt.clf()


#Plots a list of avalanches distributions simultaneously, along with fits (both optional)
#Also can perform a scale-invariance analysis, provided the appropriate scaling exponents (typically "alpha" and "beta") are provided
#Currently, only allows for different sizes to be compared
def plot_all_avalanches(size_list, common_general_params, common_misc_params, finite_size):

    #Some formatting specfications (can be adjusted/added to)
    plt.rcParams.update({'font.size': 22, 'legend.loc': 'lower left'})
    ###plt.rcParams.update({'font.size': 22, 'legend.loc': 'upper right'})
    plt.subplots_adjust(bottom=0.2, left=0.2)

    #Initializes list of arrays, providing comprehensive data on all avalanche distributions
    #Will be necessary to have this data to perform scale-invariance analysis
    list_of_arraysizes = []
    list_of_sizes = []
    list_of_percentages = []

    #Organizes parameters which are common to all ensembles (excluding size) into useful strings
    common_general_string = 'g' + str(common_general_params[0]) + '_d' + str(common_general_params[1]) + '_z' + str(common_general_params[2]) + '_x' + str(common_general_params[3]) + '_T' + str(common_general_params[4]) + '_dt' + str(common_general_params[5])
    characterizing_text = 'time_window' + str(common_misc_params[1]) + 'flip_axis' + str(common_misc_params[2]) + '_' + common_misc_params[0]
    fit = common_misc_params[3]

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
                sizes, percentages = single_avalanche_plot(get_param_string(sz, common_general_params[0], common_general_params[1], common_general_params[2], common_general_params[3], common_general_params[4], common_general_params[5], common_general_params[7]), characterizing_text, fit)
            else:
                sizes, percentages = single_avalanche_plot(get_param_string(sz, common_general_params[0], common_general_params[1], common_general_params[2], common_general_params[3], common_general_params[4], common_general_params[5], common_general_params[7]), characterizing_text, False)
        else:
            sizes, percentages = single_avalanche_plot(get_param_string(sz, common_general_params[0], common_general_params[1], common_general_params[2], common_general_params[3], common_general_params[4], common_general_params[5], common_general_params[7]), characterizing_text, fit)
        list_of_sizes.append(sizes)
        list_of_percentages.append(percentages)

        plt.scatter(sizes, percentages, label=get_sz_string(sz))

    plt.ylim(1e-6, 1)
    plt.xlabel('Avalanche Size (s)')
    plt.ylabel('P(s)')
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(common_general_params[7] + '/' + common_general_string + '/' + characterizing_text + '_all_avalanches.jpg')
    plt.clf()

    #Perform scale invariance analysis
    if finite_size[0] == True:
        ###plt.rcParams.update({'legend.loc': 'upper left'})
        for i in range(len(size_list)):
            finite_size_x = [element/((list_of_arraysizes[i][0]*list_of_arraysizes[i][1])**finite_size[2]) for element in list_of_sizes[i]]
            finite_size_y = [a*b for a,b in zip([element**finite_size[1] for element in list_of_sizes[i]],list_of_percentages[i])]

            plt.scatter(finite_size_x, finite_size_y, label=get_sz_string(size_list[i]))

        plt.xlabel(r'~$s/L^{%0.2f}$' % finite_size[2])
        plt.ylabel(r'~$s^{-%0.2f}$P(s)' % finite_size[1])
        plt.legend()
        plt.xscale('log')
        plt.yscale('log')
        plt.savefig(common_general_params[7] + '/' + common_general_string + '/' + characterizing_text + '_finite_size.jpg')
        plt.clf()
