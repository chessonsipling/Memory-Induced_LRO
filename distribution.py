import math
import numpy as np

#For fitting to a gaussian distribution
def gaussian(x, mu, sigma, A):
	return A / (sigma * np.sqrt(2 * math.pi)) * np.exp((-1/2) * ((x - mu)/sigma)**2)

#For fitting to a scale-free distribution
def scale_free(x, a, C):
	return C * x**a
