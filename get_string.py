#Returns characteristic string corresponding to this set of parameters
def get_param_string(sz, gamma, delta, zeta, xini, T, dt, sub_directory):
    return sub_directory + '/sz' + str(sz[0]) + 'x' + str(sz[1]) + '_g' + str(gamma) + '_d' + str(delta) + '_z' + str(zeta) + '_x' + str(xini) + '_T' + str(T) + '_dt' + str(dt)

#Returns lattice size information in a slightly more convenient way for plotting in figures
def get_sz_string(sz):
    return str(sz[0]) + 'x' + str(sz[1])
