#Returns characteristic string corresponding to this set of parameters
def get_param_string(sz, gamma, delta, zeta, xini, T, dt, sub_direct):
    d = len(sz)
    if d == 1:
        return sub_direct + '/sz' + str(sz[0]) + '_g' + str(gamma) + '_d' + str(delta) + '_z' + str(zeta) + '_x' + str(xini) + '_T' + str(T) + '_dt' + str(dt)
    elif d == 2:
        return sub_direct + '/sz' + str(sz[0]) + 'x' + str(sz[1]) + '_g' + str(gamma) + '_d' + str(delta) + '_z' + str(zeta) + '_x' + str(xini) + '_T' + str(T) + '_dt' + str(dt)

#Returns lattice size information in a slightly more convenient way for plotting in figures
def get_sz_string(sz):
    d = len(sz)
    if d == 1:
        return str(sz[0])
    elif d == 2:
        return str(sz[0]) + 'x' + str(sz[1])
