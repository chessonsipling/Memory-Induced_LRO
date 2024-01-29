from get_string import get_param_string

#Unwraps general parameters to make function arguments a bit more succinct
def param_unwrapper(general_params):

    sz = general_params[0]
    gamma = general_params[1]
    delta = general_params[2]
    zeta = general_params[3]
    xini = general_params[4]
    T = general_params[5]
    dt = general_params[6]
    transient = general_params[7]
    subdirec = general_params[8]
    param_string = get_param_string(sz, gamma, delta, zeta, xini, T, dt, subdirec)

    return sz, gamma, delta, zeta, xini, T, dt, transient, param_string