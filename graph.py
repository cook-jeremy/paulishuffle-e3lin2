import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys
import math

dbg = 0

def loadfile(filename):
    with open (filename, 'rb') as fp:
        result_list = pickle.load(fp)
    return result_list

def get_error(n, delta, gammas):
    errors = []    
    for gam in gammas:
        D = math.fabs(math.cos(gam)) + math.fabs(math.sin(gam))
        err = math.sqrt(4*D*D*math.log(2/delta)/(2*n))
        errors.append(err)
    return errors
    
def graph(info, data):
    if dbg: print(data)
    # info = [num_vars, num_eqns, num_samples]
    # data = [gamma, beta, <C>]
    confidence = 0.99
    delta = 1 - confidence
    error_bars = get_error(info[3], delta,  data[0])
    plt.errorbar(data[0], data[2], yerr=error_bars, fmt='ro')
    #plt.plot(data[0], data[2])
    plt.title('beta = %.2f, num_vars = %d, d = %d, num_eqns = %d, num_samples = %d' % (data[1][0], info[0], info[1], info[2], info[3]))
    plt.ylabel('<C>')
    plt.xlabel('gamma')
    plt.show()

if __name__ == '__main__':
    result_list = loadfile(sys.argv[1])
    info = result_list[0]
    eqns = result_list[1]
    data = result_list[2]
    print(eqns)
    graph(info, data)
