import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys
import math

dbg = 0

def loadfile(filename):
    with open (filename, 'rb') as fp:
        total_list = pickle.load(fp)
    return total_list

def get_error(n, delta, gammas):
    errors = []    
    for gam in gammas:
        D = math.fabs(math.cos(gam)) + math.fabs(math.sin(gam))
        #err = math.sqrt(4*D*D*math.log(2/delta)/(2*n))
        err = math.sqrt(4*D*D*math.log(1/delta)/(1*n))
        errors.append(err)
    return errors
    
def graph(info, data):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.axis([0, 5, 0, 35])
    
    if dbg: print(data)
    # info = [num_vars, num_eqns, num_samples]
    # data = [gamma, beta, <C>]
    confidence = 0.99
    delta = 1 - confidence
    for i in range(0, len(data)):
        error_bars = get_error(info[3], delta,  data[i][0])
        ax.plot(data[i][0], data[i][2], '-o')
        ax.errorbar(data[i][0], data[i][2], yerr=error_bars, fmt='-o')

    plt.title('beta = %.2f, num_vars = %d, d = %d, num_eqns = %d, num_samples = %d' % (data[0][1][0], info[0], info[1], info[2], info[3]))
    ax.set_ylabel('<C>', fontsize = 16)
    ax.set_xlabel('gamma', fontsize = 16)
    plt.show()

if __name__ == '__main__':
    total_list = loadfile(sys.argv[1])
    print(total_list)
    info = total_list[0][0]
    eqns = total_list[0][1]
    data = []
    for i in range(0, 2):# len(total_list)):
        data.append(total_list[i][2])
    print(eqns)
    graph(info, data)
