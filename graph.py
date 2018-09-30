import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys

dbg = 0

def loadfile(filename):
    with open (filename, 'rb') as fp:
        result_list = pickle.load(fp)
    return result_list
    
def graph(info, data):
    if dbg: print(data)
    plt.plot(data[0], data[2])
    plt.title('beta = %.2f, num_vars = %d, d = %d, num_eqns = %d' % (data[1][0], info[0], info[1], info[2]))
    plt.ylabel('<C>')
    plt.xlabel('gamma')
    plt.show()

if __name__ == '__main__':
    result_list = loadfile(sys.argv[1])
    info = result_list[0]
    eqns = result_list[1]
    data = result_list[2]

    graph(info, data)
