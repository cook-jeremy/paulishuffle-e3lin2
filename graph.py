import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys

dbg = 1

def loadfile(filename):
    with open (filename, 'rb') as fp:
        data = pickle.load(fp)
    return data
   
def graph(data):
    #if dbg: print(data)
    for i in range(len(data[1])):
        print(data[0][i], "\t", data[1][i], "\t", data[2][i], "\t", data[3][i])
    # data = [x1, y1, x2, y2] x1 and y1 are exact, x2 and y2 are GPU pauli shuffle
    plt.plot(data[0], data[1], '-o', label= 'exact')
    plt.plot(data[2], data[3], '-o', label= 'pauli')
    #plt.plot(data[0], data[2], '-o', label = 'D')
    #plt.plot(data[0], data[3], '-o', label = 'X')
    plt.title('Exact & Pauli')
    plt.legend()
    plt.ylabel('<C>')
    plt.xlabel('gamma')
    plt.show()

if __name__ == '__main__':
    graph(loadfile(sys.argv[1]))
