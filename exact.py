import numpy as np
import math, sys, random

def inner_x_w(x, num_vars, bit_eqns, gamma):
    total_phase = 1/math.sqrt(2**num_vars)
    # for each equation, calculate overlap
    for i in range(len(bit_eqns)): 
        total_phase *= math.e**(-complex(0,1)*gamma*0.5*((-1)**(bit_eqns[i][1] + num_ones(x & bit_eqns[i][0]))))
    return total_phase

def convert_to_bin(num_vars, eqn_set):
    n = num_vars - 1
    bit_eqns = []
    for eqn in eqn_set:
        y = 0
        for j in range(3):
            y += 2**(n - eqn[j])
        bit_eqns.append([y, eqn[3]])
    return bit_eqns

def num_ones(n): 
    count = 0
    while (n): 
        count += n & 1
        n >>= 1 
    return count

def e3lin2_exact(i, eqns_location, gamma):
    all_eqns = []
    f = open(eqns_location, "r")
    line = f.readline()
    while line:
        all_eqns.append(list(map(int, line.split(","))))
        line = f.readline()

    base_eqn = all_eqns[i]
    qubits = set()
    for j in range(len(all_eqns)):
        qubits |= set(all_eqns[j][:-1])
    num_vars = len(qubits)
    bit_eqns = convert_to_bin(num_vars, all_eqns)
    ys = bit_eqns[i][0]
    
    total = 0
    num_samples = 1000
    for i in range(num_samples):
        x = random.randint(0, 2**(num_vars)-1)
        alpha = -complex(0,1)*((-1)**(num_ones(x & ys)))
        row_index = x ^ ys
        inner1 = np.conj(inner_x_w(row_index, num_vars, bit_eqns, gamma))
        inner2 = inner_x_w(x, num_vars, bit_eqns, gamma)
        total += 2**num_vars*inner1*inner2*alpha
    total /= num_samples

    '''
    # exact calculation, sum over basis states instead of sampling
    total = 0
    for i in range(2**num_vars):
        alpha = -complex(0,1)*((-1)**(num_ones(i & ys)))
        row_index = i ^ ys
        inner1 = np.conj(inner_x_w(row_index, num_vars, bit_eqns, gamma))
        inner2 = inner_x_w(i, num_vars, bit_eqns, gamma)
        total += 2**num_vars*inner1*inner2*alpha

    total /= 2**num_vars    
    #total *= (-1)**base_eqn[3] why don't we need this?
    '''

    # error calculation
    delta = 0.001
    eps = math.sqrt(2*math.log(2/delta)/num_samples)
    error = eps*len(all_eqns)/2

    return np.real(total), error

if __name__ == '__main__':
    if(len(sys.argv) != 4):
        print('Please specify <eqn_number>, <eqns_location>, <gamma>')
    else:
        eqn_number = int(sys.argv[1])
        eqns_location = sys.argv[2]
        gamma = float(sys.argv[3])

        print(e3lin2_exact(eqn_number, eqns_location, gamma))
