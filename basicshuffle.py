import numpy as np
import math

I = np.matrix('1 0; 0 1')
X = np.matrix('0 1; 1 0')
Y = np.matrix('0 -1;1 0')*complex(0,1)
Z = np.matrix('1 0; 0 -1')
paulis = [I,X,Y,Z]

# decompose for n=1: \rho = \sigma_j tr(\rho\sigma_j)
def decompose_1d(R):
    c = []
    for P in paulis:
        c.append(np.asscalar((0.5*R*P).trace()))
    return c
'''
def pick_pauli(R):
    # decompose the matrix into sum of paulis
    consts = decompose_1d(R)

    # find the normalized weights
    consts_sum = 0
    for i in consts: consts_sum += abs(i)
    weights = []
    for i in consts: weights.append(abs(i)/consts_sum)

    # pick a specific pauli based on the weights
    return paulis[np.random.choice(4,p=weights)]
'''

def pick_pauli(R):
    # decompose the matrix into sum of paulis
    consts = decompose_1d(R)

    # find the normalized weights
    consts_sum = 0
    for i in consts: consts_sum += abs(i)
    #print(consts)
    weights = []
    for i in consts: weights.append(abs(i)/consts_sum)
    #print(weights)

    # pick a specific pauli based on the weights
    choice = np.random.choice(4,p=weights)
    return [paulis[choice],consts[choice]/weights[choice]]

def classical_circuit(INIT):
    #R2 = X*INIT*X
    R2 = INIT
    result = (Z*R2).trace()
    print(result)

def circuit(INIT):
    #R2 = pick_pauli(INIT)
    results = []
    samples = 5000
    for i in range(0,samples):
        r = pick_pauli_last(INIT)
        #print(r)
        p_hat = r[1]*(r[0]*Z).trace()
        #print(p_hat)
        results.append(p_hat)

    avg = sum(results)/samples
    print(np.asscalar(avg))


if __name__ == '__main__':
    # input state is a\ket{0} + b\ket{1}
    a = 1/2
    b = math.sqrt(3)/2
    INIT = np.matrix([[a*a,a*b],[b*a,b*b]])
    classical_circuit(INIT)
    circuit(INIT)
    #result = circuit(INIT)
