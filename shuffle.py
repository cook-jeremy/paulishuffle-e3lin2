import numpy as np
import math
import itertools
from scipy.linalg import expm

I = np.matrix('1 0; 0 1')
X = np.matrix('0 1; 1 0')
Y = np.matrix('0 -1;1 0')*complex(0,1)
Z = np.matrix('1 0; 0 -1')
paulis = [I,X,Y,Z]

H = np.matrix('1 1; 1 -1')*(1/math.sqrt(2))

# take the kronecker (tensor) product of a list of len(m) matrices
def kron(m):
    total = np.kron(m[0], m[1])
    if len(m) > 2:
        for i in range(2, len(m)):
            total = np.kron(total, m[i])
    return total

# decompose a 2^n x 2^n matrix into n-qubit tensor products of Pauli matrices
def decompose(R):
    # take the cartesian product of the pauli matrices n times and tensor them together
    dim = int(np.log2(R.shape[0]))
    cart_prod = itertools.product(paulis, repeat=dim)
    basis = []
    if dim > 1:
        for prod in cart_prod:
            basis.append(kron(prod))
    else:
        basis = paulis
   
    # decompose matrix R into its components
    const = 1/(math.pow(2,dim))
    c = []
    for B in basis: 
        c.append(np.asscalar((const*R*B).trace()))
    return c

# from a state R, pick a pauli matrix from the probability distribution
def pick_pauli(R):
    # decompose the matrix into sum of paulis and get dimension
    consts = decompose(R)
    dim = int(np.log2(R.shape[0]))

    # find the normalized weights
    consts_sum = 0
    for i in consts: 
        consts_sum += abs(i)
    weights = []
    for i in consts: 
        weights.append(abs(i)/consts_sum)

    # pick a random weight
    choice = np.random.choice(int(math.pow(4,dim)), p=weights)

    # pick a pauli from the weights
    P_CHOICE = 0
    for i in range(0,dim):
        pos = int((choice/math.pow(4,i)) % 4)
        if i == 0:
            P_CHOICE = paulis[pos]
        else:
            P_CHOICE = np.kron(paulis[pos],P_CHOICE)
    
    return [P_CHOICE, consts[choice], weights[choice]]

def classical_circuit(INIT):
    result = (Z*X*H*INIT*H*X).trace()
    print(result)

def circuit(INIT):
    results = []
    samples = 10000
    for i in range(0,samples):
        # pick a pauli for INIT
        r1 = pick_pauli(INIT)
        STATE = r1[0]

        STATE = H*STATE*H
        r2 = pick_pauli(STATE)
        STATE = r2[0]

        STATE = X*STATE*X
        r3 = pick_pauli(STATE)
        STATE = r3[0]
        
        rf = r1[1]*r2[1]*r3[1]
        pf = r1[2]*r2[2]*r3[2]
        
        p_hat = (rf/pf)*(STATE*Z).trace()
        results.append(p_hat)

    avg = sum(results)/samples
    print(np.asscalar(avg))

if __name__ == '__main__':
    # input state is a\ket{0} + b\ket{1}
    a = 1/2
    b = math.sqrt(3)/2
    #a = 1
    #b = 0
    INIT = np.matrix([[a*a,a*b],[b*a,b*b]])
    #classical_circuit(INIT)
    #circuit(INIT)
    pick_pauli(np.kron(X,Y))
