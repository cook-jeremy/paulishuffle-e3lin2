import numpy as np
import math

I = np.matrix('1 0; 0 1')
X = np.matrix('0 1; 1 0')
Y = np.matrix('0 -1;1 0')*complex(0,1)
Z = np.matrix('1 0; 0 -1')
paulis = [I,X,Y,Z]

H = np.matrix('1 1; 1 -1')*(1/math.sqrt(2))

# decompose for n=1: \rho = \sigma_j tr(\rho\sigma_j)
def decompose(R):
    c = []
    for P in paulis:
        c.append(np.asscalar((0.5*R*P).trace()))
    return c

def pick_pauli(R):
    # decompose the matrix into sum of paulis
    consts = decompose(R)
    #print(consts)

    # find the normalized weights
    consts_sum = 0
    for i in consts: consts_sum += abs(i)
    #print(consts)
    weights = []
    for i in consts: weights.append(abs(i)/consts_sum)
    #print(weights)

    # pick a specific pauli based on the weights
    choice = np.random.choice(4,p=weights)
    return [paulis[choice],consts[choice],weights[choice]]

def classical_circuit(INIT):
    result = (Z*H*INIT*H).trace()
    print(result)

def circuit(INIT):
    results = []
    samples = 10000
    for i in range(0,samples):
        # circuit: |0> -->  <Z>

        # pick a pauli for INIT
        r1 = pick_pauli(INIT)
        STATE = r1[0]

        STATE = H*STATE*H
        r2 = pick_pauli(STATE)
        STATE = r2[0]
        
        rf = r1[1]*r2[1]
        pf = r1[2]*r2[2]
        
        p_hat = (rf/pf)*(STATE*Z).trace()
        results.append(p_hat)

    avg = sum(results)/samples
    print(np.asscalar(avg))

def circuit2(INIT):
    results = []
    samples = 1
    for i in range(0,samples):
        # circuit: |0> -->  H  --> H  -->  <Z>

        # pick a pauli for H
        r1 = pick_pauli(H)
        #print('r1[1] = ' + str(r1[1]))
        # apply it to our input state
        print(r1[0])
        R1 = r1[0]*INIT*r1[0]

        # pick again for H
        r2 = pick_pauli(H)
        # apply it to the current state
        R2 = r2[0]*R1*r2[0]
        print(r2[0])

        rf = r1[1]*r2[1]
        print('rf = ' + str(rf))
        pf = r1[2]*r2[2]
        print('pf = ' + str(pf))

        #print(r)
        p_hat = (rf/pf)*(R2*Z).trace()
        print('p_hat = ' + str(p_hat))
        print('R2 = ' + str(R2))
        #print(p_hat)
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
    #print(decompose(INIT))
    classical_circuit(INIT)
    circuit(INIT)
    #result = circuit(INIT)
