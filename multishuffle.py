import numpy as np
import math
import itertools
from scipy.linalg import expm
import gen_equations
import datetime
import os
import pickle

I = np.matrix('1 0; 0 1') 
X = np.matrix('0 1; 1 0')
Y = np.matrix('0 -1;1 0')*complex(0,1)
Z = np.matrix('1 0; 0 -1')
paulis = [I,X,Y,Z]

f_results = []
dbg = 0

# apply e^{i\beta X} to each qubit
def apply_B(beta, states):
    eiB = np.matrix([[math.cos(beta),complex(0,math.sin(beta))],[complex(0,math.sin(beta)),math.cos(beta)]])
    result = [] 
    for state in states[0]:
        # apply the B operator locally
        decomp = pick_pauli(eiB*state*np.conj(eiB))
        result.append(decomp[0][0])
        states[1] *= decomp[1]
        states[2] *= decomp[2]
    return [result, states[1], states[2]]

# apply e^{i\gamma Z_a Z_b Z_c} to equation containing variables a,b,c
def apply_C(num_vars, equations, gamma, states):
    for i in range(0, len(equations)):
        index0 = equations[i][0]
        index1 = equations[i][1]
        index2 = equations[i][2]
        solution = equations[i][3]
        dabc = 1 - 2*solution # 0 --> 1 and 1 --> -1
        local_state = [states[0][index0], states[0][index1], states[0][index2]]
        local_state = kron(local_state)
        eiC = np.asmatrix(expm(complex(0,1)*gamma*dabc*kron([Z,Z,Z])))
        decomp = pick_pauli(eiC*local_state*np.conj(eiC))
        states[1] *= decomp[1]
        states[2] *= decomp[2]
        states[0][index0] = decomp[0][0]
        states[0][index1] = decomp[0][1]
        states[0][index2] = decomp[0][2]
    return states

# create the C observable
def create_C(num_vars, equations):
    # if the equation x_a + x_b + x_c exists then add Z_a * Z_b * Z_c to our observable
    C = 0
    for i in range(0, len(equations)):
        local_list = [I]*num_vars
        solution = equations[i][3]
        dabc = 1 - 2*solution # 0 --> 1 and 1 --> -1
        for j in range(0,3):
            local_list[equations[i][j]] = Z
        L = dabc*kron(local_list)
        C += L

    C = (1/2)*C
    # return the final observable
    return C

# take the kronecker (tensor) product of a list of len(m) matrices
def kron(m):
    if len(m) == 1:
        return m
    total = np.kron(m[0], m[1])
    if len(m) > 2:
        for i in range(2, len(m)):
            total = np.kron(total, m[i])
    return total

# decompose a 2^n x 2^n matrix into n-qubit tensor products of Pauli matrices
def decompose(R):
    # take the cartesian product of the pauli matrices n times and tensor them together
    R = np.asmatrix(R)
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
    R = np.asmatrix(R)
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

    # pick a pauli for each qubit from the weights
    pauli_choices = []
    for i in range(dim-1,-1,-1):
        pos = int((choice/math.pow(4,i)) % 4)
        pauli_choices.append(paulis[pos])

    return [pauli_choices, consts[choice], weights[choice]]

# where the magic happens
def e3lin2(num_vars, num_samples, num_steps, input_eqns):
    # create the observable C
    C = create_C(num_vars, input_eqns)

    # create our input state |+>|+>...|+>
    input_state = [0.5*(I+X)]*num_vars
    rho = kron(input_state)
    
    gammas = []
    betas = []
    expectations = []

    gamma = 0
    beta = math.pi/4
    for g in range(0, num_steps):
        results = []
        for i in range(0, num_samples):
            # pick pauli, pass into B then C
            init = pick_pauli(C)
            op1 = apply_B(beta, init)
            op2 = apply_C(num_vars, input_eqns, gamma, op1)
            
            #print(kron(op2[0]))
            p_hat = (op2[1]/op2[2])*(kron(op2[0])*rho).trace()
            results.append(p_hat)

        avg = sum(results) / num_samples
        expectation = np.real(np.asscalar(avg))
        print('gamma = %.2f, beta = %.2f, <C> = %.4f' % (gamma, beta, expectation))
        gammas.append(gamma)
        betas.append(beta)
        expectations.append(expectation)
        gamma += math.pi/num_steps

    f_results.append([gammas, betas, expectations])

def convert_eqns(input_eqns):
    EQNS = input_eqns[0]
    b = input_eqns[1]
    num_rows = EQNS.shape[0]
    num_cols = EQNS.shape[1]
    eqn_list = []
    for i in range(0, num_rows):
        eqn = []
        for j in range(0, num_cols):
            if EQNS[i, j] == 1:
                eqn.append(j)
        eqn.append(int(np.asscalar(b[i])))
        eqn_list.append(eqn)
    return eqn_list


if __name__ == '__main__':
    num_vars = 5
    d_constraint = 9
    num_eqns = 10
    num_samples = 300
    num_steps = 50
    total_results = []

    input_eqns = gen_equations.gen_eqns(num_vars, d_constraint, num_eqns)
    for i in range(0,5):
        f_results = []
        # get input matrix and result, Ax = b
        f_results.append([num_vars, d_constraint, num_eqns, num_samples])
        f_results.append(input_eqns)
        print('num_vars: %d, d_constraint: %d, num_eqns: %d, num_samples: %d' % (num_vars, d_constraint, num_eqns, num_samples))
        #print('equations: \n %s x = %s' % (str(input_eqns[0]), str(input_eqns[1])))
        print(np.c_[input_eqns[0], input_eqns[1]])
        # turn input eqns matrix into a list
        input_eqns_list = convert_eqns(input_eqns)
        e3lin2(num_vars, num_samples, num_steps, input_eqns_list)
        total_results.append(f_results)
        
    # print results to file
    now = datetime.datetime.now()
    directory = 'results/' + str(now.strftime('%m-%d-%Y')) + '/'
    if not os.path.exists(directory):
            os.makedirs(directory)
    filename = directory + str(now.strftime('%H:%M:%S%p'))

    with open(filename, 'wb') as fp:
        pickle.dump(total_results, fp)
