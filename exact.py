import numpy as np
import math
import sys

I = np.matrix('1 0; 0 1') 
X = np.matrix('0 1; 1 0')
Y = np.matrix('0 -1;1 0')*complex(0,1)
Z = np.matrix('1 0; 0 -1')
paulis = [I,X,Y,Z]

f_results = []
basis_dict = {}
dbg = 0

# apply e^{i\gamma Z_a Z_b Z_c} to equation containing variables a,b,c
def apply_C(num_vars, rho, overlap_eqns, gamma):
    for i in range(0, len(overlap_eqns)):
        solution = overlap_eqns[i][3]
        dabc = 1 - 2*solution # 0 --> 1 and 1 --> -1
        state = [I]*num_vars
        for j in range(0, 3):
            state[overlap_eqns[i][j]] = Z
        exp = kron(state)
        eiC = np.asmatrix(np.diag(np.exp(np.diag(np.complex(0,1)*0.5*gamma*dabc*exp))))
        rho =  eiC*rho*np.conj(eiC)
    return rho

# take the kronecker (tensor) product of a list of len(m) matrices
def kron(m):
    if len(m) == 1:
        return m
    total = np.kron(m[0], m[1])
    if len(m) > 2:
        for i in range(2, len(m)):
            total = np.kron(total, m[i])
    return total

def init_Y(num_vars, picked_eqn):
    input_state = [I]*num_vars
    for i in range(0, 3):
        input_state[picked_eqn[i]] = Y
    rho = kron(input_state)
    return rho

# given initially picked equation, 
def e3lin2_exact_helper(num_vars, picked_eqn, overlap_eqns, gamma):
    beta = math.pi/4

    # create our input state |+>|+>...|+>
    f_state = [0.5*(I+X)]*num_vars
    f_rho = kron(f_state)

    # create our input state I...IYYYI...I
    init_state = init_Y(num_vars, picked_eqn)
    op2 = apply_C(num_vars, init_state, overlap_eqns, gamma)
    expec = np.asscalar((op2*f_rho).trace()).real
    return expec

#Returns the set of equations with any shared variables
def shared_eqns(eqns, i):
    eqn = set(eqns[i][0:3])
    num_eqns = len(eqns)
    indices = []
    for j in range(num_eqns):
        if j == i:
            pass
        eqn_j = set(eqns[j][0:3])
        if len(eqn & eqn_j) != 0:
            indices.append(j)

    return np.array([eqns[j] for j in indices])

def e3lin2_exact(num_vars, num_eqns, i, eqns_location, gamma):
  
    eqns = dict()

    f = open(eqns_location, "r")
    for j in range(num_eqns):
        eqns[j] = map(int, f.readline().split(","))

    return e3lin2_exact_helper(num_vars, eqns[i], shared_eqns(eqns, i), gamma)


if __name__ == '__main__':
    if(len(sys.argv) != 6):
        print('Please specify <num_vars>, <num_eqns>, <eqn_number>, <eqns_location>, <gamma>')
    num_vars = int(sys.argv[1])
    num_eqns = int(sys.argv[2])
    eqn_number = int(sys.argv[3])
    eqns_location = sys.argv[4]
    gamma = float(sys.argv[5])
    print(e3lin2_exact(num_vars, num_eqns, eqn_number, eqns_location, gamma))
