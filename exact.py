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

        exp = np.array(kron(state)).diagonal()
        eiC = np.asmatrix(np.diag(np.exp(np.complex(0,1)*0.5*gamma*dabc*exp)))
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
    # create our input state |+>|+>...|+>
    f_state = [0.5*(I+X)]*num_vars
    f_rho = kron(f_state)

    # create our input state I...IYYYI...I
    init_state = init_Y(num_vars, picked_eqn)
    op2 = apply_C(num_vars, init_state, overlap_eqns, gamma)
    expec = np.asscalar((op2*f_rho).trace()).real
    return expec


def e3lin2_exact(i, eqns_location, gamma):
    all_eqns = []
    f = open(eqns_location, "r")
    line = f.readline()
    while line:
        all_eqns.append(map(int, line.split(",")))
        line = f.readline()

    base = set(all_eqns[i][:-1])  # qubits of base equation
    qubits = set(all_eqns[i][:-1])  # qubits in subcircuit
    neighbor_eqns = []

    for j in range(len(all_eqns)):
        if len(base & set(all_eqns[j][:-1])) > 0:
            qubits |= set(all_eqns[j][:-1])
            neighbor_eqns.append(all_eqns[j])

    qubits = list(qubits)
    num_vars = len(qubits)

    base_eqn = all_eqns[i]
    for k in range(3):
        base_eqn[k] = qubits.index(base_eqn[k])

    for j in range(len(neighbor_eqns)):
        for k in range(3):
            neighbor_eqns[j][k] = qubits.index(neighbor_eqns[j][k])

    print(num_vars)
    print(base_eqn)
    print(neighbor_eqns)
    return e3lin2_exact_helper(num_vars, base_eqn, neighbor_eqns, gamma)


if __name__ == '__main__':
    if(len(sys.argv) != 4):
        print('Please specify <eqn_number>, <eqns_location>, <gamma>')
    else:
        eqn_number = int(sys.argv[1])
        eqns_location = sys.argv[2]
        gamma = float(sys.argv[3])

        print(e3lin2_exact(eqn_number, eqns_location, gamma))
