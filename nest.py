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

# calculate the coefficient of <x|w>
def inner_x_w(x, all_eqns, gamma):
    n = int(math.log(x,2))
    bit_pos = {}
    # get positions of 1's in bitstring of x
    for i in range(n, -1, -1):
        if (x - 2**i) > 0:
            bit_pos[n-i] = n-i
            x = x - 2**i
        elif (x - 2**i) == 0:
            bit_pos[n-i] = n-i
            break

    total_phase = 1/math.sqrt(2**n)
    # for each equation, calculate overlap
    for i in range(len(all_eqns)): 
        num_ones = 0
        for j in range(0,3):
            if all_eqns[i][j] in bit_pos: num_ones += 1
        total_phase *= math.e**(-complex(0,1)*gamma*0.5*((-1)**(all_eqns[i][3] + num_ones)))

    return total_phase

# get alpha (entry) and row index (of first nonzero entry in column x)
def get_entry(x, tensor_list, size, level, alpha, row_index):

    if size == 1:
        print('alpha: '),
        print(alpha),
        print(' row_index: %d' % row_index)
        return 0

    if x < size/2:
        # left side of sub matrix
        if tensor_list[level] == 1:
            alpha *= complex(0,1)
            row_index += size/2
        level += 1
        size /= 2
        get_entry(x, tensor_list, size, level, alpha, row_index)
    else: 
        # right side of sub matrix
        if tensor_list[level] == 1:
            alpha *= -complex(0,1)
        else: 
            row_index += size/2
        level += 1
        x -= size/2
        size /= 2
        get_entry(x, tensor_list, size, level, alpha, row_index)

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
    base_eqn = all_eqns[i]

    # 0 is I, 1 is Y
    ten = [1,0,0]
    get_entry(0, ten, 2**(len(ten)), 0, 1, 0)
    get_entry(1, ten, 2**(len(ten)), 0, 1, 0)
    get_entry(2, ten, 2**(len(ten)), 0, 1, 0)
    get_entry(3, ten, 2**(len(ten)), 0, 1, 0)
    get_entry(4, ten, 2**(len(ten)), 0, 1, 0)
    get_entry(5, ten, 2**(len(ten)), 0, 1, 0)
    get_entry(6, ten, 2**(len(ten)), 0, 1, 0)
    get_entry(7, ten, 2**(len(ten)), 0, 1, 0)

    #print(inner_x_w(23, all_eqns, gamma))
    #return e3lin2_exact_helper(num_vars, base_eqn, neighbor_eqns, gamma)

if __name__ == '__main__':
    if(len(sys.argv) != 4):
        print('Please specify <eqn_number>, <eqns_location>, <gamma>')
    else:
        eqn_number = int(sys.argv[1])
        eqns_location = sys.argv[2]
        gamma = float(sys.argv[3])

        print(e3lin2_exact(eqn_number, eqns_location, gamma))
