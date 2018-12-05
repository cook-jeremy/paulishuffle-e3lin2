import numpy as np
import math
from scipy.linalg import expm
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
        eiC = np.asmatrix(expm(np.complex(0,1)*0.5*gamma*dabc*exp))
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
def e3lin2_exact(num_vars, picked_eqn, overlap_eqns, gamma):
    beta = math.pi/4

    # create our input state |+>|+>...|+>
    f_state = [0.5*(I+X)]*num_vars
    f_rho = kron(f_state)

    # create our input state I...IYYYI...I
    init_state = init_Y(num_vars, picked_eqn)
    op2 = apply_C(num_vars, init_state, overlap_eqns, gamma)
    expec = np.asscalar((op2*f_rho).trace()).real
    print(expec)

if __name__ == '__main__':
    if(len(sys.argv) != 5):
        print('Please specify <num_vars>, <picked_eqn>, <overlapping_eqns>, <gamma>')
    #picked_eqn = [0,1,2,0]
    #overlapping_eqns = [[0,3,4,1],[0,1,2,0]]
    #e3lin2_exact(5, picked_eqn, overlapping_eqns, 1.3)
    num_vars = int(sys.argv[1])
    picked_eqn = np.array(eval(sys.argv[2]))
    overlap_eqns = np.array(eval(sys.argv[3]))
    gamma = float(sys.argv[4])
    e3lin2_exact(num_vars, picked_eqn, overlap_eqns, gamma)
