import random
import sys
import numpy as np

dbg = 1

def fully_connected(eqns):
    # given a set of eqns, ensure the graph where nodes are variables and edges are equations is fully connected
    return True

def contains_all_vars(EQNS):
    # check if all columns have at least one 1 in them
    num_rows = EQNS.shape[0]
    num_cols = EQNS.shape[1]

    for j in range(0, num_cols):
        flag = False
        for i in range(0, num_rows):
            if EQNS[i, j] == 1:
                flag = True
        if not flag:
            # found a column without a 1 in it
            return False

    # otherwise
    return True

def init_eqn_select(n, d, d_constraint):
    # given how many times a variable has already shown up (d_constraint), return list of useable variables
    eqn_select = []
    for i in range(0,n):
        # if we haven't used d of this variable, add it to the list of potential eqns to be picked
        if d_constraint[i] < d:
            eqn_select.append(i)
    return eqn_select

def gen_eqns(n, d, f):
    satisfied = False
    while not satisfied:
        # create A and b in Ax = b
        sys = create_eqns(n, d, f)
        
        # check if we have used all variables
        if not contains_all_vars(sys[0]):
            continue

        # check if we are fully connected
        if not fully_connected(sys[0]):
            continue

        # check if we have determinant 0
        # row reduce and check if system is solvable
        return sys
        

def create_eqns(n, d, f):
    EQNS = np.asmatrix(np.zeros((f, n)))
    SOLS = np.asmatrix(np.zeros((f, 1)))

    if f < n/3:
        print('ERROR(gen_equations): too few equations to use all variables!')
        sys.exit(1)

    if 3*f > n*d:
        print('ERROR(gen_equations): 3f <= nd')
        sys.exit(1)

    d_constraint = [0]*n

    # loop over equations left to create 
    while f > 0:
        # get variables to be used [1,2,...,n]
        eqn_select = init_eqn_select(n, d, d_constraint)
        if(len(eqn_select) < 3): 
            print('ERROR(gen_equations): too many constraints (d is too small)')
            print('d_constraint = ' + str(d_constraint))
            sys.exit(1)

        for i in range(0,3):
            # pick a random variable
            s = random.randint(0,len(eqn_select)-1)
            var = eqn_select[s]

            # add to matrix
            EQNS[f-1, var] = 1
            d_constraint[var] += 1
            eqn_select.remove(var)

        # pick solution
        SOLS[f-1, 0] = random.randint(0,1)
        f = f - 1

    # make the equation set degenerate so we dont have exact solutions (is there a better way to do this?)
    return [EQNS, SOLS]

if __name__ == '__main__':
    ret = gen_eqns(6, 10, 10)
    print(ret[0])
    print(ret[1])
