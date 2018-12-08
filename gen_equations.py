import random
import sys
import numpy as np
import numpy.linalg

dbg = 1

def fully_connected(EQNS):
    # given a set of eqns, ensure the graph where nodes are variables and edges are equations is fully connected
    num_rows = EQNS.shape[0]
    num_cols = EQNS.shape[1]

    # get first row variables
    pot_vars = [] # variables still needed to search
    seen_vars = [] # variables where we found all their equations
    for j in range(0, num_cols):
        if EQNS[0, j] == 1:
            pot_vars.append(j)     

    # while pot_vars is not empty and size of seen_vars is < n
    flag = True
    while flag:
        var = pot_vars[0]
        seen_vars.append(var)
        pot_vars.remove(var)

        for i in range(0, num_rows):
            if EQNS[i, var] == 1:
                for j in range(0, num_cols):
                    if EQNS[i, j] == 1 and j != var:
                        pot_vars.append(j)

        # remove duplicates 
        pot_vars = list(set(pot_vars))
        # remove those already in seen_vars
        for elem in seen_vars:
            if elem in pot_vars:
                pot_vars.remove(elem)
        # look at all variables that are connected by some equation
        merged_list = pot_vars + seen_vars

        if len(pot_vars) > 0 and len(merged_list) < num_cols:
            # not done yet searching
            continue
        elif len(merged_list) == num_cols:
            # fully connected because we've seen every variable
            return True
        elif len(pot_vars) == 0:
            # seen_vars is not full and we are out of things to search, so not connected
            flag = False
    return False

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

def solvable(EQNS, SOLS):
    # check if the system is inconsistent via Rouche Capelli theorem
    AUGMENTED = np.concatenate((EQNS, SOLS), 1)
    rank_AUG = sum(np.linalg.svd(AUGMENTED)[1] != 0)
    rank_COEFF = sum(np.linalg.svd(EQNS)[1] != 0)
    #rank_AUG = matrix_rank(AUGMENTED)
    #rank_COEFF = matrix_rank(EQNS)
    if rank_COEFF < rank_AUG:
        return False
    else:
        return True

def unique_eqns(EQNS, SOLS):
    # check if all the eqns are unique
    num_rows = EQNS.shape[0]
    A = np.concatenate((EQNS, SOLS), 1)

    for i in range(len(A)): #generate pairs
        for j in range(i+1,len(A)): 
            if np.array_equal(A[i],A[j]): #compare rows
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
    num_failures = 0
    while not satisfied:
        # create A and b in Ax = b
        sys = create_eqns(n, d, f)
        
        # check if we have used all variables
        if not contains_all_vars(sys[0]):
            #print('doesn\'t contain all variables, trying again...')
            num_failures += 1
            continue

        # check for the same equations
        if not unique_eqns(sys[0], sys[1]):
            num_failures += 1
            continue

        # check if we are fully connected
        if not fully_connected(sys[0]):
            #print('ins\'t fully connected, trying again...')
            num_failures += 1
            continue

        # row reduce and check if system is solvable
        if solvable(sys[0], sys[1]):
            num_failures += 1
            continue

        # passed all conditions
        break

    print('num failures: %s' % num_failures)
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
    
    if f <= n/3:
        print('ERROR(gen_equations): not enough equations to have fully connected graph')
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
    ret = gen_eqns(9, 100, 9)
    print(ret[0])
    print(ret[1])
