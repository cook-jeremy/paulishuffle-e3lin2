import random
import sys

def init_eqn_select(n, d, d_constraint):
    # init random equation selecter
    eqn_select = []

    all_vars = True
    for i in range(0,n):
        # haven't used all variables yet
        if d_constraint[i] == 0:
            all_vars = False
            eqn_select.append(i)

    # create a list with all variables that haven't been used yet
    if not all_vars:
        while len(eqn_select) < 3:
            flag = False
            add = random.randint(0,n-1)
            for i in range(0,len(eqn_select)):
                if eqn_select[i] == add:
                    flag = True
            if not flag:
                eqn_select.append(add)
    else: 
        for i in range(0,n):
            # if we haven't used d of this variable, add it to the list of potential eqns to be picked
            if d_constraint[i] < d:
                eqn_select.append(i)
    return eqn_select

# given n (num of vars), d (max eqs for a single var), g (num of equations), create a system of linear equations that uses exactly n vars in g eqns with NO exact solutions
def create_eqn_list(n, d, f):
    equations = []

    if f < n/3:
        print('ERROR too few equations to use all variables!')
        sys.exit(1)

    # only make f - 1 equations for degeneracy later
    '''
    f = 1
    if g > 1:
        f = g - 1
    '''

    d_constraint = [0]*n
    #print('d_constra = ' + str(d_constraint) + '\n')

    # loop over equations left to create 
    while f > 0:
        # create init eqn selecter [1,2,...,n]
        eqn_select = init_eqn_select(n, d, d_constraint)
        if(len(eqn_select) < 3): 
            print('ERROR too many constraints (d is too small)')
            print('d_constraint = ' + str(d_constraint))
            sys.exit(1)

        #print('eqn start = ' + str(eqn_select))
        local_eqn = []

        for i in range(0,3):
            s = random.randint(0,len(eqn_select)-1)
            var = eqn_select[s]
            local_eqn.append(var)
            d_constraint[var] += 1
            eqn_select.remove(var)

        local_eqn.sort()
        # pick solution
        local_eqn.append(random.randint(0,1))

        #print('local_eqn = ' + str(local_eqn))
        #print('d_constra = ' + str(d_constraint))
        #print('')
        equations.append(local_eqn)
        f = f - 1

    # make the equation set degenerate so we dont have exact solutions (is there a better way to do this?)
    '''
    if g > 1:
        opposite = 0
        if equations[0][3] == 0:
            opposite = 1
        degen_eqn = [equations[0][0],equations[0][1],equations[0][2],opposite]
        equations.append(degen_eqn)
    '''
    return equations

if __name__ == '__main__':
    print(create_eqn_list(9, 11, 30))
