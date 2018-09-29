import random
import sys

def init_eqn_select(n, d, d_constraint):
    # init random equation selecter
    eqn_select = []
    for i in range(0,n):
        # if we haven't used d of this variable, add it to the list of potential eqns to be picked
        if d_constraint[i] < d:
            eqn_select.append(i)
    return eqn_select

# given n (number of variables), d (max eqs a single variable can appear in), g (number of equations)
def create_eqn_list(n, d, g):
    equations = []

    # only make f - 1 equations for degeneracy later
    f = 1
    if g > 1:
        f = g - 1

    d_constraint = [0]*n
    #print('d_constra = ' + str(d_constraint) + '\n')

    # loop over equations left to create 
    while f > 0:
        # create init eqn selecter [1,2,...,n]
        eqn_select = init_eqn_select(n, d, d_constraint)
        if(len(eqn_select) < 3): 
            print('ERROR too many constraints (d is too small)')
            sys.exit(1)

        #print('eqn start = ' + str(eqn_select))
        local_eqn = []

        for i in range(0,3):
            s = random.randint(0,len(eqn_select)-1)
            var = eqn_select[s]
            local_eqn.append(var)
            d_constraint[var] += 1
            eqn_select.remove(var)

        # pick solution
        local_eqn.append(random.randint(0,1))

        #print('local_eqn = ' + str(local_eqn))
        #print('d_constra = ' + str(d_constraint))
        #print('')
        equations.append(local_eqn)
        f = f - 1

    # make the equation set degenerate so we dont have exact solutions (is there a better way to do this?)
    if g > 1:
        opposite = 0
        if equations[0][3] == 0:
            equations[1][3] = 1
        degen_eqn = [equations[0][0],equations[0][1],equations[0][2],opposite]
        equations.append(degen_eqn)

    return equations

if __name__ == '__main__':
    print(create_eqn_list(10, 6, 20))
