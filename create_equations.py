import gen_equations
import datetime
import os
import numpy as np

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

def create_eqns(filename, num_vars = 20, d_constraint = 30, num_eqns = 45):
    # get input matrix and result, Ax = b
    input_eqns = gen_equations.gen_eqns(num_vars, d_constraint, num_eqns)
    # equation list
    el = convert_eqns(input_eqns)

    #with open(filename, 'wb') as fp:
    fp = open(filename, 'w')
    for i in range(0, len(el)):
        to_write = str(el[i][0]) + "," + str(el[i][1]) + "," + str(el[i][2]) + "," + str(el[i][3]) + "\n"
        fp.write(to_write)

        #pickle.dump(f_results, fp)

    return gen_equations.get_cost_stats(input_eqns[0])

if __name__ == '__main__':
    # print eqns to file
    now = datetime.datetime.now()
    directory = 'equations/' + str(now.strftime('%m-%d-%Y')) + '/'
    if not os.path.exists(directory):
            os.makedirs(directory)
    filename = directory + str(now.strftime('%H:%M:%S%p'))
    create_eqns(filename)
