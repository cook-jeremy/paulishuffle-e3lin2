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

def create_eqns(filename, num_vars = 5, d_constraint = 4, num_eqns = 5):
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
    
    fp.close()
    fp_json = open(filename + ".json", 'w')

    to_write =  '{\"max_eqns\": ' + str(1+(d_constraint-1)*3) + ', \"num_eqns\": ' + str(num_eqns) + ', \"num_vars\": ' + str(num_vars) + ', \"d_constraint\": ' + str(d_constraint) + ', \"max_qubits\": ' + str(3+3*(d_constraint-1)*2) + '}'
    fp_json.write(to_write)
    fp_json.close()

    return gen_equations.get_cost_stats(input_eqns[0])

if __name__ == '__main__':
    # print eqns to file
    now = datetime.datetime.now()
    directory = 'equations/' + str(now.strftime('%m-%d-%Y')) + '/'
    if not os.path.exists(directory):
            os.makedirs(directory)
    filename = directory + str(now.strftime('%H:%M:%S%p'))
    create_eqns(filename)
