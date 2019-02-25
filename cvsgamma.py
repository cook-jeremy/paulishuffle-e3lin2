import pickle, datetime, os, sys, math, json
import exact, analyze, old_exact
import numpy as np
import matplotlib.pyplot as plt

def graph(exact, nest, pauli):
    #plt.plot(exact[0], exact[1], '-o', label= 'Exact', linewidth=2)

    # Plot with error bars
    plt.errorbar(nest[0], nest[1], yerr=nest[2], fmt='-o', label='Van den Nest', linewidth=2)
    plt.errorbar(pauli[0], pauli[1], yerr=pauli[2], fmt='-o', label='Pauli Shuffle', linewidth=2)

    plt.title('Estimate of <C> vs. Gamma')
    plt.legend()
    plt.ylabel('Estimate of <C>')
    plt.xlabel('Gamma')
    plt.show()

if __name__ == '__main__':
    if(len(sys.argv) < 2):
        print('missing <eqn_file>')
        sys.exit()

    # Van den Nest data
    nest_x = []
    nest_y = []
    nest_err = []
    nest = [nest_x, nest_y, nest_err]

    # Pauli Shuffle data
    pauli_x = []
    pauli_y = []
    pauli_err = []
    pauli = [pauli_x, pauli_y, pauli_err]

    num_steps = 25 
    PI = 3.141592653

    f = open(sys.argv[1] + ".json")
    eqn_params = json.loads(f.read())
    f.close()
    num_vars = eqn_params["num_vars"]
    d_constr = eqn_params["d_constraint"]
    num_eqns = eqn_params["num_eqns"]

    print('\nnest', end='')
    gamma = 0
    log_num_samples = 0
    # nest gpu
    for i in range(0, num_steps+1):
        nest_x.append(gamma)
        os.system('./exact %s %d %d %f > results/tmp' % (sys.argv[1], num_eqns, num_vars, gamma))
        tmp = open('results/tmp', 'r') 
        estimate = float(tmp.readline())
        error = float(tmp.readline())
        tmp.close()
        nest_y.append(estimate)
        nest_err.append(error)
        gamma += PI/num_steps
        sys.stdout.write('.')
        sys.stdout.flush()

    print('\npauli', end='')
    gamma = 0
    # gpu
    for i in range(0, num_steps+1):
        pauli_x.append(gamma)
        samples = {}
        os.system('./sample %s %d %d %f > results/tmp' % (sys.argv[1], num_eqns, d_constr, gamma))
        tmp = open('results/tmp', 'r') 
        log_num_samples = int(tmp.readline())
        for line in tmp:
            results = list(map(int, line.split(",")))
            key = results[0]
            if key not in samples:
                samples[key] = (results[1], results[2])
            else:
                samples[key] = (samples[key][0] + results[1],\
                        + samples[key][1] + results[2])
        tmp.close()
        estimate, error, hoeffding = analyze.get_value(num_eqns, d_constr, log_num_samples, samples, gamma)
        pauli_y.append(estimate)
        pauli_err.append(error)
        gamma += PI/num_steps
        sys.stdout.write('.')
        sys.stdout.flush()
    print('')

    # print results to file
    now = datetime.datetime.now()
    directory = 'results/' + str(now.strftime('%m-%d-%Y')) + '/'
    if not os.path.exists(directory):
        os.makedirs(directory)
    filename = directory + str(now.strftime('%H:%M:%S%p'))

    results = [nest, pauli]
    with open(filename, 'wb') as fp:
        pickle.dump(results, fp)

    graph(exact, nest, pauli)
