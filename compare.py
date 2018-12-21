import pickle, datetime, os, sys, math
import exact, analyze, nest
import numpy as np
import matplotlib.pyplot as plt

def graph(data):
    #if dbg: print(data)
    for i in range(len(data[1])):
        print(data[0][i], "\t", data[1][i], "\t", data[2][i], "\t", data[3][i])
    plt.plot(data[0], data[1], '-o', label= 'exact', linewidth=2)
    plt.plot(data[2], data[3], '-o', label= 'nest', linewidth=2)
    #plt.errorbar(data[2], data[3], yerr=data[4], fmt='-o', label='nest', linewidth=2)
    #plt.plot(data[0], data[3], '-o', label= 'pauli', linewidth=2)
    plt.errorbar(data[5], data[6], yerr=data[7], fmt='-o', label='pauli', linewidth=2)
    plt.title('Exact, Nest, Pauli')
    plt.legend()
    plt.ylabel('<C>')
    plt.xlabel('gamma')
    plt.show()

if __name__ == '__main__':
    if(len(sys.argv) < 2):
        print('missing <eqn_file>')
        sys.exit()
        
    params_location = sys.argv[1] + '.params'
    x1 = []
    x2 = []
    x3 = []
    y1 = []
    y2 = []
    y2err = []
    y3 = []
    y3err = []
    f_results = [x1, y1, x2, y2, y2err, x3, y3, y3err]
    num_steps = 25
    PI = 3.141592653

    params = open(params_location, "r").readline().rstrip().split(',')
    for i in range(len(params)): params[i] = int(params[i])
    num_vars = params[0]
    d_constr = params[1]
    num_eqns = params[2]

    gamma = 0
    # exact
    print('exact', end='')
    for i in range(0, num_steps+1):
        x1.append(gamma)
        exact_sum = 0
        for i in range(num_eqns):
            exact_sum += float(exact.e3lin2_exact(i, sys.argv[1], gamma))/2
        y1.append(exact_sum)
        gamma += PI/num_steps
        sys.stdout.write('.')
        sys.stdout.flush()

    gamma = 0.01
    # nest
    print('\nnest', end='')
    for i in range(0, num_steps+1):
        x2.append(gamma)
        exact_sum = 0
        for i in range(num_eqns):
            estimate, error = nest.e3lin2_exact(i, sys.argv[1], gamma)
            exact_sum += estimate/2
        #print('%f %f' % (gamma, exact_sum))
        y2.append(exact_sum)
        y2err.append(error)
        gamma += PI/num_steps
        sys.stdout.write('.')
        sys.stdout.flush()

    print('\ngpu', end='')
    gamma = 0.02
    log_num_samples = 0
    # gpu
    for i in range(0, num_steps+1):
        x3.append(gamma)
        samples = {}
        os.system('./sample %s %f > results/tmp' % (sys.argv[1], gamma))
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
        y3.append(estimate)
        y3err.append(error)
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

    with open(filename, 'wb') as fp:
        pickle.dump(f_results, fp)

    graph(f_results)
