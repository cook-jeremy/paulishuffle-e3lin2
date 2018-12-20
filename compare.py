import pickle, datetime, os, sys, math
import exact, analyze

if __name__ == '__main__':
    if(len(sys.argv) < 2):
        print('missing <eqn_file>')
        sys.exit()
    x1 = []
    y1 = []
    x2 = []
    y2 = []
    f_results = [x1, y1, x2, y2]
    num_steps = 50
    gamma = 0
    PI = 3.141592653
    num_eqns = 6
    d_constr = 4
    num_vars = 5

    #num_eqns = 1
    #d_constr = 1
    #num_vars = 3

    #num_eqns = 10
    #d_constr = 7
    #num_vars = 6

    # exact
    for i in range(0, num_steps):
        x1.append(gamma)
        exact_sum = 0
        for i in range(num_eqns):
            exact_sum += float(exact.e3lin2_exact(i, sys.argv[1], gamma))/2
        y1.append(exact_sum)
        gamma += PI/num_steps

    # gpu
    gamma = 0
    log_num_samples = 0
    for i in range(0, num_steps):
        samples = {}
        print('%f' % gamma),
        os.system('./sample %s %f > results/tmp' % (sys.argv[1], gamma))
        tmp = open('results/tmp', 'r') 
        log_num_samples = int(tmp.readline())
        for line in tmp:
            results = map(int, line.split(","))
            key = results[0]
            if key not in samples:
                samples[key] = (results[1], results[2])
            else:
                samples[key] = (samples[key][0] + results[1],\
                        + samples[key][1] + results[2])
        tmp.close()
        #D = abs(math.cos(gamma)) + abs(math.sin(gamma))
        #x = (2**27)*abs(math.sin(gamma))/D
        #x2.append(x)
        x2.append(gamma)
        estimate, error, hoeffding = analyze.get_value(num_eqns, d_constr, log_num_samples, samples, gamma)
        #print(samples[1][0])
        #y2.append(samples[1][0])
        y2.append(estimate)
        print(estimate)
        gamma += PI/num_steps


    # print results to file
    now = datetime.datetime.now()
    directory = 'results/' + str(now.strftime('%m-%d-%Y')) + '/'
    if not os.path.exists(directory):
        os.makedirs(directory)
    filename = directory + str(now.strftime('%H:%M:%S%p'))

    with open(filename, 'wb') as fp:
        pickle.dump(f_results, fp)
