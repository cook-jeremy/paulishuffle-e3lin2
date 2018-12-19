import matplotlib.pyplot as plt
import gen_equations
import random
import numpy as np
import math

if True: # plot of family
    n = 20
    f = lambda d: 30 # lambda d: np.ceil(22.8*d - 34).astype(int)

    plt.title("Subcircuit Size\n f = 22.8*d - 34, n = 64")
    plt.xlabel("d constraint")

    d_range = [5] # range(4,15+1)
    ds = []
    max_eqns = []
    max_qubits = []


    # worst D
    log_num_samples = 30
    delta = 0.01
    eConst = math.sqrt(math.log(2/delta) / 2)
    D = abs(math.sin(math.pi/4)) + abs(math.cos(math.pi/4))

    for d in d_range:
        print("d", d, " f",f(d))
        for i in range(4):
            eqns, qubits = gen_equations.gen_eqns(n,d,f(d), cost_stats=True)
            ds.append(d)
            max_eqns.append(eqns)
            max_qubits.append(qubits)

        print("eqns:", eqns)
        hoeffding = eConst * 2 * (f(d)/2)

        for i in range(max(eqns,log_num_samples)):
            if i < eqns: hoeffding *= D
            if i < log_num_samples: hoeffding /= math.sqrt(2)

        print("err:", hoeffding)

    plt.scatter(ds, max_eqns, label="Equations")
    plt.scatter(ds, max_qubits, label="Qubits")
    plt.legend()
    plt.show()

if False: # scatter plot
    eqns_fit = [[],[]]
    vars_missed = [[],[]]
    rep_eqn = [[],[]]
    disconn = [[],[]]
    solve = [[],[]]
    good = [[],[]]

    samples = 1000

    if True:
        # n vs f, fixed d
        d = 5

        plt.title("d = " + str(d))
        plt.xlabel("n (num vars)")
        plt.ylabel("f (num eqns)")

        for i in range(samples):
            n = random.choice(range(5,60))
            f = random.choice(range(5,60))

            out = gen_equations.gen_eqns(n, d, f, fail_kind=True)
            if out == True:
                good[0].append(n)
                good[1].append(f)
            else:
                def append(l):
                    l[0].append(n)
                    l[1].append(f)

                if max(out) == out[0]: append(eqns_fit)
                elif max(out) == out[1]: append(vars_missed)
                elif max(out) == out[2]: append(rep_eqn)
                elif max(out) == out[3]: append(disconn)
                elif max(out) == out[4]: append(solve)

    if False:
        # d vs f, fixed n
        n = 64

        plt.title("n = " + str(n))
        plt.xlabel("d constraint")
        plt.ylabel("f (num eqns)")

        for i in range(samples):
            d = random.choice(range(3,25))
            f = random.choice(range(40,350))

            out = gen_equations.gen_eqns(n, d, f, fail_kind=True)
            if out == True:
                good[0].append(d)
                good[1].append(f)
            else:
                def append(l):
                    l[0].append(d)
                    l[1].append(f)

                if max(out) == out[0]: append(eqns_fit)
                elif max(out) == out[1]: append(vars_missed)
                elif max(out) == out[2]: append(rep_eqn)
                elif max(out) == out[3]: append(disconn)
                elif max(out) == out[4]: append(solve)

    s = 100

    if len(good[0]) > 0:
        coeffs = np.polyfit(good[0],good[1],1)
        print("y = ",np.round(coeffs[0],2), "+ x * ",np.round(coeffs[1],2))
    else:
        print("no successful params")

    plt.scatter(eqns_fit[0], eqns_fit[1], label="Eqns don't fit", s=s)
    plt.scatter(vars_missed[0], vars_missed[1], label="Variables missed", s=s)
    plt.scatter(rep_eqn[0], rep_eqn[1], label="Duplicater eqn", s=s)
    plt.scatter(disconn[0], disconn[1], label="Disconnected", s=s)
    plt.scatter(solve[0], solve[1], label="Solvable", s=s)
    plt.scatter(good[0], good[1], label="Good", s=s)

    plt.legend()
    plt.show()
