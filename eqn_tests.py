import matplotlib.pyplot as plt
import gen_equations
import random
import numpy as np

eqns_fit = [[],[]]
vars_missed = [[],[]]
rep_eqn = [[],[]]
disconn = [[],[]]
solve = [[],[]]
good = [[],[]]

samples = 1000

if False:
    # n vs f, fixed d
    d = 20

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

if True:
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
