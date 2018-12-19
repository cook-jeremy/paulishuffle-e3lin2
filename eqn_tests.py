import matplotlib as pyplot
import gen_equations
import random

eqns_fit = [[],[]]
vars_missed = [[],[]]
rep_eqn = [[],[]]
disconn = [[],[]]
solve = [[],[]]
good = [[],[]]

samples = 100

if True:
    # n vs f, fixed d
    d = 4

    plt.title("d = " + str(d))
    plt.xlabel("n (num vars)")
    plt.ylabel("f (num eqns)")

    for i in range(samples):
        n = random.choice(range(10,64))
        f = random.choice(range(20,50))

        out = gen_equations.gen_eqns(n, d, f, fail_kind = False)
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
    pass

plt.scatter(eqns_fit[0], eqns_fit[1], label="Eqns don't fit")
plt.scatter(vars_missed[0], vars_missed[1], label="Variables missed")
plt.scatter(rep_eqn[0], rep_eqn[1], label="Duplicater eqn")
plt.scatter(disconn[0], disconn[1], label="Disconnected")
plt.scatter(solve[0], solve[1], label="Solvable")
plt.scatter(good[0], good[1], label="Good")

plt.legend()
plt.show()
