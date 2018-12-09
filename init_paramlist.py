import create_equations
import numpy as np

username = "dliang27"

num_vars = 20
d_constraint = 30
num_eqns = 45
noise_p = 0
nodes = 1

eqns_location = "equations/figure_eqns/equations"

def init_sample():
    create_equations.create_eqns(eqns_location, num_vars, d_constraint, num_eqns, noise_p)
    gammalist = [0.2]
    for gamma in gammalist:
        f = open("paramlist", "w")
        for i in range(nodes):
	    f.write("sample.exe equations/figure_eqns/equations %g > output/tally%g_%d.o 2> output/tally%g_%d.e\n" % (gamma, gamma, i+1, gamma, i+1))
            f.write("python exact.py %d %d %s %g > output/exact%g.o 2> output/exact%g.e\n" % (num_vars, num_eqns, eqns_location, gamma, gamma, gamma))
	f.close()

def sample():
    pass

def main():
    init_sample()


if __name__ == '__main__':
    main()
