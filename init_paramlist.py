import create_equations
import sys, os
import json, datetime

gamma = 0.12
noise_p = 0.1 # unused at the moment
log_num_gpus = 0

# eqns_location = "equations/test"
eqns_location = "none" # creates new equations

if eqns_location == "none":
    num_vars = 21
    d_constraint = 4
    num_eqns = 28

def init_sample():
    f = open("paramlist_gpu", "w")

    for i in range(2**log_num_gpus):
        f.write("sample.exe %s %g > output/tally_%d.o 2> output/tally_%d.e\n"\
                % (eqns_location, gamma, i, i))

    f.close()

def init_exact():
    f = open("paramlist_exact", "w")

    for i in range(num_eqns):
        f.write("python exact.py %d %s %g > output/exact_%d.o 2> output/exact_%d.e\n"\
                % (i, eqns_location, gamma, i, i))

    f.close()

if __name__ == '__main__':
    # clear output dir
    for fname in os.listdir("output/"):
        os.unlink(os.path.join("output/", fname))

    if eqns_location == "none":
        now = datetime.datetime.now()
        eqns_location = str(now.strftime('equations/%H-%M-%S-%m-%d-%Y'))
        print("generating new equations at " + eqns_location)

        create_equations.create_eqns(eqns_location, num_vars, d_constraint, num_eqns)

        f = open(eqns_location+".json", "w")
        f.write(json.dumps({
            "num_vars": num_vars,
            "d_constraint": d_constraint,
            "num_eqns": num_eqns,
        }))
        f.close()
    else:
        f = open(eqns_location+".json")
        data = json.loads(f.read())
        f.close()
        num_eqns = data["num_eqns"]

    f = open("output/run_params", "w")
    f.write(json.dumps({
        "gamma":gamma,
        "noise_p": noise_p,
        "log_num_gpus": log_num_gpus,
        "eqns_location": eqns_location,
    }))
    f.close()

    init_sample()
    init_exact()

