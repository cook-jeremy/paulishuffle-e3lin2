import os, sys, json, datetime, time
import exact, analyze, old_exact

if __name__ == '__main__':
    f = open("output/run_params")
    run_params = json.loads(f.read())
    f.close()

    f = open(run_params["eqns_location"] + ".json")
    eqn_params = json.loads(f.read())
    f.close()

    for i in range(eqn_params["num_eqns"]):
        st = datetime.datetime.fromtimestamp(time.time()).strftime('%H:%M:%S')
        print("eqn %d of %d, time: %s" % (i+1, eqn_params["num_eqns"], st))
        os.system('./exact %d %s %f > output/exact_%d.o' % (i, run_params["eqns_location"], run_params["gamma"], i))

    print("pauli")
    os.system('./sample %s %f > output/tally_0.o' % (run_params["eqns_location"], run_params["gamma"]))
