import os, sys, json
import analyze

if __name__ == '__main__':
    f = open("output/run_params")
    run_params = json.loads(f.read())
    f.close()

    f = open(run_params["eqns_location"]+".json")
    eqn_params = json.loads(f.read())
    f.close()

    ################# GPU Data

    log_num_samples = None
    samples = {}

    for i in range(2**run_params["log_num_gpus"]):
        f = open("output/tally_%d.o" % (i,), "r")

        file_samples = int(f.readline())
        if log_num_samples is None: log_num_samples = file_samples
        elif log_num_samples != file_samples:
            raise ValueError("Files differ in number of samples taken.")

        while True:
            line = f.readline()
            if not line: break

            results = map(int, line.split(","))
            key = results[0]
            if key not in samples:
                samples[key] = (results[1], results[2])
            else:
                samples[key] = (samples[key][0] + results[1],\
                        + samples[key][1] + results[2])
        f.close()

    log_num_samples += run_params["log_num_gpus"]

    print("log samples: "+str(log_num_samples))
    for i in range(len(list(samples.keys()))):
        print(str(i) + " - "+str(samples[i]))

    estimate, error, hoeffding = analyze.get_value(eqn_params["num_eqns"],\
            eqn_params["d_constraint"], log_num_samples, samples, run_params["gamma"])

    ################ Exact data

    exact_sum = 0

    for i in range(eqn_params["num_eqns"]):
        f = open("output/exact_%d.o" % (i,), "r")
        exact_sum += float(f.read())/2
        f.close()

    ############## Export output

    print("Estimate:\t"+str(estimate))
    print("Exact answer:\t"+str(exact_sum))
    print("Actual error:\t"+str(abs(estimate-exact_sum)))
    print("Error bound:\t"+str(error))
    print("Hoeffding:\t"+str(hoeffding))

