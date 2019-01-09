import os, sys, json
import analyze

if __name__ == '__main__':
    f = open("output/run_params")
    run_params = json.loads(f.read())
    f.close()

    f = open(run_params["eqns_location"]+".json")
    eqn_params = json.loads(f.read())
    f.close()

    max_qubits = eqn_params["max_qubits"]
    max_eqns = eqn_params["max_eqns"]

    print("eqn:"+ str(run_params["eqns_location"]))
    print("d:"+ str(eqn_params["d_constraint"]))
    print("gamma:"+str(run_params["gamma"]))


    ################ Exact data

    exact_sum = 0
    exact_error = None

    for i in range(eqn_params["num_eqns"]):
        try:
            f = open("output/exact_%d.o" % (i,), "r")

            read = f.readline()
            print(read)
            if read == "": continue
            value = float(read)
            while value == int(value):
                value = float(f.readline())

            exact_sum += value/2

            file_error = float(f.readline())
            if exact_error is None: exact_error = file_error
            elif exact_error != file_error:
                raise ValueError("Files differ in van den nest error")

            f.close()
        except Exception as e:
            print("Error in output/exact_%d.o" % (i,))
            raise e

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

    # print("log samples: "+str(log_num_samples))
    # for i in range(len(list(samples.keys()))):
    #     print(str(i) + " - "+str(samples[i]))

    estimate, error, hoeffding = analyze.get_value(eqn_params["num_eqns"],\
            eqn_params["d_constraint"], log_num_samples, samples, run_params["gamma"], max_bin=max_eqns)

    ############## Export output

    print("GPU Estimate:\t"+str(estimate))
    print("Nest Estimate:\t"+str(exact_sum))

    print("Nest Error:\t"+str(abs(exact_error)))
    print("GPU Error:\t"+str(hoeffding))
    print("Smart GPU Error:\t"+str(error))

    print("GPU Error > \t"+str(abs(estimate-exact_sum)-exact_error))

