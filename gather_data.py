import os, sys


def gather_data_helper(samples, f):
    log_num_sample = int(f.readline())
    while True:
        line = f.readline()
        if not line:
            break
        results = map(int, line.split(","))
        key = results[0]
        if key not in samples:
            samples[key] = (results[1], results[2])
        else:
            samples[key] = (samples[key][0] + results[1], + samples[key][1] + results[2])

def gather_data(gamma):
    samples = dict() 
    i = 0
    file_location = "output/tally%g_%d.o" % (gamma, i+1)
    while os.path.exists(file_location):
        f = open(file_location, "r")
        gather_data_helper(samples, f)
        f.close()
        i += 1
        file_location = "output/tally%g_%d.o" % (gamma, i+1)
    return samples
      

def main():
    print(gather_data(0.2))

if __name__ == '__main__':
    main()
