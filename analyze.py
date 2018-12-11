import os, math, json, shutil



def get_value(max_bin, log_num_samples, sample_dict, gamma):
    D = abs(math.sin(gamma)) + abs(math.cos(gamma))

    ####  Get estimate.
    estimate = 0

    for i in range(max_bin+1):
        estimate += sample_dict[i][1]*D**i

    for i in range(log_num_samples):
        estimate /= 2

    #### Get Hoeffding error bound.

    delta = 0.01
    eConst = math.sqrt(math.log(2/delta) / 2)
    # hoeffding error = eConst * range / sqrt(numer of sam[ples)

    # this is unstable because big numbers are big
    # hoeffding = eConst * 2 * D**max_bin / 2**(log_num_samples/2)

    hoeffding = eConst * 2

    for i in range(max(max_bin,log_num_samples)):
        if i < max_bin: hoeffding *= D
        if i < log_num_samples: hoeffding /= math.sqrt(2)

    #### Get smart error bound

    delta = 1 - (1 - delta)*(1 - delta) # need two variables to be accurate w.p. 1-delta
    eConst = math.sqrt(math.log(2/delta) / 2)

    num_nonzero_samples = 0
    for i in range(max_bin+1): num_nonzero_samples += float(sample_dict[i][0])

    # probability of not sampling 0
    p = num_nonzero_samples
    for i in range(log_num_samples): p /= 2

    # error on that probability
    p_error = eConst
    for i in range(log_num_samples): p_error /= math.sqrt(2)

    # mean of nonzero samples
    nonzero_estimate = 0
    for i in range(max_bin+1):
        nonzero_estimate += sample_dict[i][1]*D**i
    nonzero_estimate /= num_nonzero_samples

    # error on nonzero samples mean
    nonzero_error = eConst * 2 / math.sqrt(num_nonzero_samples)
    for i in range(max_bin): nonzero_error *= D

    # combine error on p and error on nonzero samples
    error = p_error*(nonzero_estimate + nonzero_error)
    error += nonzero_error*(p + p_error)
    error += nonzero_error * p_error

    return estimate, error, hoeffding



if __name__ == "__main__":
    max_bin = 5 # 1 + (d-1)*3
    log_num_samples = 20
    sample_dict = {0:(12,12), 1:(42,12), 2:(0,0), 3:(4,10), 4:(12,10), 5:(12,0)}
    gamma = 0.12

    estimate, error, hoeffding = get_value(max_bin, log_num_samples, sample_dict, signed_sample_dict, gamma)

    print(estimate, error, hoeffding)


