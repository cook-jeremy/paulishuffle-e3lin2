import os, math, json, shutil



def get_value(num_eqns, d_constraint, log_num_samples, sample_dict, gamma, max_bin=None):

    if max_bin is None:
        # get bound on number of quasiprobability samples
        max_bin = 3*(d_constraint-1) + 1
        #if (num_eqns < max_bin): max_bin = num_eqns
        num_eqns = float(num_eqns)

    D = abs(math.sin(gamma)) + abs(math.cos(gamma))

    ####  Get estimate.
    estimate = 0

    for i in range(max_bin+1):
        estimate += sample_dict[i][1]*D**i

    for i in range(log_num_samples):
        estimate /= 2

    estimate *= (num_eqns/2)

    #### Get Hoeffding error bound.

    delta = 0.01
    eConst = math.sqrt(math.log(2/delta) / 2)
    # hoeffding error = eConst * range / sqrt(numer of samples)

    # this is unstable because big numbers are big
    # hoeffding = eConst * 2 * D**max_bin / 2**(log_num_samples/2)

    hoeffding = eConst * 2 * (num_eqns/2)

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
        
    if num_nonzero_samples != 0:  
        nonzero_estimate /= num_nonzero_samples

    # error on nonzero samples mean
    nonzero_error = 0
    if num_nonzero_samples != 0:
        nonzero_error = eConst * 2 / math.sqrt(num_nonzero_samples)
    nonzero_error *= (num_eqns/2)
    for i in range(max_bin): nonzero_error *= D

    # combine error on p and error on nonzero samples
    error = p_error*(abs(nonzero_estimate) + nonzero_error)
    error += nonzero_error*(p + p_error)
    error += nonzero_error * p_error

    return estimate, error, hoeffding
