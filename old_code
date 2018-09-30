H = np.matrix('1 1; 1 -1')*(1/math.sqrt(2))

def classical_circuit(INIT):
    result = (Z*X*H*INIT*H*X).trace()
    print(result)

def circuit(INIT):
    results = []
    samples = 10000
    for i in range(0,samples):
        # pick a pauli for INIT
        r1 = pick_pauli(INIT)
        STATE = r1[0]

        STATE = H*STATE*H
        r2 = pick_pauli(STATE)
        STATE = r2[0]

        STATE = X*STATE*X
        r3 = pick_pauli(STATE)
        STATE = r3[0]
        
        rf = r1[1]*r2[1]*r3[1]
        pf = r1[2]*r2[2]*r3[2]
        
        p_hat = (rf/pf)*(STATE*Z).trace()
        results.append(p_hat)

    avg = sum(results)/samples
    print(np.asscalar(avg))

 '''
    local_state = kron(states[0])
    eiC = np.asmatrix(expm(complex(0,1)*gamma*kron([Z,Z,Z])))
    decomp = pick_pauli(eiC*local_state*np.conj(eiC))
    states[1] *= decomp[1]
    states[2] *= decomp[2]
    return [decomp[0], states[1], states[2]]
    '''

