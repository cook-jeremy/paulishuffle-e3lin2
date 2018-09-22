import numpy as np
import math

I = np.matrix('1 0; 0 1')
X = np.matrix('0 1; 1 0')
Y = np.matrix('0 -1;1 0')*complex(0,1)
Z = np.matrix('1 0; 0 -1')
paulis = [I,X,Y,Z]

# decompose \rho = e^{i\beta\sigma_x} \sigma_j e^{-i\beta\sigma_x}
def decompose(beta, pauli):
    B = np.matrix([[math.cos(beta),complex(0,math.sin(beta))],[complex(0,math.sin(beta)),math.cos(beta)]])
    nB = np.matrix([[math.cos(beta),complex(0,-math.sin(beta))],[complex(0,-math.sin(beta)),math.cos(beta)]])
    R = B*pauli*nB
    c = []
    for P in paulis:
        c_p = (0.5*P*R).trace()
        c.append(np.asscalar(c_p))
    return c

#if __name__ == '__main__':
    #decompose(math.pi/2,Y)
