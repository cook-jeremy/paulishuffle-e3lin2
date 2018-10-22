#include <stdio.h>

// return array of coefficient matrix and solution vector
int** gen_eqns(int n, int d, int f) {
    bool satisfied = false;
    int num_failures = 0;
    /**
    while(!satisfied) {
        // create A and b in Ax = b
        //sys = create_eqns(n, d, f)
        
        // check if we have used all variables
        if(!contains_all_vars(sys[0])) {
            // print('doesn\'t contain all variables, trying again...')
            num_failures += 1;
            continue
        }

        // check for the same equations
        if(!unique_eqns(sys[0], sys[1])) {
            num_failures += 1;
            continue
        }

        // check if we are fully connected
        if(!fully_connected(sys[0])) {
            // print('ins\'t fully connected, trying again...')
            num_failures += 1;
            continue
        }

        // row reduce and check if system is solvable
        if(solvable(sys[0], sys[1])) {
            num_failures += 1;
            continue
        }

        // passed all conditions
        break
    }

    //print('num failures: %s' % num_failures)
    return sys
    **/
    int sys[2][3] = {{1,1,1},{0,0,0}};
    return sys;
}

int main(void) {
    int ret[][] = gen_eqns(9, 100, 9)
    //print(ret[0])
    //print(ret[1])
    return 0;
}
