// reading a text file
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string>
#include <stdint.h>
#include <vector>
#include <sstream>
using namespace std;

struct matrix {
    vector< vector<int> > coeff;
    vector<int> sol;
}eqn_matrix;

// array of equations in bitmask form, i.e. x_2 + x_3 + x_4 for 5 variables is 01110
uint64_t *equations;

/**
int main(void) {
    // generate equation set
    

    // do ALL the gates in the circuit

    // pick random equation to create ...ZZZ...
    // apply e^{i\beta X} to above state (I -> I, Z -> Y)
    // apply e^{i\gamma C} to above state


    // pick random equation to create ...YYY...
    // apply e^{i\gamma C} to all equations that act on ...YYY...
    // take the trace with the |+> state to get observable
}
**/

void read_file(char* filename) {
    string line;
    ifstream myfile(filename);
    vector< vector<int> > eqns;
    vector<int> sol;
    uint64_t b_eqn = 0;
    if (myfile.is_open()) {
        while(getline(myfile,line)) {
            //cout << line << '\n';
            vector<int> eqn;
            stringstream ss(line);
            int i;
            int counter = 0;
            while (ss >> i) {
                if(counter < 3) eqn.push_back(i);
                else sol.push_back(i);
                if(ss.peek() == ',') ss.ignore();
                counter++;
            }
            eqns.push_back(eqn);
        }
        myfile.close();
    }
    else {
        cout << "Unable to open file" << endl; 
    }
    eqn_matrix.coeff = eqns;
    eqn_matrix.sol = sol;
}

int main(int argc, char **argv) {
    // first arugment is equation file
    if(argc < 2) {
        cout << "not enough arguments, please specify equation file" << endl;
        return 0;
    }
    cout << "eqn file: " << argv[1] << endl;
    // read from file
    read_file(argv[1]);

    for(int i = 0; i < eqn_matrix.coeff.size(); i++) {
        cout << eqn_matrix.coeff.at(i).at(0) << " " << eqn_matrix.coeff.at(i).at(1) << " " << eqn_matrix.coeff.at(i).at(2) << " | " << eqn_matrix.sol.at(i) << endl;
    }

    return 0;
}
