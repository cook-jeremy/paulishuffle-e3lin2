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
int count_lines(char *filename) {
    // count the number of lines in the file called filename                                    
    FILE *fp = fopen(filename,"r");
    int ch=0;
    int lines=0;

    if (fp == NULL) return 0;

    while(!feof(fp)) {
        ch = fgetc(fp);
        if(ch == '\n') {
            lines++;
        }
    }
    fclose(fp);
    return lines;
}

void read_file(char* filename) {
    // get total number of equations and malloc bitmask array
    int num_lines = count_lines(filename);
    equations = (uint64_t *) malloc(num_lines*sizeof(uint64_t));

    //FILE *eqn_file = fopen(filename, "r");
    cout << "number of lines: " << num_lines << endl;

    
    string line;
    ifstream myfile(filename);

    //vector< vector<int> > eqns;
    //vector<int> sol;
    uint64_t b_eqn = 0;

    if (myfile.is_open()) {
        while(getline(myfile,line)) {
            //cout << line << '\n';
            vector<int> eqn;
            stringstream ss(line);
            int i;
            int counter = -1;
            while (ss >> i) {
                if(counter < 3) {
                    eqn.push_back(i);
                    b_eqn += pow(2,i);
                } else { 
                    sol.push_back(i);
                }
                if(ss.peek() == ',') ss.ignore();
                counter++;
            }
            //eqns.push_back(eqn);
            equations[counter] = b_eqn;
        }
        myfile.close();
    }
    else {
        cout << "Unable to open file" << endl; 
    }
    //eqn_matrix.coeff = eqns;
    //eqn_matrix.sol = sol;
    
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
