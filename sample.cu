// reading a text file
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

struct matrix {
    vector< vector<int> > coeff;
    vector<int> sol;
};

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

vector<int>* read_file(char* filename) {
    string line;
    ifstream myfile(filename);
    vector< vector<int> > eqns;
    vector<int> sol;
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
    for(int i = 0; i < eqns.size(); i++) {
        cout << eqns.at(i).at(0) << " " << eqns.at(i).at(1) << " " << eqns.at(i).at(2) << " | " << sol.at(i) << endl;
    }
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
    return 0;
}
