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
        if(ch == '\n') lines++;
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
    
    // Create bitmasks    
    FILE *fp = fopen(filename, "r");
    for(int i = 0; i < num_lines; i++) {
        char buff[255];
        fscanf(fp, "%s", buff);
        //cout << "buff: " << buff << endl;
        char *pt;
        pt = strtok(buff, ",");
        int counter = 0;
        uint64_t b_eqn = 0;
        while (pt != NULL) {
            int a = atoi(pt);
            if(counter < 3) {
                b_eqn += pow(2,a);
            }
            pt = strtok(NULL, ",");
            counter++;
        }
        // add to bitmask array
        equations[i] = b_eqn;
        b_eqn = 0;
    }
}

int main(int argc, char **argv) {
    // first arugment is equation file
    if(argc < 2) {
        cout << "not enough arguments, please specify equation file" << endl;
        return 0;
    }
    cout << "eqn file: " << argv[1] << endl;
    // read from file and create bitmask array
    read_file(argv[1]);
    return 0;
}
