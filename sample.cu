// reading a text file
#include <stdio.h>
#include <fstream>
//#include <string>

#include <iostream>
#include <limits>
#include <iomanip>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <curand.h>
#include <curand_kernel.h>
#include <experimental/filesystem>

using namespace std;

// array of equations in bitmask form, i.e. x_2 + x_3 + x_4 for 5 variables is 01110
uint64_t *eqn_masks;
int num_eqns;

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
    num_eqns = count_lines(filename);
    eqn_masks = (uint64_t *) malloc(num_eqns*sizeof(uint64_t));

    //FILE *eqn_file = fopen(filename, "r");
    cout << "number of lines: " << num_eqns << endl;
    
    // Create bitmasks    
    FILE *fp = fopen(filename, "r");
    for(int i = 0; i < num_eqns; i++) {
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
        eqn_masks[i] = b_eqn;
        b_eqn = 0;
    }
    fclose(fp);
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

    // Read circuit. This should be in constant memory.
    //Gate* h_circ = (Gate*)malloc(circuitSize); 
    //FILE *circFile = fopen("circuit.circ", "rb");
    //fread(h_circ, sizeof(char), circuitSize, circFile);
    //fclose(circFile);
    cudaMemcpyToSymbol(c_circ, h_circ, circuitSize);
    free(h_circ);  // we don't need the circuit on the host.

    // Host memory for tallying output.
    tallyType* h_chunkTally = (tallyType*)malloc(tallyArraySize*sizeof(tallyType));
    tallyType* outputTally = (tallyType*)malloc(tallyArraySize*sizeof(tallyType));
    for (int i = 0; i < tallyArraySize; i++) outputTally[i] = 0;

    // Random seed
    // time_t t;
    // srand((unsigned) time(&t)); // random version
    srand(0);

    for (int j = 0; j < (1 << numChunks); j++) {
        //std::cout << "Running chunk " << (j+1) << " of " << (1 << numChunks) << std::endl;

        // Take samples
        sample<<<(1 << blocksPerChunk), (1 << threadsPerBlock)>>>(time(0)); //random version
        //sample<<<blocksPerChunk, threadsPerBlock>>>(j);//deterministic version

        // Wait for GPU to finish before accessing on host
        cudaDeviceSynchronize();

        // Copy samples to host, zero out device data
        cudaMemcpyFromSymbol(h_chunkTally, d_chunkTally, tallyArraySize*sizeof(tallyType));
        cudaMemset(d_chunkTally, 0, tallyArraySize*sizeof(tallyType));

        // Add chunk tally to overall tally
        for (int i = 0; i < tallyArraySize; i++) outputTally[i] += h_chunkTally[i];
    }

    // print output
    std::cout << nsamples << std::endl;
    for (int i = 0; 2*i < tallyArraySize; i+=1) {
        std::cout << i << "," << outputTally[2*i] << "," << outputTally[2*i+1] << std::endl;
    }

    // Free memory
    free(h_chunkTally);
    free(outputTally);

    return 0;
}


__device__ void parity(int* x, uint64_t data) {
    /* Alternate method? https://graphics.stanford.edu/~seander/bithacks.html#ParityWith64Bits
    unsigned char b;  // byte value to compute the parity of
    bool parity = (((b * 0x0101010101010101ULL) & 0x8040201008040201ULL) % 0x1FF) & 1;
    // multiplication may be slower than just iterating through?
    */

    while (data) {
        *x *= -1;
        data = data & (data - 1);
    }

}

__device__ void random(curandState_t* state, uint64_t* output) {
    *output =  (((uint64_t)curand(state)) << 32) | (uint64_t)curand(state);
}

__global__ void sample(int seed) {

    // Initialize curand
    curandState_t state;
    curand_init(seed, blockIdx.x, threadIdx.x, &state);

    // Per thread local memory. Can probably make this smaller with uglier code.
    uint64_t xs, zs, rand_data, toFlip, scales, upshiftedxs, upshiftedzs;
    tallyType out; // keep temporary answer in local memory
    int sign; // store separately

    //int index = blockIdx.x * blockDim.x + threadIdx.x;
    //int stride = blockDim.x * gridDim.x;
    
    for(int j = 0; j < (1 << samplesPerThread); j++) {
        
        // Pick a random equation from eqn_masks
        uint64_t init_mask = eqn_masks[randint between 0 and num_eqns];
        xs = mask;
        zs = mask;
        
        for(int i = 0; i < num_eqns; i++) {
            uint64_t mask = eqn_masks[i];
            int *test = 1;
            parity(test, mask & xs);
            if(*test == -1) {
                cout << "doesn't commute" << endl;
                // where is gamma in all of this?
            }
        }




        /**
        // either X or I at random
        random(&state, &xs);
        xs &= (0xffffffffffffffff >> (64 - width*height));
        zs = 0;

        out = 0;
        sign = 1;

        for (int i = 0; i < depth*gatesPerLayer; i++) {
            switch(c_circ[i].type) {
                case ROOTX:
                    parity(&sign, ~xs & zs & c_circ[i].data);
                    xs ^= c_circ[i].data & zs;
                    break;
                case ROOTY:
                    parity(&sign, ~xs & zs & c_circ[i].data);
                    toFlip = (zs ^ xs) & c_circ[i].data;
                    xs ^= toFlip;
                    zs ^= toFlip;
                    break;     
                case TGATE: // Only gate that is random
                    random(&state, &rand_data);
                    rand_data &= c_circ[i].data;

                    parity(&sign, rand_data & xs & zs);
                    
                    scales = xs & c_circ[i].data;
                    while (scales) {
                        out += 1;
                        //out *= sqrt2;
                        scales = scales & (scales - 1);
                    }
                    zs ^= xs & rand_data;

                    break;
                case CPHASERIGHT:
                case CPHASEDOWN:
                    if (c_circ[i].type == CPHASERIGHT) {
                        upshiftedxs = ((c_circ[i].data << 1) & xs) >> 1;
                        upshiftedzs = ((c_circ[i].data << 1) & zs) >> 1;
                        zs ^= (c_circ[i].data & xs) << 1;
                    } else {
                        upshiftedxs = ((c_circ[i].data << width) & xs) >> width;
                        upshiftedzs = ((c_circ[i].data << width) & zs) >> width;
                        zs ^= (c_circ[i].data & xs) << width;
                    }
                    
                    zs ^= upshiftedxs;
                    parity(&sign, c_circ[i].data & xs & upshiftedxs & ((c_circ[i].data & zs) ^ upshiftedzs));
                    break; 
                default:
                    break; // silence warning
            }
        }
        **/
        // Compare to output string to get probability
		if (xs == 0) { // <0|X|0> = <1|X|1> = <0|Y|0> = <1|Y|1> = 0
            parity(&sign, zs & c_circ[depth*gatesPerLayer].data); // Only <1|Z|1> = -1

            // Write to global output memory. Use atomic add to avoid sync issues.
            atomicAdd(&d_chunkTally[out*2], (tallyType)1);
            atomicAdd(&d_chunkTally[out*2+1], (tallyType)sign);
        }
    }
}

