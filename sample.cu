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
//#include <experimental/filesystem>

#include "types.h" // circuit size specified here

#define PI 3.14159265358979323846

//Defined as powers of 2
#define samplesPerThread (long int)3  // Number of samples generated per thread.
#define threadsPerBlock (long int)10  // Number of threads per block.
#define blocksPerChunk (long int)10   // Number of blocks per output array.
#define numChunks (long int) 2        // Do the whole thing each time for a new gamma

//Also powers of 2
#define samplesPerChunk samplesPerThread + threadsPerBlock + blocksPerChunk
#define nsamples numChunks + samplesPerChunk

// define tally type
#define tally_t int // Maybe switch to long int?

using namespace std;

__global__ void sample(int seed);

// Device Memory
__device__ __constant__ uint64_t d_eqn_masks[tally_array_size];
__device__ tally_type d_chunk_tally[tally_array_size]; // Global memory
// Chunktally stores the number of samples with n encountered D's
// in d_chunk_tally[2*t], and the tally taking into acount sign in d_chunk_tally[2*t+1]
   
// array of equations in bitmask form, i.e. x_2 + x_3 + x_4 for 5 variables is 01110
uint64_t *h_eqn_masks;
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
    
    cout << "finished reading file" << endl;

    //if (!print_device_properties()) return 0;

    // First malloc space on device for d_eqn_masks, since it depends on num_eqns
    cudaMalloc(d_eqn_masks, num_eqns*sizeof(uint64_t));
    // Then copy equation bit mask array to device
    cudaMemcpyToSymbol(d_eqn_masks, h_eqn_masks, num_eqns*sizeof(uint64_t));
    // We don't need the masks on the host.
    free(h_eqn_masks);  

    // Host memory for tallying output.
    int tally_array_size = 2*num_eqns;
    tally_t* h_chunk_tally = (tally_t*) malloc(tally_array_size*sizeof(tally_t));
    tally_t* output_tally = (tally_t*)malloc(tally_array_size*sizeof(tally_t));
    // Initialize both arrays to 0
    for (int i = 0; i < tally_array_size; i++) {
        h_chunk_tally[i] = 0;
        output_tally[i] = 0;
    }

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
        cudaMemcpyFromSymbol(h_chunk_tally, d_chunk_tally, tally_array_size*sizeof(tally_t));
        cudaMemset(d_chunk_tally, 0, tally_array_size*sizeof(tally_t));

        // Add chunk tally to overall tally
        for (int i = 0; i < tally_array_size; i++) output_tally[i] += h_chunk_tally[i];
    }

    // print output
    std::cout << nsamples << std::endl;
    for (int i = 0; 2*i < tally_array_size; i+=1) {
        std::cout << i << "," << output_tally[2*i] << "," << output_tally[2*i+1] << std::endl;
    }

    // Free memory
    free(h_chunk_tally);
    free(output_tally);

    return 0;
}


__device__ void parity(int* x, uint64_t data) {
    // Flip the sign of x if data has an odd number of 1's in it
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
    uint64_t xs, zs;
    tally_t num_D; // number of D(\gamma)
    int sign; // store separately

    //int index = blockIdx.x * blockDim.x + threadIdx.x;
    //int stride = blockDim.x * gridDim.x;
    
    for(int j = 0; j < (1 << samplesPerThread); j++) {
        // Pick a random equation from eqn_masks
        int rand = 0; // TODO how to get a rand int on device from 0 to num_eqns - 1
        uint64_t init_mask = eqn_masks[rand];
        xs = mask;
        zs = mask;
         
        for(int i = 0; i < num_eqns; i++) {
            uint64_t mask = eqn_masks[i];
            int *test = 1;
            parity(test, mask & xs);
            if(*test == -1) {
                // Doesn't commute
                rand = 0; // TODO generate rand float between 0 and 1

            }

            /** TODO FINISH THIS FOR LOOP **/
        }

        if (zs == 0 && (xs & zs) == 0) { // <0|X|0> = <1|X|1> = <0|Y|0> = <1|Y|1> = 0
            // Write to global output memory. Use atomic add to avoid sync issues.
            atomicAdd(&d_chunkTally[num_D*2], (tally_t)1);
            atomicAdd(&d_chunkTally[num_D*2+1], (tally_t)sign);
        }
    }
}
