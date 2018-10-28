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
#include <inttypes.h>
//#include <experimental/filesystem>

#include "types.h" // circuit size specified here

#define PI 3.14159265358979323846

//Defined as powers of 2
#define samplesPerThread (long int)1  // Number of samples generated per thread.
#define threadsPerBlock (long int)1   // Number of threads per block.
#define blocksPerChunk (long int)1    // Number of blocks per output array.
#define numChunks (long int) 1        // Do the whole thing each time for a new gamma

//Also powers of 2
#define samplesPerChunk samplesPerThread + threadsPerBlock + blocksPerChunk
#define nsamples numChunks + samplesPerChunk

// define tally type
#define tally_t int // Maybe switch to long int?

using namespace std;

__global__ void sample(int seed, int num_eqns, uint64_t *d);

// Device Memory
#define SN 3
__device__ __constant__ double s_array[SN];
__device__ tally_t *d_chunk_tally; // Global memory
// Chunktally stores the number of samples with n encountered D's
// in d_chunk_tally[2*t], and the tally taking into acount sign in d_chunk_tally[2*t+1]
   
// array of equations in bitmask form, i.e. x_2 + x_3 + x_4 for 5 variables is 01110
uint64_t *h_eqn_masks;
int num_eqns;
//double gamma;

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
    h_eqn_masks = (uint64_t *) malloc(num_eqns*sizeof(uint64_t));

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
        h_eqn_masks[i] = b_eqn;
        b_eqn = 0;
    }
    fclose(fp);
}

int main(int argc, char **argv) {
    // first arugment is equation file
    if(argc < 3) {
        cout << "not enough arguments, please specify <equation file> and <gamma>" << endl;
        return 0;
    }
    cout << "eqn file: " << argv[1] << endl;
    cout << "gamma: " << argv[2] << endl;
    // TODO convert char array to float
    //gamma = argv[2];
    double gamma = strtod(argv[2],NULL);

    // read from file and create bitmask array
    read_file(argv[1]);
    
    cout << "finished reading file" << endl;

    //if (!print_device_properties()) return 0;

    // First malloc space on device for d_eqn_masks, since it depends on num_eqns
    printf("mallocing memory for eqn bit masks\n");

    uint64_t *d;
    cudaMalloc((void **)&d, num_eqns*sizeof(uint64_t));
    // Then copy equation bit mask array to device
    cudaMemcpy(d, h_eqn_masks, num_eqns*sizeof(uint64_t), cudaMemcpyHostToDevice);


    double c_array[SN];
    double c_tot = abs(sin(gamma)) / (abs(cos(gamma)) + abs(sin(gamma)));
    double sign_s = 1;
    double sign_c = 1;
    if(sin(gamma) < 0) sign_s = -1;
    if(cos(gamma) < 0) sign_c = -1;
    c_array[0] = c_tot;
    c_array[1] = sign_s;
    c_array[2] = sign_c;
    //printf("|s|/(|s|+|c|) = %f\n", s);
    cudaMemcpyToSymbol(s_array, c_array, SN*sizeof(double));
    //double *s_ptr;
    //cudaMalloc((void **) &s_ptr, SN*sizeof(double));
    //cudaMemcpy(d_s, &s, sizeof(double), cudaMemcpyHostToDevice);

    /**
    for(int i = 0; i < num_eqns; i++) {
        //printf("h_mask %d:%ud\n", i, h_eqn_masks[i]);
        printf("eqn %d: %" PRIu64 "\n", i, h_eqn_masks[i]);
    }
    **/
    // We don't need the masks on the host.
    free(h_eqn_masks);

    // Host memory for tallying output.
    int tally_size = 2*num_eqns;
    tally_t* h_chunk_tally = (tally_t*) malloc(tally_size*sizeof(tally_t));
    tally_t* output_tally = (tally_t*)malloc(tally_size*sizeof(tally_t));

    // Initialize both arrays to 0
    for (int i = 0; i < tally_size; i++) {
        h_chunk_tally[i] = 0;
        output_tally[i] = 0;
    }

    // Random seed
    // time_t t;
    // srand((unsigned) time(&t)); // random version
    srand(0);

    //for (int j = 0; j < (1 << numChunks); j++) {
    for(int j = 0; j < numChunks; j++) {
        //std::cout << "Running chunk " << (j+1) << " of " << (1 << numChunks) << std::endl;
        std::cout << "Running chunk " << (j+1) << " of " << numChunks << std::endl;

        // Take samples
        //sample<<<(1 << blocksPerChunk), (1 << threadsPerBlock)>>>(time(0)); //random version
        sample<<<1, threadsPerBlock>>>(time(0), num_eqns, d); //random version
        //sample<<<blocksPerChunk, threadsPerBlock>>>(j);//deterministic version

        // Wait for GPU to finish before accessing on host
        cudaDeviceSynchronize();

        // Copy samples to host, zero out device data
        //cudaMemcpyFromSymbol(h_chunk_tally, d_chunk_tally, tally_size*sizeof(tally_t));
        //cudaMemset(d_chunk_tally, 0, tally_size*sizeof(tally_t));

        // Add chunk tally to overall tally
        //for (int i = 0; i < tally_size; i++) output_tally[i] += h_chunk_tally[i];
    }

    // print output
    /**
    std::cout << nsamples << std::endl;
    for (int i = 0; 2*i < tally_size; i+=1) {
        std::cout << i << "," << output_tally[2*i] << "," << output_tally[2*i+1] << std::endl;
    }

    // Free memory
    free(h_chunk_tally);
    free(output_tally);
    */

    return 0;
}

__device__ void printb(size_t const size, void const * const ptr)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;

    for (i=size-1;i>=0;i--) {
        for (j=7;j>=0;j--) {
            byte = (b[i] >> j) & 1;
            printf("%u", byte);
        }
    }
    printf("");
    //puts("");
}

__device__ void print_xs_zs(uint64_t xs, uint64_t zs) {
    printf("xs: ");
    printb(sizeof(xs), &xs);
    printf("\n");
    printf("zs: ");
    printb(sizeof(zs), &zs);
    printf("\n");
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

// get a uniformly random integer inclusively between min and max
__device__ int get_rand_int(curandState_t state, int min, int max) {
    float rand_f = curand_uniform(&state);
    rand_f *= (max - min + 0.999999);
    rand_f += min;
    return (int)truncf(rand_f);
}

__global__ void sample(int seed, int num_eqns, uint64_t *d_eqn_masks) {
    // Initialize curand
    curandState_t state;
    curand_init(seed, blockIdx.x, threadIdx.x, &state);

    printf("block, thread %d, %d\n", blockIdx.x, threadIdx.x);
    
    printf("STARTING IN DEVICE CODE NOW\n");
    for(int i = 0; i < 3; i++) {
        printf("s_array[%d]: %f\n", i, s_array[i]);
    }
    /**
    for(int i = 0; i < num_eqns; i++) {
        printf("eqn %d: %" PRIu64 "\n", i, d[i]);
    }
    **/

    //printf("random: %d\n", get_rand_int(state, 0, 10));
    printf("num eqns: %d\n", num_eqns);
    
    // Per thread local memory. Can probably make this smaller with uglier code.
    uint64_t xs, zs;
    tally_t num_D; // number of D(\gamma)
    int sign = 1; // store separately
    
    //for(int j = 0; j < (1 << samplesPerThread); j++) {
    for(int j = 0; j < samplesPerThread; j++) {
        // Pick a random equation from eqn_masks
        int rand = get_rand_int(state, 0, num_eqns - 1);
        printf("rand: %d\n", rand);
        printf("--------- INIT ---------\n");
        uint64_t init_mask = d_eqn_masks[rand];
        xs = init_mask;
        zs = init_mask;

        print_xs_zs(xs, zs);
        printf("-------- Applying e^{i gamma C} --------\n"); 
        for(int i = 0; i < num_eqns; i++) {
            uint64_t mask = d_eqn_masks[i];
            printf("pq: ");
            printb(sizeof(uint64_t), &mask);
            printf("\n");
            print_xs_zs(xs, zs);
            int test = 1;
            parity(&test, mask & xs);
            if(test == -1) {
                // Doesn't commute
                float rand_f = curand_uniform(&state);
                printf("rand float: %f\n", rand_f);
                if(rand_f <= s_array[0]) {
                    printf("in S branch\n");
                    // apply iZZZ\rho TODO
                    if(s_array[1] < 0) sign = -1*sign;
                } else if(s_array[2] < 0) sign = -1*sign;
                num_D += 1;            
            }
        }
        
        // Because <+|Y|+> = <+|Z|+> = 0, we only care if both of these don't happen
        if (zs == 0 && (xs & zs) == 0) { 
            // Write to global output memory. Use atomic add to avoid sync issues.
            atomicAdd(&d_chunk_tally[num_D*2], (tally_t)1);
            atomicAdd(&d_chunk_tally[num_D*2+1], (tally_t)sign);
        }
        
    }
}
