// must compile with flag: -arch=sm_60
#include <stdio.h>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h> 
#include <curand.h>
#include <curand_kernel.h>
#include <inttypes.h>
//#include <complex.h>
#include <thrust/complex.h>
#define PI 3.14159265358979323846
#define E  2.71828182845904523536
using namespace std;

//Defined as powers of 2 
#define samplesPerThread (long int) 10 // Number of samples generated per thread.
#define threadsPerBlock (long int) 10   // Number of threads per block.
#define blocksPerChunk (long int) 5    // Number of blocks per output array.
#define numChunks (long int) 6         // Do the whole thing each time for a new gamma
#define samplesPerChunk samplesPerThread + threadsPerBlock + blocksPerChunk
#define nsamples numChunks + samplesPerChunk

// Device Memory
__device__ __constant__ uint64_t *d_eqn_masks;
__device__ __constant__ bool *d_sols;
__device__ __constant__ int d_num_eqns;
__device__ __constant__ int d_num_vars;
__device__ __constant__ int d_picked_eqn;
__device__ __constant__ double d_gamma;
__device__ double d_total = 0;
   
// Host Memory
uint64_t *h_eqn_masks; // array of equations in bitmask form, i.e. x_2 + x_3 + x_4 for 5 variables is 01110
bool *h_sols; // solutions to each equation, either 0 or 1
int num_eqns;
int num_vars;

__global__ void sample(int seed);

// Count number of lines in file, which indicates number of equations
int count_lines(char *filename) {
    FILE *fp = fopen(filename,"r");
    int ch = 0;
    int lines = 0;
    if (fp == NULL) return 0;
    while(!feof(fp)) {
        ch = fgetc(fp);
        if(ch == '\n') lines++;
    }
    fclose(fp);
    return lines;
}

void count_vars(int *vars, int var) {
    bool flag = false;
    for(int i = 0; i < num_vars; i++) {
        if(vars[i] == var) flag = true;
    }
    // add var to list
    if(!flag) {
        vars[num_vars] = var;
        num_vars++;
    }
}

void read_file(char* filename) {
    num_eqns = count_lines(filename);
    h_eqn_masks = (uint64_t *) malloc(num_eqns*sizeof(uint64_t));
    h_sols = (bool *) malloc(num_eqns*sizeof(bool));
    int *vars = (int *) malloc(3*num_eqns*sizeof(int));

    // initialize to all -1
    for(int i = 0; i < 3*num_eqns; i++) {
        vars[i] = -1;
    }
    
    // Create bitmasks    
    FILE *fp = fopen(filename, "r");
    for(int i = 0; i < num_eqns; i++) {
        char buff[255];
        fscanf(fp, "%s", buff);
        char *pt;
        pt = strtok(buff, ",");
        int counter = 0;
        uint64_t b_eqn = 0;
        while (pt != NULL) {
            int a = atoi(pt);
            if(counter < 3) {
                count_vars(vars, a); 
                b_eqn += pow(2,a);
            } else {
                h_sols[i] = a;
            }
            pt = strtok(NULL, ",");
            counter++;
        }
        // add to bitmask array
        h_eqn_masks[i] = b_eqn;
        b_eqn = 0;
    }
    free(vars);
    fclose(fp);
}

int main(int argc, char **argv) {
    // first arugment is equation file, second is gamma
    if(argc != 4) {
        cout << "not enough arguments, please specify <equation number> <equation file> <gamma>" << endl;
        return 0;
    }

    int picked_eqn = strtod(argv[1], NULL);
    double gamma = strtod(argv[3], NULL);
    read_file(argv[2]);

    // Copy bit mask array to device
    uint64_t *d_ptr;
    cudaMalloc((void **)&d_ptr, num_eqns*sizeof(uint64_t));
    cudaMemcpyToSymbol(d_eqn_masks, &d_ptr, sizeof(uint64_t *));
    cudaMemcpy(d_ptr, h_eqn_masks, num_eqns*sizeof(uint64_t), cudaMemcpyHostToDevice);
    // Copy solutions to equations to device
    bool *sol_ptr;
    cudaMalloc((void **)&sol_ptr, num_eqns*sizeof(bool));
    cudaMemcpyToSymbol(d_sols, &sol_ptr, sizeof(bool *));
    cudaMemcpy(sol_ptr, h_sols, num_eqns*sizeof(bool), cudaMemcpyHostToDevice);
    // Copy num equations and vars to device
    cudaMemcpyToSymbol(d_num_eqns, &num_eqns, sizeof(int));
    cudaMemcpyToSymbol(d_num_vars, &num_vars, sizeof(int));
    cudaMemcpyToSymbol(d_picked_eqn, &picked_eqn, sizeof(int));
    cudaMemcpyToSymbol(d_gamma, &gamma, sizeof(double));
    
    // We don't need the masks or sols on the host.
    free(h_eqn_masks);
    free(h_sols);

    double h_total = 0;
    double total = 0;
    double *total_ptr;
    cudaGetSymbolAddress((void **)&total_ptr, d_total);

    for (int j = 0; j < (1 << numChunks); j++) {
        //std::cout << "Running chunk " << (j+1) << " of " << (1 << numChunks) << std::endl;
        
        // Get current nanosecond for seed
        struct timespec now;
        clock_gettime(CLOCK_MONOTONIC_RAW, &now);

        // Take samples
        sample<<<(1 << blocksPerChunk), (1 << threadsPerBlock)>>>((double) now.tv_nsec); //random version

        // Wait for GPU to finish before accessing on host
        cudaDeviceSynchronize();
        
        h_total = 0;
        // Copy total to host
        cudaMemcpy(&h_total, total_ptr, sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemset(total_ptr, 0, sizeof(double));

        h_total /= (1 << samplesPerThread) * (1 << threadsPerBlock) * (1 << blocksPerChunk);
        total += h_total;
    }

    total /= (1 << numChunks);

    //printf("%ld\n", nsamples);
    printf("%f\n", total);

    float delta = 0.001;
    float eps = sqrt(2*log(2/delta)/pow(2,nsamples));
    //float error = eps*num_eqns/2;
    printf("%f\n", eps);

    return 0;
}

// Get a uniformly random integer inclusively between min and max
__device__ int get_rand_int(curandState_t *state, int min, int max) {
    int out = curand(state) % (max-min + 1);
    return out + 0;
}

__device__ int num_ones(uint64_t n) {
    int count = 0;
    while(n) {
        count += n & 1;
        n >>= 1;
    }
    return count;
}


__device__ thrust::complex<float> inner_x_w(uint64_t x)  {
    thrust::complex<float> total_phase = 1/sqrt(pow(2.0f, d_num_vars));
    // for each equation, calculate overlap
    for(int i = 0; i < d_num_eqns; i++) { 
        double arg1 = pow(-1.0f, d_sols[i] + num_ones(x & d_eqn_masks[i]));
        thrust::complex<float> arg2 = -thrust::complex<float>(0.0f, 1.0f)*d_gamma*0.5*arg1;
        total_phase *= pow((float) E, arg2);
    }
    return total_phase;
}


__global__ void sample(int seed) {
    // Initialize curand
    curandState_t state;
    curand_init(seed, blockIdx.x, threadIdx.x, &state);

    uint64_t y = d_eqn_masks[d_picked_eqn];
    uint64_t x;
    thrust::complex<float> alpha, inner1, inner2, final;
    
    for(int j = 0; j < (1 << samplesPerThread); j++) {
        x = get_rand_int(&state, 0, pow(2.0f, d_num_vars) - 1);
        alpha = thrust::complex<float>(0.0f, -1.0f)*pow(-1.0f, num_ones(x & y));
        inner1 = thrust::conj(inner_x_w(x ^ y));
        inner2 = inner_x_w(x);
        final = pow(2.0f, d_num_vars)*inner1*inner2*alpha;
        atomicAdd(&d_total, final.real());
    }
}
