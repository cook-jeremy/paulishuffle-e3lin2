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
#include <thrust/complex.h>
#define PI 3.14159265358979323846
#define E  2.71828182845904523536
using namespace std;

//Defined as powers of 2 
#define samplesPerThread (long int) 5 // Number of samples generated per thread.
#define threadsPerBlock (long int) 5   // Number of threads per block.
#define blocksPerChunk (long int) 5   // Number of blocks per output array.
#define numChunks (long int) 5         // Do the whole thing each time for a new gamma
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

void read_file(char* filename) {
    h_eqn_masks = (uint64_t *) malloc(num_eqns*sizeof(uint64_t));
    h_sols = (bool *) malloc(num_eqns*sizeof(bool));
        
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
                b_eqn += pow(2, a);
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
    fclose(fp);
}

int main(int argc, char **argv) {
    if(argc != 5) {
        cout << "not enough arguments, please specify\n <equation file> <num_eqns> <num_vars> <gamma>" << endl;
        return 0;
    }

    // Read off arguments and eqns masks from eqn file
    num_eqns = strtod(argv[2], NULL);
    num_vars = strtod(argv[3], NULL);
    read_file(argv[1]);
    double gamma = strtod(argv[4], NULL);

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
    // Copy num equations, vars, and gamma to device
    cudaMemcpyToSymbol(d_num_eqns, &num_eqns, sizeof(int));
    cudaMemcpyToSymbol(d_num_vars, &num_vars, sizeof(int));
    cudaMemcpyToSymbol(d_gamma, &gamma, sizeof(double));
    
    // We don't need the masks or sols on the host.
    free(h_eqn_masks);
    free(h_sols);

    double estimate = 0;

    for(int i = 0; i < num_eqns; i++) {
        // copy over eqn #
        cudaMemcpyToSymbol(d_picked_eqn, &i, sizeof(int));

        double h_total = 0;
        double *total_ptr;
        cudaGetSymbolAddress((void **)&total_ptr, d_total);

        for (int j = 0; j < (1 << numChunks); j++) {
            // Get current nanosecond for seed
            struct timespec now;
            clock_gettime(CLOCK_MONOTONIC_RAW, &now);
            sample<<<(1 << blocksPerChunk), (1 << threadsPerBlock)>>>((double) now.tv_nsec); //random version

            // Wait for GPU to finish before accessing on host
            cudaDeviceSynchronize();
            
            h_total = 0;
            // Copy total to host
            cudaMemcpy(&h_total, total_ptr, sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemset(total_ptr, 0, sizeof(double));
            estimate += h_total;
        }
    }
    
    estimate /= 2*(1 << nsamples);
    // Print estimate and error
    printf("%lf\n", estimate);
    double delta = 0.05;
    double error = num_eqns*sqrt(4*log(2/delta)/ (1 << nsamples))/2;
    printf("%lf\n", error);
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
    thrust::complex<float> total_phase = 1.0;
    // for each equation, calculate overlap
    for(int i = 0; i < d_num_eqns; i++) { 
        float arg1 = pow(-1.0f, (float) (d_sols[i] + num_ones(x & d_eqn_masks[i])));
        thrust::complex<float> arg2 = -thrust::complex<float>(0.0f, 1.0f)*d_gamma*0.5*arg1;
        total_phase *= thrust::exp(arg2);
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
        x = get_rand_int(&state, 0, pow(2.0f, (float) d_num_vars) - 1);
        alpha = -thrust::complex<float>(0.0f, 1.0f)*pow(-1.0f, (float) (d_sols[d_picked_eqn] + num_ones(x & y)));
        //alpha = -thrust::complex<float>(0.0f, 1.0f)*pow(-1.0f, (float) (num_ones(x & y)));
        inner1 = thrust::conj(inner_x_w(x ^ y));
        inner2 = inner_x_w(x);
        final = inner1*inner2*alpha;
        atomicAdd(&d_total, final.real());
    }
}
