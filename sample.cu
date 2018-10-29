#include <stdio.h>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <curand.h>
#include <curand_kernel.h>
#include <inttypes.h>
#define PI 3.14159265358979323846
using namespace std;

//Defined as powers of 2
#define samplesPerThread (long int)1  // Number of samples generated per thread.
#define threadsPerBlock (long int)1   // Number of threads per block.
#define blocksPerChunk (long int)1    // Number of blocks per output array.
#define numChunks (long int) 1        // Do the whole thing each time for a new gamma
#define samplesPerChunk samplesPerThread + threadsPerBlock + blocksPerChunk
#define nsamples numChunks + samplesPerChunk

// Device Memory
#define num_consts 3
#define tally_t int 
__device__ __constant__ uint64_t *d_eqn_masks;
__device__ __constant__ int d_num_eqns;
__device__ __constant__ double d_consts[num_consts];
// Chunktally stores the number of samples with n encountered D's in d_chunk_tally[2*t],
// and the tally taking into acount sign in d_chunk_tally[2*t+1]
__device__ tally_t *d_chunk_tally;
   
// Host Memory
// array of equations in bitmask form, i.e. x_2 + x_3 + x_4 for 5 variables is 01110
uint64_t *h_eqn_masks;
int num_eqns;

__global__ void sample(int seed);

// Count number of lines in file, which indicates number of equations
int count_lines(char *filename) {
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
    num_eqns = count_lines(filename);
    h_eqn_masks = (uint64_t *) malloc(num_eqns*sizeof(uint64_t));
    
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
            if(counter < 3) b_eqn += pow(2,a);
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
    // first arugment is equation file, second is gamma
    if(argc < 3) {
        cout << "not enough arguments, please specify <equation file> and <gamma>" << endl;
        return 0;
    }

    double gamma = strtod(argv[2],NULL);
    read_file(argv[1]);

    // Copy bit mask array to device, to be passed d pointer later
    uint64_t *d_ptr;
    cudaMalloc((void **)&d_ptr, num_eqns*sizeof(uint64_t));
    cudaMemcpyToSymbol(d_eqn_masks, &d_ptr, sizeof(uint64_t *));
    cudaMemcpy(d_ptr, h_eqn_masks, num_eqns*sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_num_eqns, &num_eqns, sizeof(int));

    // Copy relevant D(e^{i\gamma C}) constants to device
    double tot = abs(sin(gamma)) / (abs(cos(gamma)) + abs(sin(gamma)));
    double sign_s, sign_c = 1;
    if(sin(gamma) < 0) sign_s = -1;
    if(cos(gamma) < 0) sign_c = -1;
    double h_consts[num_consts] = {tot, sign_s, sign_c};
    cudaMemcpyToSymbol(d_consts, h_consts, num_consts*sizeof(double));

    // We don't need the masks on the host.
    free(h_eqn_masks);

    // Host memory for tallying output.
    int tally_size = 2*num_eqns;
    tally_t* h_chunk_tally = (tally_t*) malloc(tally_size*sizeof(tally_t));
    tally_t* output_tally = (tally_t*)malloc(tally_size*sizeof(tally_t));

    // Initialize both arrays to 0
    memset(h_chunk_tally, 0, tally_size*sizeof(tally_t));
    memset(output_tally, 0, tally_size*sizeof(tally_t));

    tally_t *d_chunk_ptr;
    cudaGetSymbolAddress((void **)&d_chunk_ptr, d_chunk_tally);

    //for (int j = 0; j < (1 << numChunks); j++) {
    for(int j = 0; j < numChunks; j++) {
        //std::cout << "Running chunk " << (j+1) << " of " << (1 << numChunks) << std::endl;
        std::cout << "Running chunk " << (j+1) << " of " << numChunks << std::endl;
        // Take samples
        //sample<<<(1 << blocksPerChunk), (1 << threadsPerBlock)>>>(time(0)); //random version
        sample<<<1, threadsPerBlock>>>(time(0)); //random version
        // Wait for GPU to finish before accessing on host
        cudaDeviceSynchronize();
        // Copy samples to host, zero out device data
        cudaMemcpyFromSymbol(h_chunk_tally, d_chunk_tally, tally_size*sizeof(tally_t));
        //cudaMemset(d_chunk_tally, 0, tally_size*sizeof(tally_t)); // are we allowed to do this?
        cudaMemset(d_chunk_ptr, 0, tally_size*sizeof(tally_t));

        // Add chunk tally to overall tally
        for (int i = 0; i < tally_size; i++) output_tally[i] += h_chunk_tally[i];
    }

    // print output
    std::cout << nsamples << std::endl;
    for (int i = 0; 2*i < tally_size; i+=1) {
        std::cout << i << "," << output_tally[2*i] << "," << output_tally[2*i+1] << std::endl;
    }
    
    // Free memory
    free(h_chunk_tally);
    free(output_tally);
    return 0;
}

// Print an integer in binary
__device__ void printb(size_t const size, void const * const ptr) {
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    for(int i = size - 1; i >= 0; i--) {
        for(int j = 7; j >= 0; j--) {
            byte = (b[i] >> j) & 1;
            printf("%u", byte);
        }
    }
    printf("");
}

__device__ void print_xs_zs(uint64_t xs, uint64_t zs) {
    printf("xs: ");
    printb(sizeof(xs), &xs);
    printf("\n");
    printf("zs: ");
    printb(sizeof(zs), &zs);
    printf("\n");
}

// Flip the sign of x if data has an odd number of 1's in it
__device__ void parity(int* x, uint64_t data) {
    while (data) {
        *x *= -1;
        data = data & (data - 1);
    }
}

// Get a uniformly random integer inclusively between min and max
__device__ int get_rand_int(curandState_t state, int min, int max) {
    float rand_f = curand_uniform(&state);
    rand_f *= (max - min + 0.999999);
    rand_f += min;
    return (int)truncf(rand_f);
}

// apply the iZZZ \rho to our state xs and zs
__device__ void apply_iZZZ(uint64_t *xs, uint64_t *zs, uint64_t mask) {
    *zs = 0;
}

__global__ void sample(int seed) {
    // Initialize curand
    curandState_t state;
    curand_init(seed, blockIdx.x, threadIdx.x, &state);

    printf("STARTING IN DEVICE CODE NOW\n");
    printf("block, thread %d, %d\n", blockIdx.x, threadIdx.x);
    printf("num eqns: %d\n", d_num_eqns);
    for(int i = 0; i < 3; i++) {
        printf("d_consts[%d]: %f\n", i, d_consts[i]);
    }
    /**
    for(int i = 0; i < d_num_eqns; i++) {
        printf("eqn %d: %" PRIu64 "\n", i, d[i]);
    }
    **/
    
    // Per thread local memory. Can probably make this smaller with uglier code.
    uint64_t xs, zs;
    tally_t num_D = 0; 
    int sign = 1;
    
    //for(int j = 0; j < (1 << samplesPerThread); j++) {
    for(int j = 0; j < samplesPerThread; j++) {
        // Pick a random equation from eqn_masks
        int rand = get_rand_int(state, 0, d_num_eqns - 1);
        printf("rand: %d\n", rand);
        printf("--------- INIT ---------\n");
        uint64_t init_mask = d_eqn_masks[rand];
        xs = init_mask;
        zs = init_mask;

        print_xs_zs(xs, zs);
        printf("-------- Applying e^{i gamma C} --------\n"); 
        for(int i = 0; i < d_num_eqns; i++) {
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
                if(rand_f <= d_consts[0]) {
                    printf("in apply branch\n");
                    apply_iZZZ(&xs, &zs, mask);
                    printf("-------------------- aaa ----------------------\n");
                    print_xs_zs(xs, zs);
                    if(d_consts[1] < 0) sign *= -1;
                } else if(d_consts[2] < 0) sign *= -1;
                num_D += 1;            
            }
        }
        
        // Because <+|Y|+> = <+|Z|+> = 0, we only care if both of these don't happen
        if (zs == 0 && (xs & zs) == 0) { 
            printf("~~~~~~~~~~~ Doing something to tally ~~~~~~~~~~~~~~~, num_D: %d\n", num_D);
            // Write to global output memory. Use atomic add to avoid sync issues.
            atomicAdd(&d_chunk_tally[num_D*2], (tally_t) 1);
            atomicAdd(&d_chunk_tally[num_D*2+1], (tally_t) sign);
        }
    }
}
