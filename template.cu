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

//Defined as powers of 2
#define samplesPerThread (long int)3 // Number of samples generated per thread.
#define threadsPerBlock (long int)10  // Number of threads per block.
#define blocksPerChunk (long int)10   // Number of blocks per output array.
#define numChunks (long int)2        // Do the whole thing this many times.

//Also powers of 2
#define samplesPerChunk samplesPerThread + threadsPerBlock + blocksPerChunk
#define nsamples numChunks + samplesPerChunk


__global__ void sample(int seed);

// Device Memory
__device__ __constant__ Gate c_circ[depth*gatesPerLayer+1];
__device__ tallyType d_chunkTally[tallyArraySize]; // Global memory
// Chunktally stores the number of samples with t encountered t gates
// in d_chunkTally[2*t], and the tally taking into acount sign in d_chunkTally[2*t+1]

// question: can we gain performance by giving each block an output array?
// that way the output array can be in shared memory.

int print_device_properties(void);

int main(void) {
    //if (!print_device_properties()) return 0;
        
    // Read circuit. This should be in constant memory.
    Gate* h_circ = (Gate*)malloc(circuitSize); 
    FILE *circFile = fopen("circuit.circ", "rb");
    fread(h_circ, sizeof(char), circuitSize, circFile);
    fclose(circFile);
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
    
    for (int j = 0; j < (1 << samplesPerThread); j++) {

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

        // Compare to output string to get probability
		if (xs == 0) { // <0|X|0> = <1|X|1> = <0|Y|0> = <1|Y|1> = 0
            parity(&sign, zs & c_circ[depth*gatesPerLayer].data); // Only <1|Z|1> = -1

            // Write to global output memory. Use atomic add to avoid sync issues.
            atomicAdd(&d_chunkTally[out*2], (tallyType)1);
            atomicAdd(&d_chunkTally[out*2+1], (tallyType)sign);
        }
    }
}


int print_device_properties(void) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "Device constant mem:\t" << prop.totalConstMem << " bytes" << std::endl;
    std::cout << "Device global mem:\t" << prop.totalGlobalMem << " bytes" << std::endl;
    std::cout << "Shared mem/block:\t" << prop.sharedMemPerBlock << " bytes" << std::endl;
    std::cout << std::endl;

    std::cout << "Tally array size:\t" << tallyArraySize*sizeof(tallyType) << " bytes" << std::endl;
    std::cout << "Memory for circuit:\t" << (depth*gatesPerLayer+1)*sizeof(Gate) << " bytes" << std::endl;
    std::cout << std::endl;

    // Print sizes
    std::cout << "Samples per thread:\t2^" << samplesPerThread << " samples" << std::endl;
    std::cout << "Threads per block:\t2^" << threadsPerBlock << " samples" << std::endl;
    std::cout << "Blocks per chunk:\t2^" << blocksPerChunk << " samples" << std::endl;
    std::cout << "Total num samples:\t2^" << static_cast<long int>(nsamples) << " samples" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Threads per block:\t2^" << threadsPerBlock << " of ";
    std::cout << prop.maxThreadsPerBlock  << " threads" << std::endl;
    if ((1 << numChunks) > prop.maxThreadsPerBlock) {
        std::cout << "Too many threads per block!" << std::endl;
        return 0;
    }
    std::cout << "Blocks per chunk:\t2^" << blocksPerChunk << " blocks" << std::endl;
    std::cout << std::endl;
    return 1;
}
