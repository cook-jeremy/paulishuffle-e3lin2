#include <stdio.h>
#include <fstream>
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
using namespace std;

// Device Memory
//__device__ int *d;
   
//int *h;
int num = 3;

__global__
void sample(int *d, int num) {
    printf("ON DEVICE\n");
    for(int i = 0; i < num; i++) {
        printf("device: %d\n", d[i]);
    }
}

int main(int argc, char **argv) {
    //int *h = (int *) malloc(num*sizeof(int));
    int h[num];
    h[0] = 2;
    h[1] = 3;
    h[2] = 5;

    for(int i = 0; i < num; i++) {
        printf("host: %d\n", h[i]);
    }

    int *d;   
    cudaMalloc((void **)&d, num*sizeof(int));
    cudaMemcpy(d, h, num*sizeof(int), cudaMemcpyHostToDevice);

    sample<<<1, 1>>>(d, num); 
    cudaDeviceSynchronize();
    printf("done\n");
    cudaFree(d);
    return 0;
}
