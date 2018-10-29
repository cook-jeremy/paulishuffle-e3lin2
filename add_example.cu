#include <stdio.h>
#define N 3

__device__ int *d;

__global__ void add() {
    printf("IN DEVICE CODE\n");
    for(int i = 0; i < 3; i++) {
        printf("d[%d]: %d\n", i, d[i]);
    }
}

int main() {
    int *h = (int *) malloc(N*sizeof(int));
    for (int i = 0; i<N; ++i) {
        h[i] = i+1;
        printf("h[%d]: %d\n", i, h[i]);
    }

    int *d_ptr;
    //cudaGetSymbolAddress((void **)&d_ptr, d);
    cudaMalloc((void **) &d_ptr, N*sizeof(int));

    cudaMemcpyToSymbol(d, &d_ptr, sizeof(int *));
    cudaMemcpy(d_ptr, h, N*sizeof(int), cudaMemcpyHostToDevice);
    free(h);

    add<<<1, 1>>>();
    cudaDeviceSynchronize();

    printf("returned\n");
    //cudaFree(d_ptr);
    return 0;
}
