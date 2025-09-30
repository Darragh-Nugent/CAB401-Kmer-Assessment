#include <iostream>

__global__ void add(int *a, int *b, int *c, int N) {
    int i = threadIdx.x;
    if (i < N) c[i] = a[i] + b[i];
}

int main() {
    const int N = 5;
    int a[N] = {1, 2, 3, 4, 5};
    int b[N] = {10, 20, 30, 40, 50};
    int c[N];

    int *d_a, *d_b, *d_c;
    cudaMalloc((void**)&d_a, N * sizeof(int));
    cudaMalloc((void**)&d_b, N * sizeof(int));
    cudaMalloc((void**)&d_c, N * sizeof(int));

    cudaMemcpy(d_a, a, N * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, N * sizeof(int), cudaMemcpyHostToDevice);

    add<<<1, N>>>(d_a, d_b, d_c, N);

    cudaMemcpy(c, d_c, N * sizeof(int), cudaMemcpyDeviceToHost);

    for (int i = 0; i < N; i++)
        std::cout << c[i] << " ";
    std::cout << std::endl;

    cudaFree(d_a); cudaFree(d_b); cudaFree(d_c);
    return 0;
}
