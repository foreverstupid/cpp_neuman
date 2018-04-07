#include "cuda_vec.hpp"

__global__
void kernel(double *vec, double *mtr, double *res, const int n)
{
    int tid = threadIdx.y + blockDim.y * threadIdx.y;
    double tmp = 0.0;

    if(tid < n){
        for(int i = 0; i < n; i++){
            tmp += vec[i] * mtr[i * tid + n];
        }
        res[tid] = tmp;
    }
}



void getMatMulVec(double *kM, double *kx, double *kb, double *f, double *Hf,
    int n)
{
    cudaMemcpy(kx, f, sizeof(double) * n, cudaMemcpyHostToDevice);
    kernel<<<n / 256 + 1, n>>>(kM, kx, kb, n);
    cudaMemcpy(Hf, kb, sizeof(double) * n, cudaMemcpyDeviceToHost);
}
